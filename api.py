#!/usr/bin/env python
# -*- coding: utf-8 -*-

from flask import Flask, make_response, request, current_app, send_file, url_for
from flask_restful import Resource, Api, reqparse, abort
from flask.json import jsonify
from flasgger import Swagger
from flask_cors import CORS, cross_origin
from flask_restful.utils import cors
import logging
import base64
import uuid
import tempfile
import os
import html2text

# for emails
import smtplib
from email.message import EmailMessage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from celery import Celery
from celery.schedules import crontab

import json

#from MolD_sDNC import SortedDisplayDict
#from MolD_sDNC import Step1
#from MolD_sDNC import C_VP_PP
#from MolD_sDNC import random_position
#from MolD_sDNC import step_reduction_complist
#from MolD_sDNC import ConditionD
#from MolD_sDNC import RemoveRedundantPositions
#from MolD_sDNC import Diagnostic_combinations
#from MolD_sDNC import IndependentKey
#from MolD_sDNC import random_sequence2
#from MolD_sDNC import GenerateBarcode2
#from MolD_sDNC import Screwed_dataset31
from MolD_sDNCFASTA import mainprocessing


logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s',
                    level=logging.DEBUG, datefmt='%d.%m.%Y %I:%M:%S %p')

app = Flask(__name__)
app.config.from_envvar('APPSETTINGS')
API_VERSION = app.config.get('API_VERSION', 1)

cors = CORS(app, resources={r"/*": {"origins": "*"}}, methods=['GET', 'POST', 'PATCH', 'DELETE', 'HEAD', 'OPTIONS'])
api = Api(app, prefix=f"/api/v{API_VERSION}")

gunicorn_logger = logging.getLogger('gunicorn.error')
app.logger.handlers = gunicorn_logger.handlers
app.logger.setLevel(gunicorn_logger.level)

REDIS_HOST = app.config.get('REDIS_HOST', 'localhost')
REDIS_PORT = app.config.get('REDIS_PORT', 6379)
REDIS_DB = app.config.get('REDIS_DB', 0)

MAILUSER = app.config.get('MAILUSER')
MAILPASS = app.config.get('MAILPASS')

app.config['SWAGGER'] = {
    'uiversion': 3
}
swtemplate = {
    "info": {
        "title": "MoID API",
        "description": "MoID API is the online version of a tree independent algorithm to retrieve diagnostic nucleotide characters from monolocus datasets. BioRxiv. DOI: 10.1101/838151",
        "version": f"{API_VERSION}",
    },
    "schemes": [
        "https"
    ]
}

swagger_config = {
    "headers": [
    ],
    "specs": [
        {
            "endpoint": 'apidescr',
            "route": '/apidescr.json',
            "rule_filter": lambda rule: True,
            "model_filter": lambda tag: True,
        }
    ],
    "static_url_path": "/flasgger_static",
    "swagger_ui": True,
    "specs_route": "/docs/",
    'uiversion': 3
}

swagger = Swagger(app, config=swagger_config, template=swtemplate)

SEND_EMAILS = True

app.config['CELERY_BROKER_URL'] = 'redis://{}:{}/{}'.format(REDIS_HOST, REDIS_PORT, REDIS_DB)
app.config['CELERY_RESULT_BACKEND'] = 'redis://{}:{}/{}'.format(REDIS_HOST, REDIS_PORT, REDIS_DB)

def make_celery(app):
    celery = Celery(app.import_name, broker=app.config['CELERY_BROKER_URL'])
    celery.conf.update(app.config)
    TaskBase = celery.Task
    class ContextTask(TaskBase):
        abstract = True
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)
    celery.Task = ContextTask
    return celery

celery = make_celery(app)

#@celery.on_after_configure.connect
#def setup_periodic_tasks(sender, **kwargs):
#    sender.add_periodic_task(
#        crontab(minute=0, hour='*/3'),
#        check_pending_notifications.s(),
#    )


@celery.task
def send_email_notification(email, status, parameters, taxalist, orig_filename):
    print("Sending email")
    sender = MAILUSER
    taxa = ", ".join(taxalist)
    print(taxalist)
    taxa1 = "-".join([t.replace(" ", "_") for t in taxalist])
    msg = MIMEMultipart('related')
    msg['Subject'] = f'MolD results: {taxa}'
    msg['From'] = sender
    msg['To'] = email

    # TODO: add three txt files:
    # 
    
    email_body = """\
<html>
    <head></head>
       <body>

    Results:
     {}

    ***************

    <p>Citation: <b>Fedosov A.E.</b>, Achaz G., Puillandre N. 2019. Revisiting use of DNA characters in taxonomy with MolD – a tree independent algorithm to retrieve diagnostic nucleotide characters from monolocus datasets. BioRxiv, published online on 11.11.2019. DOI: 10.1101/838151 </p>
       </body>
</html>
    """.format(status)
    
    plain_text_results = html2text.html2text(status)
    
    message_text = MIMEText(email_body, 'html')
    msg.attach(message_text)

    text_attachment = f"""
Results: 
{plain_text_results}

***************

Citation: Fedosov A.E., Achaz G., Puillandre N. 2019. Revisiting use of DNA characters in taxonomy with MolD – a tree independent algorithm to retrieve diagnostic nucleotide characters from monolocus datasets. BioRxiv, published online on 11.11.2019. DOI: 10.1101/838151
    """
    
    mail_attachment = MIMEText(text_attachment)
    fname = f"{orig_filename}.results.txt"
    mail_attachment.add_header('Content-Disposition', 'attachment', filename=fname)
    msg.attach(mail_attachment)
    
    print("mail ready to be sent")
    s = smtplib.SMTP('smtp.gmail.com', 587)
    s.ehlo()
    s.starttls()
    s.login(MAILUSER, MAILPASS)
    s.sendmail(sender, email, msg.as_string())
    s.quit()
    print("mail sent")


@celery.task
def process_data(email, gapsaschars, taxalist, taxonrank, cutoff, numnucl, numiter, maxlenraw, maxlenrefined, pdiff, nmax, thresh, tmpfname, orig_fname):
    print("Processing data")

    results, qclades = mainprocessing(gapsaschars, taxalist, taxonrank, cutoff, numnucl, numiter, maxlenraw, maxlenrefined, pdiff, nmax, thresh, tmpfname, orig_fname)

    parameters = f"""
    List of focus taxa: {taxalist}
    Taxon rank: {taxonrank}
    Cutoff: {cutoff}
    NNNNN...: {numnucl}
    Num iterations: {numiter}
    Max length for the raw mDNCs: {maxlenraw}
    Max length for the refined mDNCs: {maxlenrefined}
    Pdiff: {pdiff}
    NMaxSeq: {nmax}
    Threshold of rDNC rating: {thresh}
    """
    
    send_email_notification.delay(email, results, parameters, qclades, orig_fname)
    
@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify({'error': 'Not found'}), 404)


class DataAPI(Resource):
    def __init__(self):
        self.method_decorators = []

    def options(self, *args, **kwargs):
        return jsonify([])


    #@token_required
    @cross_origin()
    def post(self):
        """
        POST data to analyse
        ---
        parameters:
         - in: formData
           name: file
           type: file
           required: true
           description: Data file
         - in: formData
           name: tags
           type: array
           required: false
           description: A list of image tags to link with the image
           items:
             type: string
             description: Tag text
         - in: formData
           name: width
           type: integer
           required: false
           description: Image width
         - in: formData
           name: height
           type: integer
           required: false
           description: Image height

        responses:
          200:
            description: Data processing
          400:
            description: Bad request
        """
        app.logger.debug(request.files)
        app.logger.debug(request.values)
        email = request.values.get('email', None)
        gapsaschars = request.values.get('gapsaschars', "no")
        taxalist = request.values.get('taxalist', "ALL")
        taxalist = taxalist.replace(" ", '')
        taxonrank = int(request.values.get('taxonrank', 2))
        cutoff = request.values.get('cutoff', ">0")
        numnucl = int(request.values.get('numnucl', 5))
        numiter = int(request.values.get('numites', 10000))
        maxlenraw = int(request.values.get('maxlenraw', 12))
        maxlenrefined = int(request.values.get('maxlenrefined', 7))
        pdiff = int(request.values.get('pdiff', 1))
        #prseq = float(request.values.get('prseq', 0.1))
        nmax = int(request.values.get('nmax', 20))
        thresh = int(request.values.get('thresh', 85))
        
        if not all([email, gapsaschars, taxalist, taxonrank, cutoff, numnucl, numiter, maxlenraw, maxlenrefined, pdiff, nmax, thresh]):
            return make_response(jsonify({'error': 'Parameter missing'}), 400)
            
        if not request.files:
            return make_response(jsonify({'error': 'No file provided'}), 400)

        
        data_file = [request.files.get(f) for f in request.files][-1]

        thumb_height = request.values.get('height', None)
        tmpdname = "/tmp"
        tmpfname = str(uuid.uuid4())
        tmp_upl_file = os.path.join(tmpdname, tmpfname)
        data_file.save(tmp_upl_file)
        app.logger.debug(tmp_upl_file)
        process_data.delay(email, gapsaschars, taxalist, taxonrank, cutoff, numnucl, numiter, maxlenraw, maxlenrefined, pdiff, nmax, thresh, tmp_upl_file, data_file.filename)

        return jsonify("OK"), 200


api.add_resource(DataAPI, '/data', endpoint='processdata')

if __name__ == '__main__':
    app.debug = True
    app.run(host='0.0.0.0')
