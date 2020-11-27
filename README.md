cd /home/annndrey/work/git/MolD_release
uwsgi --ini uwsgi.ini --socket :9300 --protocol=http
celery worker -E -A api.celery