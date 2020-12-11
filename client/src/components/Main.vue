<template>
<div class="hello">
  <div class="container">
    <div class="row" id="header">
      
      <div class="col-md-3 col-md-offset-3">
	<img src="@/assets/logo1.png" height="100px">
      </div>
      <div class="col-md-9">
	<h2 class="text-secondary">  {{ msg }}</h2>
	<a href="MolD_manual_Web.pdf" target="_blank" rel="noopener noreferrer">Download manual</a>
      </div>
    </div>

    <div class="row mt-3">
      <div class="col-md-9 mt-3">
	<h4 class="text-left">Data file</h4>
	<div class="form-group">
	  <div class="row mt-3">
	    <div class="col-md-1 align-self-start">
	      <span  class="align-top">
		<a role="button" data-toggle="collapse" data-target="#collapseExample" aria-expanded="false" aria-controls="collapseExample">
		  <span class="badge badge-info">
		    ?
		  </span>
		</a>
	      </span>
	    </div>
	    <div class="col-md-11 align-self-start">
	      <input type="file" ref="datafile" @change="readInputFile" class="form-control-file" id="dataUpload" aria-describedby="duHelp">
	    </div>
	  </div>
	  <div class="collapse" id="collapseExample">
	    <small id="ulHelp" class="form-text text-muted text-left">

	      The input file is in fasta format: each entry starts with the identifier line, and one or more lines of nucleotide sequence. Identifier line starts with ‘>’ and must contain two elements, separated by a pipe (‘|’) character. First is a free-style sequence identifier, without pipe symbol, the second is the taxon name of the query level. The names of the taxa to be diagnosed correspond to the second element.<br>
	      1. query ID example:
	      >GBXXXXXXX|query<br>
	      2. reference ID example:
	      >GBXXXXXXX|ref1<br>
	      	
	      <a href="Example_input_formatting-Puillandre_et_al_2014_Conidae_genera.fas" target="_blank" rel="noopener noreferrer">Example file</a>
	    </small>
	  </div>
	</div>
	  
      </div>

     
    </div>

    <div class="row mt-3">
      <div class="col-md-12">
	<flash-message ref="flashblock" transition-name="fade" class="flashmessage"></flash-message>
      </div>
    </div>
    
    <div class="row mt-3">
      <div class="col-md-12">
	<h4 class="text" style="color: #42b983">Input parameters</h4>
      </div>
    </div>
    
    <div class="row">
      <div class="col-md-8">
	<div class="form-group">
	  <label for="taxaList">List of query taxa
	    <!--<a role="button" data-toggle="collapse" data-target="#collapseExample4" aria-expanded="false" aria-controls="collapseExample4">-->
	      <span class="badge badge-info">
		?
	      </span>
	    </a>
	    <div class="collapse" id="collapseExample4">
	      The names of the taxa to be diagnosed correspond to the second column in the input file
	    </div>
	  </label>
	  <input type="test" v-model="taxalist" class="form-control" id="taxaList" aria-describedby="tlHelp">
	  <small id="tlHelp" class="form-text text-muted">comma separated list of taxa or ALL (if each taxon is to be diagnosed)</small>
	</div>
      </div>


      <div class="col-md-2">
	<div class="form-group">
	  <label for="numNuc">gaps as characters
	    <a role="button" data-toggle="collapse" data-target="#collapseExampleGaC" aria-expanded="false" aria-controls="collapseExampleGaC">
	      <span class="badge badge-info">
		?
	      </span>
	    </a>
	  </label>
	  <input type="test" v-model="gapsaschars" class="form-control" id="gapsAsChars" aria-describedby="gacHelp">
	  <small id="gacHelp" class="form-text text-muted">default 'no'</small>
	  <div class="collapse" id="collapseExampleGaC">
	    <small id="gacHelp" class="form-text text-muted">'yes' (if alignment gaps to be coded as a character), or 'no' if not.</small>
	  </div>
	</div>
      </div>

      
      <div class="col-md-2">
	<div class="form-group">
	  <label for="taxonRank">Taxon rank</label>
	  <input type="number" min="1" max="2" v-model="taxonrank" class="form-control" id="taxonRank" aria-describedby="trHelp">
	  <small id="trHelp" class="form-text text-muted">1 for species, 2 for above species</small>
	</div>
      </div>
    </div>

    <div class="row">
      <div class="col-md-12">
	<h4>
	  <a data-toggle="collapse" href="#collapseExample5" role="button" aria-expanded="false" aria-controls="collapseExample5">
	    Advanced parameters
	  </a>
	</h4>
      </div>
    </div>
    <div class="collapse" id="collapseExample5">
      <div class="row">
	<div class="col-md-12">
	  <h4 class="text-left">
	    mDNC recovery
	  </h4>
	</div>
      </div>
      <div class="row">
	<div class="col-md-2">
	  <div class="form-group">
	    <label for="cutOff">Cutoff
	      <a role="button" data-toggle="collapse" data-target="#collapseExample1" aria-expanded="false" aria-controls="collapseExample1">
	      <span class="badge badge-info">
		?
	      </span>
	      </a>
	    </label>	  
	    <input type="test" v-model="cutoff" class="form-control" id="cutOff" aria-describedby="cfHelp">
	    <small id="cfHelp" class="form-text text-muted">default 100</small>
	    <div class="collapse" id="collapseExample1">
	      <small id="cfelp" class="form-text text-muted">Number of informative positions for focus taxon to be considered, natural number</small>
	    </div>
	  </div>
	</div>
	<div class="col-md-2">
	  <div class="form-group">
	    <label for="numNuc">NNNNN...
	      <a role="button" data-toggle="collapse" data-target="#collapseExample2" aria-expanded="false" aria-controls="collapseExample2">
		<span class="badge badge-info">
		  ?
		</span>
	      </a>
	    </label>
	    <input type="test" v-model="numnucl" class="form-control" id="numNuc" aria-describedby="nnHelp">
	    <small id="nnHelp" class="form-text text-muted">default 25</small>
	    <div class="collapse" id="collapseExample2">
	      <small id="nnHelp" class="form-text text-muted">Allowed number of ambiguously called nucleotides per sequence, natural number</small>
	    </div>
	  </div>
	</div>
	
	<div class="col-md-2">
	  <div class="form-group">
	    <label for="numIter">Num search iterations</label>
	    <input type="test" v-model="numiter" class="form-control" id="numIter" aria-describedby="niHelp">
	    <small id="niHelp" class="form-text text-muted">default 10000</small>
	  </div>
	</div>
	<div class="col-md-3">
	  <div class="form-group">
	    <label for="mlRaw">Max length for the draft DNCs</label>
	    <input type="test" v-model="maxlenraw" class="form-control" id="mlRaw" aria-describedby="mlrHelp">
	    <small id="mlrHelp" class="form-text text-muted">default 12</small>
	  </div>
	</div>
	<div class="col-md-3">
	  <div class="form-group">
	    <label for="mlRef">Max length for the refined mDNCs</label>
	    <input type="test" v-model="maxlenrefined" class="form-control" id="mlRef" aria-describedby="mlrfHelp">
	    <small id="mlrfHelp" class="form-text text-muted">default 7</small>
	  </div>
	</div>
      </div>
      
      <div class="row">
	<div class="col-md-12">
	  <h4 class="text-left">rDNC scoring</h4>
	</div>
      </div>
      
      <div class="row">
	<div class="col-md-4">
	  <div class="form-group">
	    <label for="pDiff">Pdiff</label>
	    <input type="test" v-model="pdiff" class="form-control" id="pDiff" aria-describedby="pdHelp">
	    <small id="pdHelp" class="form-text text-muted">Percent difference between original and modified sequence (default 1 for species-level taxa, 3 for for supraspecific taxa)</small>
	  </div>
	</div>
	<div class="col-md-4">
	  <div class="form-group">
	    <label for="nMax">NMaxSeq</label>
	    <input type="test" v-model="nmax" class="form-control" id="nMax" aria-describedby="nmHelp">
	    <small id="nmHelp" class="form-text text-muted">Max number of sequences per taxon to modify (default 10)</small>
	  </div>
	</div>
	<div class="col-md-4">
	  <div class="form-group">
	    <label for="threshRating">Threshold of rDNC scoring
	      <a role="button" data-toggle="collapse" data-target="#collapseExample3" aria-expanded="false" aria-controls="collapseExample3">
		<span class="badge badge-info">
		  ?
		</span>
	      </a>
	    </label>
	    <input type="test" v-model="thresh" class="form-control" id="threshRating" aria-describedby="thrsHelp">
	    <small id="thrHelp" class="form-text text-muted">default 75</small>
	    <div class="collapse" id="collapseExample3">
	      <small id="thrsHelp" class="form-text text-muted">
		To score each candidate rDNC, 100 simulated test datasets are created. If rDNC remains valid in a test dataset, it adds 1 to the score, so lowest possible score is 0 and highest is 100. If two consecutive scores are above the threshold value (defined here) the rDNC is output. The threshold values are like:
		<p class="text-left">lousy: 66</p>
		<p class="text-left">moderate: 75</p>
		<p class="text-left">stringent: 90</p>
		<p class="text-left">very_stringent: 95</p>
		<p class="text-left">default: moderate</p>
	      </small>
	    </div>
	  </div>
	</div>
      </div>
    </div>
    
    <div class="row">
      <div class="col-md-12">
	<h4 class="text" style="color: #9c8bb8">Your email</h4>
	<div class="form-group">
	  <!--<label for="userEmail">Your email</label>-->
	  <input type="email" v-model="email" class="form-control" id="userEmail" aria-describedby="emlHelp">
	  <small id="emlHelp" class="form-text text-muted">We will send the results to the provided email</small>
	</div>
      </div>
    </div>
    
    
    <div class="row mt-5">
      <div class="col-md-12">
	<button class="btn btn-primary" v-on:click="submitFile()">Submit</button>
      </div>
    </div>
    <div class="row" v-if="imagesource" >
      <div class="col-md-6 offset-3 mt-3">
	<figure class="figure ">
	  <img :src="imagesource.url" class="figure-img img-fluid rounded" alt="...">
	  <figcaption class="figure-caption">
	    {{imagesource.source}} · {{imagesource.width}}x{{imagesource.height}}
	  </figcaption>
	</figure>
      </div>
    </div>
    <p v-if="loading"><img class="loading" src="@/assets/loading.png" height="40px"></p>
  </div>
</div>
</template>

<script>
export default {
    name: 'Main',
    props: {
	msg: String
    },
    data () {
	return {
	    datafile: '',
	    gapsaschars: 'no',
	    taxalist: 'ALL',
	    taxonrank: 1,
	    cutoff: 100,
	    numnucl: 25,
	    numiter: 10000,
	    maxlenraw: 12,
	    maxlenrefined: 7,
	    pdiff: 1,
	    nmax: 10,
	    thresh: 75,	    
            email: '',
	    loading: false,
	    filetaxa: null
	}
    },
    watch: {
	taxonrank: function(newValue) {
	    if (newValue < 1) {
		this.taxonrank = 1
	    } else if (newValue > 2) {
		this.taxonrank = 2
	    }
	    
	    if (this.taxonrank == 1) {
		this.pdiff = 1
	    } else if (this.taxonrank == 2) {
		this.pdiff = 3
	    }
	},
	//taxalist: function(newValue) {
	//    console.log("taxalist", newValue)
	//},
	//filetaxa: function(newValue) {
	//    console.log("filetaxa", newValue)
	//}

	
    },
    methods: {
	updateFileTaxa(lst) {
	    this.filetaxa = lst
	},
	readInputFile() {
	    let inputfile = this.$refs.datafile.files[0]
	    let reader = new FileReader()
	    reader.readAsText(inputfile, "UTF-8")
	    let validateSeq = this.validateSeq
	    let validateIdentifier = this.validateIdentifier
	    let fileform = this.$refs.datafile

	    var messages = document.getElementsByClassName("flash__message")
	    if (messages.length > 0) {
		Array.from(messages).forEach( m => { 
		    m.parentNode.removeChild(m)
		})
	    }
		
	    let flash = this.flash
	    let filetaxa = this.filetaxa
	    let updateFileTaxa = this.updateFileTaxa
	    reader.onload = function() {
		let res = reader.result.split("\n")
		let linenum = 1
		let fileErrors = false

		// validate same sequences length
		//
		filetaxa = []
		let lengths = []
		let fullseq = ""
		res.forEach( (r, ind) => {
		    if (!r.includes(">") & !r.includes("|")) {
			fullseq = fullseq+r.replace("\r","")
		    } else {
			if (fullseq.length > 0) {
			    filetaxa.push(r.replace("\r","").split("|")[1])
			    lengths.push(fullseq.length)
			    fullseq = ""
			}
		    }
		})
		if (!lengths.every( (val, i, arr) => val === arr[0] )) {
		    flash("Input alignment contains sequences of different lengths", 'error', {timeout: 0})
		    fileErrors = true
		}
		if (res[res.length-1] == ""){
		    res.pop()
		}
		    
		while (res.length > 0) {
		    let fileStr = res.splice(0, 1)[0]
		    let isValid = false
		    if ( fileStr.includes(">") | fileStr.includes("|") ) {
			console.log("Val ID")
			isValid = validateIdentifier(fileStr, linenum)
		    } else {
			isValid = validateSeq(fileStr, linenum)
		    }
		    linenum = linenum + 1
		    if (!isValid) {
			fileErrors = true
			break
		    }
		}
		if (!fileErrors) {
		    flash("Input data validation passed", 'success', {timeout: 5000})
		} else {
		    fileform.type = 'text'
		    fileform.type = 'file'
		}
		updateFileTaxa(filetaxa)		
	    }

	},
	validateSeq(seq, linenum) {
	    let message = "Input data validation passed"
	    let res = true
	    if (seq) {
		seq = seq.replace("\r", "")
	    }
	    if (seq.length < 2) {
		message = "Missing data in the line " + linenum
		res = false
	    }
	    if (seq.search(/[^ACDEFGHIKLMNPQRSTUVWY-\s]/i) != -1) {
		message = "Sequence in the line " + linenum + " contains forbidden characters"
		res = false
	    }
	    if (!res) {
		this.flash(message, 'error', {timeout: 0})
	    }
	    return res
	},
	validateIdentifier(seq, linenum) {
	    let message = "Input data validation passed"
	    let res = true
	    if (seq) {
		seq = seq.replace("\r", "")
	    }

	    if (!/^\>[0-9a-zA-Z-_\.]+[^\|][\|][^\|][0-9a-zA-Z-_\.]+\r?\n?$/.test(seq)) {
		message = "Wrong identifier in the line " + linenum
		res = false
	    } 
	    if (!res) {
		this.flash(message, 'error', {timeout: 0})
	    }
	    return res
	},
	submitFile () {
	    let taxalistValidated = true
	    if (this.taxalist != "ALL") {
		let formlist = this.taxalist.split(",")
		if (formlist.length > 0) {
		    formlist.forEach(t => {
			t = t.trim()
			if (!this.filetaxa.includes(t)) {
			    taxalistValidated = false
			    this.flash(t + " not in allowed taxa list", 'error', {timeout: 0})
			}
		    })
		}
	    }
	    if (taxalistValidated) {
		this.datafile = this.$refs.datafile.files[0]
		let formData = new FormData()
		formData.append('filename', this.datafile)
		formData.append('gapsaschars', this.gapsaschars)
		formData.append('taxalist', this.taxalist)
		formData.append('taxonrank', this.taxonrank)
		formData.append('cutoff', this.cutoff)
		formData.append('numnucl', this.numnucl)
		formData.append('numiter', this.numiter)
		formData.append('maxlenraw', this.maxlenraw)
		formData.append('maxlenrefined', this.maxlenrefined)
		formData.append('pdiff', this.pdiff)
		formData.append('nmax', this.nmax)
		formData.append('thresh', this.thresh)
		formData.append('email', this.email)
		this.loading = true
		let token = this.token
		const headers = { 'Content-Type': 'multipart/form-data',}
		
		this.$axios.post('https://mold.testapi.me:6060/api/v2/data', formData, {headers: headers})
                    .then(request => {// this.parseResponse(request)
			this.loading = false
			this.flash('Your data is loaded, the results would be sent to the email provided', 'success', {timeout: 5000})
                    })
                    .catch(request => {//this.failedResponse(request)
			this.loading = false
			this.flash('Something went wrong, we cannot process the data', 'error', {timeout: 0})
                    })
	    }
	}
    },
}

</script>

<!-- Add "scoped" attribute to limit CSS to this component only -->
<style >
  .fade-enter-active, .fade-leave-active {
  transition: opacity .5s;
  }
  .fade-enter, .fade-leave-to /* .fade-leave-active до версии 2.1.8 */ {
  opacity: 0;
  }

h4.text
  
h3 {
  margin: 40px 0 0;
}
ul {
  list-style-type: none;
  padding: 0;
}
li {
  display: inline-block;
  margin: 0 10px;
}
a {
  color: #42b983;
}

p {
    margin: 0px 0;
}
</style>
