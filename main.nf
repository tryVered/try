$HOSTNAME = ""
params.outdir = 'results'  

	// Add for each process an option to change the parameters. Default is the set params
//* params.nproc =  1  //* @input @description:"How many processes to use for each step. Default 1"
//* params.edit_collapse_seq_params =  "no"   //* @dropdown @options:"yes","no" @show_settings:"collapse_seq"
//* params.edit_split_seq_params =  "no"   //* @dropdown @options:"yes","no" @show_settings:"split_seq"
//* params.edit_parse_header_table =  "no"   //* @dropdown @options:"yes","no" @show_settings:"parse_headers" @description:"If edit true, then parametrs in the parse_header_table tab should be edited"
//* params.edit_parse_header1 =  "no"   //* @dropdown @options:"yes","no" @show_settings:"parse_headers" @description:"If edit true, then parametrs in the parse_header_table tab should be edited"
//* params.edit_parse_header2 =  "no"   //* @dropdown @options:"yes","no" @show_settings:"parse_headers" @description:"If edit true, then parametrs in the parse_header_table tab should be edited"

//* autofill
if ($HOSTNAME == "default"){
    $DOCKER_IMAGE = "immcantation/suite:4.3.0"
    $DOCKER_OPTIONS = "-v /work:/work"

}

//* platform
if ($HOSTNAME == "ig03.lnx.biu.ac.il"){
    $DOCKER_IMAGE = "immcantation/suite:4.3.0"
    $DOCKER_OPTIONS = "-v /work:/work"
	$CPU  = 48
    $MEMORY = 300 
}
//* platform


//* autofill



if((params.edit_collapse_seq_params && (params.edit_collapse_seq_params == "no"))){
    // Process Parameters for params.edit_collapse_seq_params:
    params.collapse_sequences_collapse_seq.act = "sum"
    params.collapse_sequences_collapse_seq.max_missing = 20
    params.collapse_sequences_collapse_seq.inner = "true"
    params.collapse_sequences_collapse_seq.uf = "CPRIMER"
    params.collapse_sequences_collapse_seq.cf = "CONSCOUNT"
    params.collapse_sequences_collapse_seq.nproc = params.nproc
}

if((params.edit_split_seq_params && (params.edit_split_seq_params == "no"))){
    // Process Parameters for params.Parse_header_parse_headers:
    params.split_sequences_split_seq.field = "CONSCOUNT"
    params.split_sequences_split_seq.num = 2
}

if((params.edit_parse_header_table && (params.edit_parse_header_table == "no"))){
    // Process Parameters for params.Parse_header_parse_headers:
    params.Parse_header_parse_headers.method = "table"
    params.Parse_header_parse_headers.act = ""
    params.Parse_header_parse_headers.args = "-f ID CREGION CONSCOUNT DUPCOUNT"
}

if((params.edit_parse_header1 && (params.edit_parse_header1 == "no"))){
    // Process Parameters for params.Parse_header1_parse_headers:
    params.Parse_header1_parse_headers.method = "table"
    params.Parse_header1_parse_headers.act = ""
    params.Parse_header1_parse_headers.args = "-f ID CREGION CONSCOUNT DUPCOUNT"
}

if((params.edit_parse_header2 && (params.edit_parse_header2 == "no"))){
    // Process Parameters for params.Parse_header2_parse_headers:
    params.Parse_header2_parse_headers.method = "table"
    params.Parse_header2_parse_headers.act = ""
    params.Parse_header2_parse_headers.args = "-f ID CREGION CONSCOUNT DUPCOUNT"
}


if (!params.reads){params.reads = ""} 

Channel.fromPath(params.reads, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_24_reads_g35_16;g_24_reads_g44_15}


process collapse_sequences_collapse_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_collapse-unique.fast.*$/) "reads_collapse_unique/$filename"}
input:
 set val(name), file(reads) from g_24_reads_g35_16

output:
 set val(name),  file("*_collapse-unique.fast*")  into g35_16_reads0_g36_20, g35_16_reads0_g43_15
 set val(name),  file("*_collapse-duplicate.fast*") optional true  into g35_16_reads_duplicate11
 set val(name),  file("*_collapse-undetermined.fast*") optional true  into g35_16_reads_undetermined22
 file "CS_*"  into g35_16_logFile33

script:
max_missing = params.collapse_sequences_collapse_seq.max_missing
inner = params.collapse_sequences_collapse_seq.inner
fasta = params.collapse_sequences_collapse_seq.fasta
act = params.collapse_sequences_collapse_seq.act
uf = params.collapse_sequences_collapse_seq.uf
cf = params.collapse_sequences_collapse_seq.cf
nproc = params.collapse_sequences_collapse_seq.nproc
failed = params.collapse_sequences_collapse_seq.failed

inner = (inner=="true") ? "--inner" : ""
fasta = (fasta=="true") ? "--fasta" : ""
act = (act=="none") ? "" : "--act ${act}"
cf = (cf=="") ? "" : "--cf ${cf}"
uf = (uf=="") ? "" : "--uf ${uf}"
failed = (failed=="false") ? "" : "--failed"

"""
CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act} --log CS_${name}.log ${failed}
"""

}


process Parse_header2_parse_headers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${out}$/) "outputparam/$filename"}
input:
 set val(name), file(reads) from g35_16_reads0_g43_15

output:
 set val(name),file("*${out}")  into g43_15_reads00

script:
method = params.Parse_header2_parse_headers.method
act = params.Parse_header2_parse_headers.act
args = params.Parse_header2_parse_headers.args

out="_reheader.fastq"
if(method=="collapse" || method=="add" || method=="copy" || method=="rename" || method=="merge"){
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} --act ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} -o ${name}.tab ${args}
			"""	
	}else{
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}



}


process split_sequences_split_seq {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_atleast-.*.fastq$/) "outputparam/$filename"}
input:
 set val(name),file(reads) from g35_16_reads0_g36_20

output:
 set val(name), file("*_atleast-*.fastq")  into g36_20_reads0_g37_15

script:
field = params.split_sequences_split_seq.field
num = params.split_sequences_split_seq.num

readArray = reads.toString()	
if(num!=0){
	num = " --num ${num}"
}else{
	num = ""
}

"""
SplitSeq.py group -s ${readArray} -f ${field} ${num}
"""

}


process Parse_header_parse_headers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${out}$/) "outputparam/$filename"}
input:
 set val(name), file(reads) from g36_20_reads0_g37_15

output:
 set val(name),file("*${out}")  into g37_15_reads00

script:
method = params.Parse_header_parse_headers.method
act = params.Parse_header_parse_headers.act
args = params.Parse_header_parse_headers.args

out="_reheader.fastq"
if(method=="collapse" || method=="add" || method=="copy" || method=="rename" || method=="merge"){
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} --act ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} -o ${name}.tab ${args}
			"""	
	}else{
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}



}


process Parse_header1_parse_headers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${out}$/) "outputparam/$filename"}
input:
 set val(name), file(reads) from g_24_reads_g44_15

output:
 set val(name),file("*${out}")  into g44_15_reads00

script:
method = params.Parse_header1_parse_headers.method
act = params.Parse_header1_parse_headers.act
args = params.Parse_header1_parse_headers.args

out="_reheader.fastq"
if(method=="collapse" || method=="add" || method=="copy" || method=="rename" || method=="merge"){
	"""
	ParseHeaders.py  ${method} -s ${reads} ${args} --act ${act}
	"""
}else{
	if(method=="table"){
			out=".tab"
			"""
			ParseHeaders.py ${method} -s ${reads} -o ${name}.tab ${args}
			"""	
	}else{
		"""
		ParseHeaders.py ${method} -s ${reads} ${args}
		"""		
	}
}



}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
