###################################################
######## shortpants configuration file ############
###################################################

###################################################
### Date: 2012-01-05                            ###
### Auther: Vikas Gupta, Stig Anderson          ###
###################################################

###################################################
### importing modules							###
import time,os,datetime							###
work_dir=os.getcwd()							###
miRNA_dir=work_dir+"/miRDP1.3"					###
###################################################


### define a hash for storing variables
var={}
f=open('configuration.txt','r')
for line in f:
	line = line.strip()
	if(len(line)>0):
		if(line[0] != '#'):
			token=line.split('=') ### take value after '='
			token[0]=token[0].replace(' ','')
			tokens=token[1].split('#') ### remove everything after #
			if(tokens[0]=="''"):
				var[token[0]]=''
			else:
				tokens[0]=tokens[0].replace("'",'')
				tokens[0]=tokens[0].replace('"','')
				tokens[0]=tokens[0].replace(' ','')
				tokens[0]=tokens[0].replace('\t','')
				var[token[0]]=tokens[0]

##########################################################
### please read and modify following wherever required ###
### only the lines not ending with "###" should be     ###
### modified                                           ###
##########################################################

##########################################################
### 			Dependencies						   ###
##########################################################
### please make sure following 						   ###
### 1. perl version v5.12.3 or higher installed        ###
### 2. python version 2.7.2 or higer installed         ###
### 3. matlibplot v1.0.1, scipy v0.9.0 and numpy 1.6.1 ###
###    or higer installed							   ###
### 4. bowtie version 0.12.7 or higer installed 	   ###
##########################################################

##########################################################
### define folder		###
if len(var['date'])==0:
	date=str(datetime.date.today()) 
else:
	date = var['date']
### if you want to use previous output folders then    ###	
### replace date varible by privous date               ###
### by dfault it takes current date and makes new      ###
### folder if unavailable                              ###
##########################################################

##########################################################
### define size range of small RNAs to considered for  ###
### creating profiles                                  ###
min_seq_size=int(var['min_seq_size'])
max_seq_size=int(var['max_seq_size'])
##########################################################


#### Mapping parameter
mismatches=var['mismatches']
max_maq_err=var['max_maq_err']
seedlen=var['seedlen']
no_of_cores=var['no_of_cores']

##########################################################
### read abundace cut-off for profiles, cut-off refers ### 
### to sum of read abundance is all libraries		   ###
cutoff=int(var['cutoff'])
##########################################################

##########################################################
### libraries genotypic headers	for each library	   ###		
### which can used in plots legends					   ###	
genotype=var['genotype'] 
genotype=[i for i in genotype.split(',')]
##########################################################


##########################################################
### Following is a list of Modules 					   ###
### to use (add "True" if using else "False")		   ###
### please remember that "True" and "False" are case   ###
### sensitive										   ###
##########################################################

##########################################################
################### Module 1 #############################
### Filtering and Sequencing Error correction          ###
##########################################################
### inputs are raw fastq files 						   ###
### here you should provided each library with all     ###
### sizes in seperate folder check example dataset for ###
### more details for demo							   ###
### for example: fq_file_dir=$HOME/libraries           ###
##########################################################

fq_file_dir=var['fq_file_dir']

#########################################################
######## adapter_filtering                     ##########
#########################################################
adapter_filtering=var['adapter_filtering']
r=adapter_filtering
if (r==True) or (r=='True') or r=='T' or r=='t' :
	adapter_filtering=True
adapter=var['adapter']
trim=var['trim']

keep_with_adapter_only=var['keep_with_adapter_only']
r=keep_with_adapter_only
if (r==True) or (r=='True') or r=='T' or r=='t' :
	keep_with_adapter_only=True



#########################################################
######## fastq size fractionation              ##########
#########################################################
size_fractionation=var['size_fractionation']
r=size_fractionation
if (r==True) or (r=='True') or r=='T' or r=='t' :
	size_fractionation=True



##########################################################
### please put folder name seperated by commas         ###
### for example: lib= "lib-1","lib-2"                  ###
##########################################################

lib=var['lib']
lib=[i for i in lib.split(',')]

infile=var['infile'] 

######################################################
### if True provide following parameters           ###
### for example server='dna-54-144.mb.au.dk' 	   ###	
### for example user= $HOME						   ###
### for example password=12345					   ###
### database must be named as 'profiles'		   ###
######################################################

server=var['server']
user=var['user']
password=var['password']
database=var['database']

genome_file=var['genome_file']

############### Use or not only mapped reads
use_only_mapped_reads=var['use_only_mapped_reads']
r=use_only_mapped_reads
if (r==True) or (r=='True') or r=='T' or r=='t' :
	use_only_mapped_reads=True
##########################################################
### True if module needed to be used                   ###
##########################################################
Raw_data_filtering=var['Raw_data_filtering']


######################################################
### if sequencing reads contains many unidentified ###
### based 'N', an Error correction is recommonded  ###
### please remember it is a fairly time consuming  ###
### step										   ###
### True if sub-module needed to be used           ###
######################################################
ErrorCorrection=var['ErrorCorrection']
r=ErrorCorrection
if (r==True) or (r=='True') or r=='T' or r=='t' :
	ErrorCorrection=True

Path_2_ECHO=var['Path_2_ECHO']
r=Path_2_ECHO
if (r==True) or (r=='True') or r=='T' or r=='t' :
	Path_2_ECHO=True



r=Raw_data_filtering
if (r==True) or (r=='True') or r=='T' or r=='t' :                    ###
	######################################################
	### if True provide following parameters           ###
	### provide file for repeat filtering              ###
	### for example repeat_database=$HOME/repeat.fa    ###
	######################################################
	
	repeat_database=var['filtering_file']
##########################################################
###################### Module 1.2 ########################
### Normalization and profiling                        ###
### True if module needed to be used                   ###
### submodule generates a normalized profile    	   ###
### table from the fastq files                         ###
### True if module needed to be used			       ###
##########################################################
Make_profile_files=var['Make_profile_files']
score_regulation_cutoff=int(var['score_regulation_cutoff'])
##########################################################
################### Module 1.3 ###########################
### submodule connects with the MySQL server    	   ###
### and stores data in tables with name same as 	   ###
### the profile output file                     	   ###
##########################################################

Add_output_to_mysql_database=var['Add_output_to_mysql_database']

r=Add_output_to_mysql_database
if (r==True) or (r=='True') or r=='T' or r=='t' :
	######################################################
	### if True provide following parameters           ###
	### for example server='dna-54-144.mb.au.dk' 	   ###	
	### for example user= $HOME						   ###
	### for example password=12345					   ###
	### database must be named as 'profiles'		   ###
	######################################################
	
	server=var['server']
	user=var['user']
	password=var['password']
	database=var['database']
	
	######################################################
	### make sure that normalized expression values    ###
	### file is in output foler						   ###	
	######################################################



##########################################################
###################### Module 2 ##########################
### Genomic Mapping                                    ###
### True if module needed to be used                   ###
##########################################################

Genome_mapping=var['Genome_mapping']
r=Genome_mapping
if (r==True) or (r=='True') or r=='T' or r=='t' :		###
	######################################################
	### if True provide following parameters           ###
	### for example genome_file=$HOME/reference.fa     ###
	######################################################
	
	
	######################################################
	### provide a input fasta file if it defers from   ###
	### the previous modules else leave it empty       ### 
	### for example infile=$HOME/reads.fa              ###
	######################################################
	
	infile=var['infile'] 
		

##########################################################
###################### Module 3 ##########################
### Clustering 				                           ###
### True if module needed to be used                   ###
##########################################################


##########################################################
###################### Module 3.1 ########################
### outputs small RNAs clusters 					   ###  
##########################################################

Cluster_prediction=var['Cluster_prediction']

r=Cluster_prediction
if (r==True) or (r=='True') or r=='T' or r=='t' :
	
	### for NiBLS
	
	NiBLS=var['NiBLS']
	r=NiBLS
	if (r==True) or (r=='True') or r=='T' or r=='t' :
		NiBLS=True
	
	######################################################
	### if True provide following parameters           ###
	### devide_by_size variable if True, invokes a size###
	### fractionated clustering algorithm else it 	   ###
	### considers all the sequences in the mapped .igv ###
	### file										   ###
	### abundance cut-off is the minimal read abundance###
	### in a gemonic region which doesn't contain a gap###
	### larger than 50 base-pair and has atleast 5     ###
	### unique sequences							   ###
	######################################################
	
	divide_by_size=var['divide_by_size']
	
	######################################################
	### abundance cut-off is the minimal read abundance###
	### in a gemonic region which doesn't contain a gap###
	### larger than 50 base-pair and has atleast 5     ###
	### unique sequences							   ###
	######################################################
	
	abundance_cutoff=int(var['abundance_cutoff'])
	max_gap_size=int(var['max_gap_size'])
	min_unique_reads=int(var['min_unique_reads'])
	
	######################################################
	### igv_infile is mapped file produced in module-2 ###
	### leave igv_infile='' 						   ###
	### if file used is same from module-2 			   ###
	### else provide a igv format mapped file		   ###
	######################################################
	
	igv_infile=var['igv_infile'] 
	
	######################################################
	### True if module needed to be used               ###
	### it will setup table for clusters in the MySQL  ###
	######################################################
	Add_cluster_to_mySQL=var['Add_cluster_to_mySQL']
	
	r=Add_cluster_to_mySQL
	if (r==True) or (r=='True') or r=='T' or r=='t' :	###
		Add_cluster_to_mySQL=True
		##################################################
		### if True provide following parameters       ###
		### for example server='dna-54-229.mb.au.dk'   ###	
		### for example user= $HOME					   ###
		### for example password=12345				   ###
		### database must be named as 'profiles'	   ###
		##################################################

		server=var['server']
		user=var['user']
		password=var['password']
		database=var['database']

'''
##########################################################
###################### Module 3.x ########################
### it will calculate amount unique sequences present  ###
### in form of clusters								   ###
### True if module needed to be used                   ###
##########################################################

Calculate_unique_sequences_in_cluster=var['Calculate_unique_sequences_in_cluster']

r=Calculate_unique_sequences_in_cluster
if (r==True) or (r=='True') or r=='T' or r=='t' :	   ###
	######################################################
	### if True provide following parameters       	   ###
	### use cluster_infile='' if cluster were produced ###
	### from module 4.1 else use the full path of 	   ###
	### cluster file								   ###
	######################################################
	
	cluster_infile=var['cluster_infile']


##########################################################
###################### Module 3.x ########################
### it calculates amount of differentailly regulated   ###
### sequences in the cluster filtered by score         ###
### where score is standard_deviation/square_root(avg) ###
### True if module needed to be used                   ###
##########################################################

Regulated_sequences_in_cluster=(var['Regulated_sequences_in_cluster'])

r=Regulated_sequences_in_cluster
if (r==True) or (r=='True') or r=='T' or r=='t' :			   ###
	######################################################
	### if True provide following parameters       	   ###
	### score cut-off should be higher to capture      ###
	### highly differentially regulated sequences      ###
	######################################################
	
	score_regulation_cutoff=int(var['score_regulation_cutoff'])
	
	######################################################
	### use cluster_infile='' if cluster were produced ###
	### from module 3.1 else use the full path of 	   ###
	### cluster file								   ###
	######################################################

	cluster_infile=var['cluster_infile']
	
	######################################################
	### regulation_file should be file containing      ###
	### differentially regulated sequences it should be###
	### left empty if we need to generate regulation   ###
	### based on regulation score cut-off              ###
	######################################################
	
	regulation_file=var['regulation_file'] 
'''
##########################################################
####################### Module 3.2 #######################
### this modules should ebe invoked if genomic region  ###
### for cluster is need to be annotated, it also       ###
### it also calculates number of exonic, intronic and  ###
### intergenic clusters								   ###
### True if module needed to be used                   ###
##########################################################

Cluster_genomic_region_prediction=var['Cluster_genomic_region_prediction']

r=Cluster_genomic_region_prediction
if (r==True) or (r=='True') or r=='T' or r=='t' :		   ###
	######################################################
	### if True provide following parameters       	   ###
	### gtf_file is a genemodel where exons introns are###
	### given with genomic co-ordinate check demo file ###
	### for more details                               ###
	######################################################
	
	gtf_file=var['gtf_file'] 
	
	######################################################
	### use cluster_infile='' if cluster were produced ###
	### from module 3.1 else use the full path of 	   ###
	### cluster file								   ###
	######################################################

	cluster_infile=var['cluster_infile']
	
##########################################################
###################### Module 4 ##########################
### True if module needed to be used                   ###
##########################################################


#########################################################
###################### Module 4.1 #######################
### miRNA prediction 				                  ###
#########################################################

MiRNA_prediction=var['MiRNA_prediction']

r=MiRNA_prediction
if (r==True) or (r=='True') or r=='T' or r=='t' :							   ###
	######################################################
	### if True provide following parameters       	   ###
	### mi_genome_file is reference file as used before###
	### for example mi_genome_file=$HOME/reference.fa  ###
	######################################################
	
	mi_genome_file=var['genome_file']
	
	######################################################
	### mi_sequence_file='' if sequences are to be     ###
	### considered directly from the profiles 		   ###
	### else provide properly formatted multi fasta    ###
	### file, read mirDeep-P for format details        ###
	######################################################
	
	mi_sequence_file=var['mi_sequence_file']
	
	######################################################
	### gff_file is an annotation file consisting other###
	### types of folding RNAs than miRNAs. It could be ###
	### left empty if no annotation file is available  ###
	######################################################
	
	gff_file=var['gff_file']

	mirDeep_2=var['mirDeep_2']
	r=mirDeep_2
	if (r==True) or (r=='True') or r=='T' or r=='t' :
		mirDeep_2 = True
		
	miR_known_mature=var['miR_known_mature']
	miR_related_mature=var['miR_related_mature']
	miR_known_precursor=var['miR_known_precursor']
	specie_name=var['specie_name']
	
##########################################################
###################### Module 4.2 ########################
### this module maps the pedicted miRNAs to conserved  ###
### and know miRNAs 								   ###
##########################################################

MiRNA_dataset_mapping=var['MiRNA_dataset_mapping']

r=MiRNA_dataset_mapping
if (r==True) or (r=='True') or r=='T' or r=='t' :
	######################################################
	### if True provide following parameters       	   ###
	### mapping_infile='' if miRNA prediction file is  ###
	### generated from module 4.1 else provide a 	   ###
	### miRDeep-P format outfile 					   ###
	######################################################
	
	mapping_infile=var['mapping_infile']
	
	######################################################
	### miRNA_database is set of conserved miRNAs      ###
	### miRBase release 18 is attached with demo       ###
	######################################################
	
	miRNA_database=var['miRNA_database']

##########################################################
###################### Module 5 ##########################
### this module predicts the ta-siRNAs                 ###
##########################################################

TasiRNA_prediction=var['TasiRNA_prediction']

r=TasiRNA_prediction
if (r==True) or (r=='True') or r=='T' or r=='t' :						   ###
	######################################################
	### if True provide following parameters       	   ###
	### genome_file is reference file as used before   ###
	### for example genome_file=$HOME/reference.fa  ###
	######################################################

	genome_file=var['genome_file']
	
	######################################################
	### tasi_mapping_infile='' if mapped file is       ###
	### generated in module-2                          ###
	######################################################
	
	tasi_mapping_infile=var['igv_infile']
	
	######################################################
	### tasi_cluster_infile='' if the cluster file is  ###
	### same as generated in module -3.1 			   ###
	######################################################
	
	
	tasi_cluster_infile=var['cluster_infile']



##########################################################
####################### Module 6 #########################
### module adds annotations to small RNA profiles      ###
### store in the MySQL table						   ###
##########################################################

Add_annotation_to_mysql_database=var['Add_annotation_to_mysql_database']

r=Add_annotation_to_mysql_database

if (r==True) or (r=='True') or r=='T' or r=='t' :			   ###
	######################################################
	### if True provide following parameters           ###
	### for example server='dna-54-229.mb.au.dk'       ###	
	### for example user= $HOME					       ###
	### for example password=12345				       ###
	### database must be named as 'profiles'	       ###
	######################################################

	server=var['server']
	user=var['user']
	password=var['password']
	database=var['database']
	
	######################################################
	### provide annotation files, separated by commas  ###
	### these file could host genome, homologous genome###
	### repeat dataset, miRNA dataset or any other     ###
	### fatsa files, annotation coloumn is added and   ###
	### header represent name of the file			   ###
	######################################################
	annotation_file=[anno for anno in var['annotation_file'].split(',')]
	
	######################################################
	### this module adds the annotation as exonic,     ###
	### intronic or intergenic depending upon the      ###
	### mapping position 							   ###
	######################################################
	
	add_genomic_region=var['add_genomic_region']
	
	r=add_genomic_region
	if (r==True) or (r=='True') or r=='T' or r=='t' :						   ###
		######################################################
		### if True provide following parameters       	   ###
		### gtf file is gene model and must be given       ###
		######################################################
		add_genomic_region=True
		gtf_file=var['gtf_file']


##############################################################
######################## Module 7 ############################
### generates plots 						               ###
##############################################################


##############################################################
######################### Module 7.1 #########################
### generates library wise comparison plots                ###
### by default all the comparison plots will be made       ###
### by default every second library will consider the first### 
### one as reference, otherwise please change files        ###
### compare_libraries and libraries_reference, after the   ###
### first run			                                   ###
##############################################################
 
 
Pattern_plot=var['Pattern_plot']

r=Pattern_plot
if (r==True) or (r=='True') or r=='T' or r=='t' :
	Pattern_plot=True
	##########################################################
	### if True provide following parameters       	       ###
	### score_regulation_cutoff, to plot differentially    ###
	### regulated sequences score is as defined in         ###
	### module-4, put 0 to consider all the sequences      ###
	##########################################################
		
	score_regulation_cutoff=int(var['score_regulation_cutoff'])
	
	##########################################################
	### pattern_plot_file='' if same as produced in        ###
	### module-3, else provide a profile table             ###
	##########################################################
	
	pattern_plot_file=var['pattern_plot_file']



##############################################################
######################## Module 8  ###########################
### Pattern search 						                   ###
##############################################################

 
##############################################################
######################## Module 8.1  #########################
### this module uses hidden morkov model to calcuate       ###
### frequency of next nucleotides given previous one       ###
### a morkov probability output file is generated          ###
##############################################################

Find_pattern=var['Find_pattern']

r=Find_pattern
if (r==True) or (r=='True') or r=='T' or r=='t' :
	Find_pattern=True
	##########################################################
	### if True provide following parameters       	       ###
	### infile_for_pattern='' if same as produced in       ###
	### module-1.3, else provide a profile table           ###
	##########################################################

	infile_for_pattern=var['infile_for_pattern']
	
	##########################################################
	## ref_file_for_pattern is reference file as used before #
	## for example ref_file_for_pattern=$HOME/reference.fa  ##
	##########################################################
	
	ref_file_for_pattern=var['genome_file']


##############################################################
######################## Module 8.2  #########################
### map all the sequences with a starting pattern          ###
##############################################################

Map_pattern=var['Map_pattern']

r=Map_pattern
if (r==True) or (r=='True') or r=='T' or r=='t' :									   ###
	##########################################################
	### if True provide following parameters       	       ###
	### genome_file is reference file as used before       ###
	### for example genome_file=$HOME/reference.fa         ###
	##########################################################
	
	genome_file=var['genome_file']
	
	##########################################################
	### infile_for_pattern='' if same as produced in       ###
	### module-3, else provide a profile table             ###
	##########################################################
	
	infile_for_pattern_map=var['infile_for_pattern_map']
	
	###########################################################
	### pattern is the set of nucleotide which are present  ###
	### in the start of sequence                            ###
	###########################################################
	
	pattern=var['pattern']


##############################################################
######################## Module 9  ###########################
### making mySQL queries 						           ###
##############################################################


MySQL_query=var['MySQL_query']

r=MySQL_query
if (r==True) or (r=='True') or r=='T' or r=='t' :									   ###
	##########################################################
	### if True provide following parameters       	       ###
	### make a file with name 'query_infile' and write     ###
	### MySQL syntax query, it will be processed as output ###
	### will be send to text file                          ###
	##########################################################
	
	query_infile=var['query_infile']
	
	######################################################
	### if True provide following parameters           ###
	### for example server='dna-54-229.mb.au.dk'       ###	
	### for example user= $HOME					       ###
	### for example password=12345				       ###
	### database must be named as 'profiles'	       ###
	######################################################

	server=var['server']
	user=var['user']
	password=var['password']
	database=var['database']


