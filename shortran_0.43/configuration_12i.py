###################################################
##### shortran interactive configuration file #####
###################################################

###################################################
### Date: 2012-06-11                            ###
### Auther: Vikas Gupta, Stig Anderson          ###
###################################################

###################################################
### importing modules							###
import sys,time,os,datetime						###
work_dir=os.getcwd()							###
miRNA_dir=work_dir+"/miRDP1.3"					###
###################################################


##########################################################
### Following is a list of Modules 					   ###
### to use (add "True" if using else "False")		   ###
### please remember that "True" and "False" are case   ###
### sensitive										   ###
##########################################################

### check if number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def to_run(module):
	Raw_data_filtering=False
	ErrorCorrection=False
	Genome_mapping=False
	Make_profile_files=False
	Add_output_to_mysql_database=False
	Cluster_prediction=False
	Add_cluster_to_mySQL=False
	Cluster_genomic_region_prediction=False
	MiRNA_prediction=False
	MiRNA_dataset_mapping=False
	TasiRNA_prediction=False
	Add_annotation_to_mysql_database=False
	Pattern_plot=False
	Find_pattern=False
	Map_pattern=False
	MySQL_query=False

	token=module.split(',')
	for i in range(len(token)):
		if(len(module)>0):
			try:
				m=float(token[i])
				if(m==1.1):
					Raw_data_filtering=True
				if(m==1.2):
					Make_profile_files=True
				if(m==1.3):
					Add_output_to_mysql_database=True
				if(m==2):
					Genome_mapping=True
				if(m==3.1):
					Cluster_prediction=True
				if(m==3.2):
					Cluster_genomic_region_prediction=True
				if(m==4.1):
					MiRNA_prediction=True
				if(m==4.2):
					MiRNA_dataset_mapping=True
				if(m==5):
					TasiRNA_prediction=True
				if(m==6):
					Add_annotation_to_mysql_database=True
				if(m==7):
					Pattern_plot=True
				if(m==8.1):
					Find_pattern=True
				if(m==8.2):
					Map_pattern=True
				if(m==9):
					MySQL_query=True
			except:
				print "Please enter module number for help or press E/e for exit"
				print
				take=raw_input()
				if is_number(take):
					help(float(take))
	return (Raw_data_filtering,ErrorCorrection,Genome_mapping,Make_profile_files,Add_output_to_mysql_database,Cluster_prediction,Add_cluster_to_mySQL,Cluster_genomic_region_prediction,MiRNA_prediction,MiRNA_dataset_mapping,TasiRNA_prediction,Add_annotation_to_mysql_database,Pattern_plot,Find_pattern,Map_pattern,MySQL_query)


### define help function
def module_run():
	### ask for modules needed to be run
	print "Please choose modules to run"
	print "if multiple modules need to be selected then please separate these by commas(,)"
	print
	module=raw_input("==>")
	if (module=='e' or module=='E'):
		sys.exit("Program has been terminated.")
	(Raw_data_filtering,ErrorCorrection,Genome_mapping,Make_profile_files,Add_output_to_mysql_database,Cluster_prediction,Add_cluster_to_mySQL,Cluster_genomic_region_prediction,MiRNA_prediction,MiRNA_dataset_mapping,TasiRNA_prediction,Add_annotation_to_mysql_database,Pattern_plot,Find_pattern,Map_pattern,MySQL_query)=to_run(module)
	return (Raw_data_filtering,ErrorCorrection,Genome_mapping,Make_profile_files,Add_output_to_mysql_database,Cluster_prediction,Add_cluster_to_mySQL,Cluster_genomic_region_prediction,MiRNA_prediction,MiRNA_dataset_mapping,TasiRNA_prediction,Add_annotation_to_mysql_database,Pattern_plot,Find_pattern,Map_pattern,MySQL_query)

def help(module):
	if(module==1.1):
		print '''
		##########################################################
		################### Module 1 #############################
		## Adapter filtering, Read trimming, Size fractionation ##
		### Sequencing error correction, homology filtering    ###
		##########################################################
		### Here the input is raw fastq files                  ###
		## You should provide each library in a separate folder ##
		### for example: $HOME/libraries                       ###
		### For more details see the demo data                 ###
		##########################################################
		
		##########################################################
		### Please write folder names separated by commas      ###
		### for example: "lib-1","lib-2"                       ###
		##########################################################
		
		##########################################################
		### Define size range for small RNAs to considered for ###
		### creating profile tables                            ###
		##########################################################
		
		##########################################################
		### Define read abundance cut-off for profile tables   ###
		### Here read abundance cut-off refers to the minimum  ###
		### sum of read abundances across all libraries        ###
		##########################################################
		
		##########################################################
		### Define header for each library.                    ###		
		### It will be used in plots and MySQL table legends.  ###	
		##########################################################
		
		##########################################################
		### Provide file for homology filtering                ###
		### for example filtering_file=$HOME/repeat.fa         ###
		##########################################################

		##########################################################
		### If reads contain many unidentified bases 'N'       ### 
		### an error correction is recommended.                ###
		### Please remember that it takes a significant amount ###
		### of time (see manual for more details)              ###
		##########################################################

		'''
		
	
	if (module==1.2):
		print '''
		##########################################################
		###################### Module 1.2 ########################
		### Normalization and profiling                        ###
		### This sub-module generates a normalized profile     ###
		### table from the fastq files                         ###
		##########################################################
		'''
		
	if (module==1.3):
		print '''
		##########################################################
		###################### Module 1.3 ########################
		### Normalization and profiling                        ###
		### This sub-module connects to the MySQL server       ###
		### and stores the data in tables with the same        ###
		### name as the profile table                          ###
		##########################################################
		
		'''
		
	if (module==2):
		
		print '''
		##########################################################
		###################### Module 2 ##########################
		###                Genomic Mapping                     ###
		##########################################################
		
		##########################################################
		### Provide a genome file                              ###
		### For example: $HOME/reference.fa                    ###
		##########################################################
		
		##########################################################
		### If the genomic mapping module needs to be used     ###
		### independent of module-1 output, please provide     ###
		### a fasta file to be mapped                          ###
		### For example: $HOME/reads.fa                        ###
		##########################################################
		'''
		
	if (module==3.1):
		print '''
		##########################################################
		###################### Module 3 ##########################
		###                   Clustering                       ###
		##########################################################
		
		###################### Module 3.1 ########################
		### Outputs small RNA clusters                         ###
		##########################################################
		
		##########################################################
		### If divide_by_size is True, it invokes a            ###
		### size fractionated clustering algorithm. Otherwise  ###
		### it considers all the sequences in the mapped .igv  ###
		### file. Additional variables can be provided         ###
		### abundance_cutoff: minimal sum of read abundances   ###
		### max_gap_size: Maximum length of a gap in cluster   ###
		### min_unique_reads: Minimum number of unique reads   ###
		##########################################################
		
		##########################################################
		### igv_infile is the file produced by Module 2        ###
		##########################################################
		
		##########################################################
		### Add_cluster_to_mySQL                               ###
		### Adds clusters information to the MySQL table       ###
		##########################################################
		'''		
		
	if (module==3.2):
		print '''
		
		##########################################################
		####################### Module 3.2 #######################
		### This modules annotates clusters with genomic       ###
		### region information.                                ###
		### It also calculates the number of exonic, intronic  ###
		### and intergenic clusters				        	   ###
		##########################################################

		'''
		
	if (module==4.1):
		print '''
		
		##########################################################
		###################### Module 4.1 ########################
		###                 miRNA prediction 		           ###
		##########################################################
		
		##########################################################
		### mi_genome_file is reference genome file            ###
		### For example: $HOME/reference.fa                    ###
		##########################################################
		
		##########################################################
		### mi_sequence_file='' if sequences are to be         ###
		### analyzed directly from the profiles.               ###
		### Else provide properly formatted multi-fasta        ###
		### file, read mirDeep-P/mirDeep2 documentation        ###
		### for format details                                 ###
		##########################################################
		
		##########################################################
		### gff_file is an annotation file consisting other    ###
		### types of folding RNAs than miRNAs. It could be     ###
		### left empty if no annotation file is available      ###
		##########################################################
		'''
		
	if (module==4.2):
		print '''
		##########################################################
		###################### Module 4.2 ########################
		### This module maps the pedicted miRNAs to conserved  ###
		### and known miRNAs                                   ###
		##########################################################
		
		##########################################################
		### if True provide the following parameters           ###
		### mapping_infile='' if miRNA prediction file is      ###
		### generated from module 4.1 else provide a           ###
		### miRDeep-P format outfile                           ###
		##########################################################
		
		##########################################################
		### miRNA_database is set of conserved miRNAs          ###
		### miRBase release 18 is included in the demo data    ###
		##########################################################
		'''
	if (module==5):
		print '''
		##########################################################
		###################### Module 5 ##########################
		### This module predicts plant Ta-siRNAs               ###
		##########################################################
		
		######################################################
		### if True provide the following parameters       ###
		### genome_file is the reference genome file       ###
		### for example genome_file=$HOME/reference.fa     ###
		######################################################
		
		######################################################
		### tasi_mapping_infile='' to use mapping file     ###
		### generated by Module 2                          ###
		######################################################
		
		######################################################
		### tasi_cluster_infile='' to use the cluster file ###
		### generated by Module 3                          ###
		######################################################
        '''
			
	if (module==6):
		print '''
		##########################################################
		####################### Module 6 #########################
		### Module adds annotations to small RNA profiles  and ###
		### stores them in the MySQL table                     ###
		##########################################################
		
		##########################################################
		### Provide fasta files for annotation,                ###
		### separated by commas.                               ###
		### These files could be gene sequences,               ###
		### repeats, miRNAs or any other fasta files.          ###
		### For each file, an annotation column is added with  ###
		### file names as column headers                       ###
		##########################################################
		
		##########################################################
		### This module also adds genomic region annotation -  ###
		### exonic, intronic or intergenic - depending         ###
		### upon the mapping positions                         ###
		##########################################################
		
		######################################################
		### if True provide following parameters:          ###
		### gtf_file with gene model information           ###
		######################################################
		'''
	
	if (module==7):
		print '''
		##############################################################
		######################## Module 7 ############################
		### Generates library comparison plots                     ###
		##############################################################
		
		##############################################################
		######################### Module 7.1 #########################
		### generates library wise comparison plots                ###
		### by default all the comparison plots will be made       ###
		### by default every second library will consider the prior### 
		### one as reference, otherwise please change files        ###
		### compare_libraries and libraries_reference, after       ###
		### running module 1                                       ###
		##############################################################
		
		##########################################################
		### if True provide following parameters:      	       ###
		### score_regulation_cutoff, to plot differentially    ###
		### regulated sequences. Score is as defined in        ###
		### Module 3, put 0 to consider all the sequences      ###
		##########################################################
		
		##########################################################
		### pattern_plot_file='' if same as produced in        ###
		### module-1, else provide a profile table             ###
		##########################################################
		'''
	if (module==7.2):
		print '''
		##############################################################
		######################## Module 7 ############################
		### generates plots                                        ###
		##############################################################
		
		##############################################################
		######################### Module 7.2 #########################
		### generates plot of different size fraction of reads     ###
		### present in the dataset                                 ###
		### you can provide a profile table, if you do not wan to  ###
		### use the default one generated by Module 1              ### 
		##############################################################
		'''
	
	if (module==8.1):
		print '''
		##############################################################
		######################## Module 8  ###########################
		### Pattern search                                         ###
		##############################################################
		
		 
		###############################################################
		######################## Module 8.1  ##########################
		### This module uses a hidden Markov model to calculate     ###
		### the frequency of the next nucleotide given the previous ###
		### A markov probability output file is generated           ###
		###############################################################
		
		##########################################################
		### if True provide the following parameters           ###
		### infile_for_pattern='' to use the one generated by  ###
		### Module 3, else provide a profile table             ###
		##########################################################
		
		##########################################################
		## ref_file_for_pattern is the genome reference file     #
		## for example ref_file_for_pattern=$HOME/reference.fa  ##
		##########################################################
		'''
	if (module==8.2):
		print '''
		##############################################################
		######################## Module 8  ###########################
		### Pattern search                                         ###
		##############################################################
		
		
		##############################################################
		######################## Module 8.2  #########################
		### map all the sequences with a starting pattern          ###
		##############################################################
		
		##########################################################
		### if True provide following parameters       	       ###
		### genome_file is the genome reference file           ###
		### for example genome_file=$HOME/reference.fa         ###
		##########################################################
		
		##########################################################
		### infile_for_pattern='' to use the one generated by  ###
		### Module 3, else provide a profile table             ###
		##########################################################
		
		###########################################################
		### pattern is a set of nucleotides which are present   ###
		### at the 5' end of the sequence                       ###
		###########################################################
		'''
		
	if (module==9):
		print '''
		##############################################################
		######################## Module 9  ###########################
		### making mySQL queries                                   ###
		##############################################################
		
		##########################################################
		### if True provide following parameters               ###
		### make a file with name 'query_infile' and write     ###
		### MySQL syntax query. It will be processed and       ###
		### the output will be sent to a text file             ###
		##########################################################
		'''
		
	(Raw_data_filtering,ErrorCorrection,Genome_mapping,Make_profile_files,Add_output_to_mysql_database,Cluster_prediction,Add_cluster_to_mySQL,Cluster_genomic_region_prediction,MiRNA_prediction,MiRNA_dataset_mapping,TasiRNA_prediction,Add_annotation_to_mysql_database,Pattern_plot,Find_pattern,Map_pattern,MySQL_query)=module_run()
	
	return (Raw_data_filtering,ErrorCorrection,Genome_mapping,Make_profile_files,Add_output_to_mysql_database,Cluster_prediction,Add_cluster_to_mySQL,Cluster_genomic_region_prediction,MiRNA_prediction,MiRNA_dataset_mapping,TasiRNA_prediction,Add_annotation_to_mysql_database,Pattern_plot,Find_pattern,Map_pattern,MySQL_query)
	
def to_run(module):
	Raw_data_filtering=False
	ErrorCorrection=False
	Genome_mapping=False
	Make_profile_files=False
	Add_output_to_mysql_database=False
	Cluster_prediction=False
	Add_cluster_to_mySQL=False
	Cluster_genomic_region_prediction=False
	MiRNA_prediction=False
	MiRNA_dataset_mapping=False
	TasiRNA_prediction=False
	Add_annotation_to_mysql_database=False
	Pattern_plot=False
	Find_pattern=False
	Map_pattern=False
	MySQL_query=False

	token=module.split(',')
	for i in range(len(token)):
		if(len(module)>0):
			try:
				m=float(token[i])
				if(m==1.1):
					Raw_data_filtering=True
				if(m==1.2):
					Make_profile_files=True
				if(m==1.3):
					Add_output_to_mysql_database=True
				if(m==2):
					Genome_mapping=True
				if(m==3.1):
					Cluster_prediction=True
				if(m==3.2):
					Cluster_genomic_region_prediction=True
				if(m==4.1):
					MiRNA_prediction=True
				if(m==4.2):
					MiRNA_dataset_mapping=True
				if(m==5):
					TasiRNA_prediction=True
				if(m==6):
					Add_annotation_to_mysql_database=True
				if(m==7):
					Pattern_plot=True
				if(m==8.1):
					Find_pattern=True
				if(m==8.2):
					Map_pattern=True
				if(m==9):
					MySQL_query=True
			except:
				print "Please enter module number for help or E/e for exit"
				print
				take=raw_input()
				if is_number(take):
					(Raw_data_filtering,ErrorCorrection,Genome_mapping,Make_profile_files,Add_output_to_mysql_database,Cluster_prediction,Add_cluster_to_mySQL,Cluster_genomic_region_prediction,MiRNA_prediction,MiRNA_dataset_mapping,TasiRNA_prediction,Add_annotation_to_mysql_database,Pattern_plot,Find_pattern,Map_pattern,MySQL_query)=help(float(take))
	return (Raw_data_filtering,ErrorCorrection,Genome_mapping,Make_profile_files,Add_output_to_mysql_database,Cluster_prediction,Add_cluster_to_mySQL,Cluster_genomic_region_prediction,MiRNA_prediction,MiRNA_dataset_mapping,TasiRNA_prediction,Add_annotation_to_mysql_database,Pattern_plot,Find_pattern,Map_pattern,MySQL_query)
### ask for modules needed to be run
print "Please choose the modules you would like to run"
print "if multiple modules need to be selected then please separate these by commas(,)"
print
module=raw_input("==>")
if (module=='e' or module=='E'):
	sys.exit("Program has been terminated.")
Raw_data_filtering,ErrorCorrection,Genome_mapping,Make_profile_files,Add_output_to_mysql_database,Cluster_prediction,Add_cluster_to_mySQL,Cluster_genomic_region_prediction,MiRNA_prediction,MiRNA_dataset_mapping,TasiRNA_prediction,Add_annotation_to_mysql_database,Pattern_plot,Find_pattern,Map_pattern,MySQL_query=to_run(module)




		
### ask for mySQL server
def mysql():
	######################################################
	### if True provide following parameters           ###
	### for example server='localhost'                 ###	
	### for example user= user					       ###
	### for example password=   				       ###
	### database must be named as 'profiles'	       ###
	######################################################

	print "Provide MySQL server address:"
	print "even if installed on the same computer, it still has a server address"
	print "press (Enter/Return) for Default"
	print "else provide server address:"
	print "if defualt then server='localhost'"
	print
	take=raw_input()
	if(take==''):
		server='localhost'
	else:
		server=take
	
	print "Provide user name:"
	print "if defualt user='user'"
	print
	take=raw_input()
	if (take==''):
		user='user'
	else:
		user=take
	
	print "Provide password"
	print "if defualt password=''"
	print
	take=raw_input()
	if(take==''):
		password=''
	else:
		password=take
	
	print "Provide database"
	print "if defualt database=profiles"
	print
	take=raw_input()
	if(take==''):
		database='profiles'
	else:
		database=take

	return (server,user,password,database)

server=''
print "Would you like to configure theMySQL server ?"
print "y for Yes or n for No"
print
take=raw_input()
if(take=='y' or take=='Y'):
	server,user,password,database=mysql()


print "Provide minimum read length to be considered for profiling"
print "Press (Enter/Return) for Default"
print "Default is 19"
print
take=raw_input()
if(take==''):
	min_seq_size=19
else:
	min_seq_size=int(take)
print "Provide maximum read length to be considered for profiling"
print "press (Enter/Return) for Default"
print "Default is 24"
print
take=raw_input()
if(take==''):
	max_seq_size=24
else:
	max_seq_size=int(take)

##########################################################
### read abundace cut-off for profiles, cut-off refers ### 
### to sum of read abundance is all libraries		   ###
##########################################################
print "Provide read abundace cut-off for profile tables"
print "Press (Enter/Return) for Default"
print "Default is 1"
print
take=raw_input()
if(take==''):
	cutoff=1
else:
	cutoff=int(take)

##########################################################
### read score    cut-off for profiles, cut-off refers ### 
### to sum of read abundance is all libraries		   ###
##########################################################
print "Provide score regulation cut-off for profile tables"
print "Press (Enter/Return) for Default"
print "Default is 1"
print
take=raw_input()
if(take==''):
	score_regulation_cutoff=1
else:
	score_regulation_cutoff=int(take)

##########################################################
##################  reference genome   ###################
##########################################################

print "Provide a Reference Genome file"
print "Press (Enter/Return) for Default"
print "Default is"+ work_dir+"/demo/annotation_files/lj_r25_demo.fa"
print
take=raw_input()
if(take==''):
	genome_file=work_dir+'/demo/annotation_files/'+"lj_r25_demo.fa"
else:
	genome_file=take

##########################################################
### inputs are raw fastq files 						   ###
### here you should provided each library with all     ###
### sizes in seperate folder check example dataset for ###
### more details for demo							   ###
### for example: fq_file_dir=$HOME/libraries           ###
##########################################################

print "Provide a path to the folder, which contains all libraries"
print "Press (Enter/Return) for Default"
print "Default is ("+work_dir+"/demo)"
print
take=raw_input()
if(take==''):
	fq_file_dir=work_dir+'/demo'
else:
	fq_file_dir=take


##########################################################
### please put folder name seperated by commas         ###
### for example: lib= "lib-1","lib-2"                  ###
##########################################################

print "Provide library folder names separated by commas"
print "Press (Enter/Return) for Default"
print "Default is 'GHD-5_sample','GHD-6_sample','GHD-9_sample','GHD-10_sample'"
print
take=raw_input()
if(take==''):
	lib="GHD-5_sample","GHD-6_sample","GHD-9_sample","GHD-10_sample"
else:	
	lib=take


##########################################################
### genotypic headers	for each library	           ###		
### which can used in plots legends					   ###	
##########################################################
print "Provide headers for each library seperated by commas"
print "press (Enter/Return) for Default"
print "Default is 'Mock','wild type','Mock','NFR5'"
print
take=raw_input()
if(take==''):
	genotype='Mock','wild type','Mock','NFR5' 
else:
	genotype=take

##########################################################
################### Module 1 #############################
### Filtering and Sequencing Error correction          ###
##########################################################
### define size range of small RNAs to considered for  ###
### creating profiles                                  ###
##########################################################
def mapping_parameter():
	print "Provide mapping parameters: "
	print "Number of mismatches allowed: "
	print "Press (Enter/Return) for Default 2"
	print 
	taken=raw_input()
	if(taken==''):
		mismatches=2
	else:
		mismatches=taken
	print "Maximum Maq error allowed: "
	print "The maximum allowed sum of the quality scores of mismatched bases"
	print "A high quality base has quality 40"
	print "Press (Enter/Return) for Default 70"
	print 
	taken=raw_input()
	if taken=='':
		max_maq_err=70
	else:
		max_maq_err=int(taken)

	print "Seed length : "
	print "Press (Enter/Return) for Default 25"
	print
	taken=raw_input()
	if taken=='':
		seedlen=25
	else:
		seedlen=taken
	
	print "Number of cores to be used: "
	print "Press (Enter/Return) for Default 1"
	print
	taken=raw_input()
	if taken=='':
		no_of_cores=1
	else:
		no_of_cores=taken
		
		
	return mismatches,max_maq_err,seedlen,no_of_cores
	
def default_mapping_parameter():
	mismatches=2
	max_maq_err=70
	seedlen=25
	no_of_cores=1
	
	return mismatches,max_maq_err,seedlen,no_of_cores

use_only_mapped_reads=''
def use_only_mapped():
	print "Would you like to use only reads mapping to the filtering files for further analysis?"
	print "Press (Enter/Return) for Default [No], else press any other key"
	print 
	take=raw_input()
	if(take==''):
		use_only_mapped_reads=False
	else:
		use_only_mapped_reads=True
	return use_only_mapped_reads
	
adapter_filtering=''
size_fractionation=''
##########################################################
##########################################################
### True if module needed to be used                   ###
##########################################################

if (Raw_data_filtering==True) :                   ###
 	print "## Configuring Module 1.1 "
 	print "Would you like to keep default parameters for mapping ?"
	print "Press (enter/return) for default [yes]"
	print "else press any other key for [No]"
	print 
	
	take=raw_input()
	if(take==''): 
		mismatches,max_maq_err,seedlen,no_of_cores = default_mapping_parameter()
	else:
		mismatches,max_maq_err,seedlen,no_of_cores = mapping_parameter()
	
	print "Do Fastq files require adapter filtering?  "
	print "Press (Enter/Return) for defualt[Yes]"
	print "else press any other key for [No]"
	print 
	take=raw_input()
	if(take!=''):
		adapter_filtering=False
	else:
		adapter_filtering=True
		print "Please provide adapter sequence"
		print "if multiple then separate these by commas in the same order as the libraries"
		print "Press (Enter/Return) for default adapter [AATTGGCC] for demo data"
		
		take=raw_input()
		if(take==''):
			adapter='AATTGGCC'
		else:
			adapter=take
			
		
		if adapter != '':
			print "Would you like to keep only sequences which contained the adapter?"
			print "Press (Enter/Return) for defualt[Yes] else anyother key"
			print
			
			keep_with_adapter_only = raw_input()
			
			if(keep_with_adapter_only==''):
				keep_with_adapter_only=True
			else:
				keep_with_adapter_only = False

		
		
		print "Would you like to trim the reads from the 3' end "
		print "Press (Enter/Return) for defualt['0'] else an integer"
		print 
		
		trim = raw_input()
		if(trim==''):
			trim=0
		else:
			trim = trim
	
	print "Are Fastq files are size fractionated ?"
	print "Press (Enter/Return) for defualt[No]"
	print "else press any other key for [Yes]"
	print
	take=raw_input()
	if(take==''):
		size_fractionation=True
	else:
		size_fractionation=False
	
	######################################################
	### if True provide the following parameters       ###
	### provide file for filtering                     ###
	### for example filtering_file=$HOME/repeat.fa     ###
	######################################################
	print "Provide file for filtering"
	print "Press (Enter/Return) for Default"
	print "Press 'n' if homology filtering is not required"
	print "Default is "+work_dir+"/demo/annotation_files/lotus_rep.fa"
	print
	take=raw_input()
	if(take==''):
		repeat_database=work_dir+'/demo/annotation_files/'+'lotus_rep.fa'
	if(take=='n'):
		repeat_database=''
		use_only_mapped_reads=False
	else:
		repeat_database=take
		if use_only_mapped_reads=='':
			use_only_mapped_reads=use_only_mapped()
	
	
	######################################################
	### if sequencing reads contains many unidentified ###
	### based 'N', an Error correction is recommonded  ###
	### please remember it is a fairly time consuming  ###
	### step										   ###
	### True if sub-module needed to be used           ###
	######################################################
	print "Would you like to run the  Short Read Error Correction Algorithm?"
	print "It is very time-consuming for for large datasets"
	print "Press (Enter/Return) for Default else type any other key"
	print "Default is False"
	print
	take=raw_input()
	if(take==''):
		ErrorCorrection=False
	else:
		ErrorCorrection=True
		print 
		print "Please provide path to the executable ECHO scripts"	
		print
		take=raw_input()
		if(take==''):
			print 'Path was not configured... skipping Error correction'
			ErrorCorrection=False
		else:
			Path_2_ECHO=take
	
##########################################################
###################### Module 1.2 ########################
### Normalization and profiling                        ###
### True if module needed to be used                   ###
### submodule generates a normalized profile    	   ###
### table from the fastq files                         ###
### True if module needed to be used			       ###
##########################################################
if(Make_profile_files==True):	
	print "## Configuring Module 1.2 "
	if use_only_mapped_reads=='':
		use_only_mapped_reads=use_only_mapped()
	if adapter_filtering == '':
		adapter_filtering=False
	if size_fractionation== '':
		size_fractionation=False
	print 
##########################################################
################### Module 1.3 ###########################
### submodule connects with the MySQL server    	   ###
### and stores data in tables with name same as 	   ###
### the profile output file                     	   ###
##########################################################


if(Add_output_to_mysql_database==True):	
	######################################################
	### make sure that normalized expression values    ###
	### file is in output foler						   ###	
	######################################################
	print "## Configuring Module 1.3 ## "
	if (server==''):
		server,user,password,database=mysql()

##########################################################
###################### Module 2 ##########################
### Genomic Mapping                                    ###
### True if module needed to be used                   ###
##########################################################


if(Genome_mapping==True):							   ###
	######################################################
	### if True provide following parameters           ###
	### for example genome_file=$HOME/reference.fa     ###
	######################################################
			
	print "## Configuring Module 2 ## "	
	print "Would like to keep default parameters mapping ?"
	print "Press (enter/return) for default [yes]"
	print "else press any other key for [No]"
	print 
	
	take=raw_input()
	if(take==''): 
		mismatches,max_maq_err,seedlen,no_of_cores = default_mapping_parameter()
	else:
		mismatches,max_maq_err,seedlen,no_of_cores = mapping_parameter()


	
	######################################################
	### provide a input fasta file if it defers from   ###
	### the previous modules else leave it empty       ### 
	### for example infile=$HOME/reads.fa              ###
	######################################################
	
	print "Provide read file to mapped"
	print "if defualt then make sure that fasta file from module-2 is available"
	print "Press (Enter/Return) for Default"
	print "Default is the one generated by Module 2"
	print
	take=raw_input()
	if(take==''):
		infile=''
	else:
		infile=take
	

##########################################################
###################### Module 3 ##########################
### Clustering 				                           ###
### True if module needed to be used                   ###
##########################################################


##########################################################
###################### Module 3.1 ########################
### outputs small RNAs clusters 					   ###  
##########################################################


if(Cluster_prediction==True):
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
	print "## Configuring Module 3.1 ## "	
	print "Would you like to run NiBLS for cluster prediction"
	print "or use the default algorithm with the specified parameters"
	print "For NiBLS, press any key except (Enter/Return)"
	print
	
	take=raw_input()
	
	if(take==''):
		NiBLS=False
		
		print "Should cluster prediction need to be performed on each size separately?"
		print "if True then press any key except (Enter/Return)"
		print "Press (Enter/Return) for Default"
		print "Default is generated from module -3"
		print
		take=raw_input()
		
		if(take==''):
			divide_by_size=False
		else:
			divide_by_size=True
		######################################################
		### abundance cut-off is the minimal read abundance###
		### in a gemonic region which doesn't contain a gap###
		### larger than 50 base-pair and has atleast 5     ###
		### unique sequences							   ###
		######################################################
		
		print "Abundance cutoff for clustering algorithm"
		print "if not default provide a positive integer"
		print "Press (Enter/Return) for Default"
		print "Default is 75"
		print
		take=raw_input()
		if(take==''):
			abundance_cutoff=75 
		else:
			abundance_cutoff=int(take)
			
		print "Maximum gap size for clustering algorithm"
		print "if not default provide a positive integer"
		print "Press (Enter/Return) for Default"
		print "Default is 50"
		print
		take=raw_input()
		if(take==''):
			max_gap_size = 50
		else:
			max_gap_size = int(take)
		
		print "Minimum number of unique sRNAs for clustering algorithm"
		print "if not default provide a positive integer"
		print "Press (Enter/Return) for Default"
		print "Default is 5"
		print
		take=raw_input()
		if(take==''):
			min_unique_reads = 5
		else:
			min_unique_reads = int(take)

		
		
		
	else:
		NiBLS=True
	
	
	######################################################
	### igv_infile is mapped file produced in module-2 ###
	### leave igv_infile='' 						   ###
	### if file used is same from module-2 			   ###
	### else provide a igv format mapped file		   ###
	######################################################
	
	print "Provide a standard file with mapping co-ordinates,i.e., .bed or .igv format"
	print "press (Enter/Return) for Default"
	print "Default is file generated from module-2"
	print "For Demo, please provide file :"+work_dir+"/demo/annotation_files/sample.igv"
	print
	take=raw_input()
	if(take==''):
		igv_infile='' 
	else:
		igv_infile=take
	
	######################################################
	### True if module needed to be used               ###
	### it will setup table for clusters in the MySQL  ###
	######################################################
	
	print "Should clusters be added to MySQL Server?"
	print "if yes then press any key except (Enter/Return) "
	print "Press (Enter/Return) for no"
	print
	take=raw_input()
	if(take==''):
		Add_cluster_to_mySQL=False
	else:
		Add_cluster_to_mySQL=True
		if (server==''):
			server,user,password,database=mysql()		



		
		
##########################################################
####################### Module 3.2 #######################
### this modules should ebe invoked if genomic region  ###
### for cluster is need to be annotated, it also       ###
### it also calculates number of exonic, intronic and  ###
### intergenic clusters								   ###
### True if module needed to be used                   ###
##########################################################


if(Cluster_genomic_region_prediction==True):		   ###
	######################################################
	### if True provide following parameters       	   ###
	### gtf_file is a genemodel where exons introns are###
	### given with genomic co-ordinate check demo file ###
	### for more details                               ###
	######################################################
	print "## Configuring Module 3.2 ## "	
	print "Provide an annotation (i.e.,.gtf,.gff) file"
	print "Press (Enter/Return) for Demo"
	print "Default is "+work_dir+'/demo/annotation_files/'+"sample.gtf"
	print
	take=raw_input()
	if(take==''):
		gtf_file=work_dir+'/demo/annotation_files/'+"sample.gtf" 
	else:
		gtf_file=take
	######################################################
	### use cluster_infile='' if cluster were produced ###
	### from module 3.1 else use the full path of 	   ###
	### cluster file								   ###
	######################################################

	print "Provide a cluster file"
	print "Press (Enter/Return) for Default"
	print "Default is generated from module -3.1"
	print "For Demo, please provide file :"+work_dir+"/demo/annotation_files/sample.igv_cluster.bed"
	print
	take=raw_input()
	if(take==''):
		cluster_infile=''
	else:
		cluster_infile=take
	
##########################################################
###################### Module 4 ##########################
### miRNA prediction 				                   ###
### True if module needed to be used                   ###
##########################################################


#########################################################
###################### Module 4.1 #######################
#########################################################



if(MiRNA_prediction==True):							   ###
	######################################################
	### if True provide following parameters       	   ###
	### mi_genome_file is reference file as used before###
	### for example mi_genome_file=$HOME/reference.fa  ###
	######################################################
	print "## Configuring Module 4.1 ## "	
	print "Provide a Genome Reference file"
	print "Press (Enter/Return) for Default"
	print "Default for demo is"+ work_dir+"/demo/annotation_files/lj_r25_demo.fa"
	print
	take=raw_input()
	if(take==''):
		mi_genome_file=work_dir+'/demo/annotation_files/'+"lj_r25_demo.fa"
	else:
		mi_genome_file=take
	
	######################################################
	### mi_sequence_file='' if sequences are to be     ###
	### considered directly from the profiles 		   ###
	### else provide properly formatted multi fasta    ###
	### file, read mirDeep-P for format details        ###
	######################################################
	
	print "Provide appropriately formatted fasta file"
	print "Press (Enter/Return) for Default"
	print "Default is generated the one generated by Module 2"
	print
	take=raw_input()
	if(take==''):
		mi_sequence_file=''
	else:
		mi_sequence_file=take
	
	######################################################
	### gff_file is an annotation file consisting other###
	### types of folding RNAs than miRNAs. It could be ###
	### left empty if no annotation file is available  ###
	######################################################
	
	print "Provide an annotation (i.e.,.gtf,.gff) file"
	print "Press (Enter/Return) for Demo"
	print "Default is acceptable empty file"
	print
	take=raw_input()
	if(take==''):
		gff_file='' 
	else:
		gff_file=take

	print "Would you like to run mirDeep 2.0 for miRNA predictions"
	print "or use default mirDeep-P [specific for plant miRNAs]"
	print "For mirDeep 2.0, press any key except (Enter/Return)"
	print
	
	take=raw_input()
	
	if(take==''):
		mirDeep_2=False
	else:
		mirDeep_2=True
		
		print "Please provide the species name [i.e. C.elegans]"
		print "else press (Enter/Return) for demo"
		print
		
		take=raw_input()
		
		if (take == 'd'):
			specie_name='L. japonicus'
		else:
			specie_name=take
			
		print "Please provide a fasta file with known mature miRNAs from "+ specie_name
		print "if demo press (Enter/Return) "
		print "if file not available press 'n' "
		print
		
		take=raw_input()
		
		if take == '':
			miR_known_mature=work_dir+'/demo/annotation_files/'+"lja_mature.fa"
		elif take == 'n':
			miR_known_mature='none'
		else:
			miR_known_mature=take
	
		print "Please provide a fasta file with closely related species mature miRNAs for"+ specie_name
		print "if demo press (Enter/Return) "
		print "if file not available press 'n' "
		print
		
		take=raw_input()
		
		if take == '':
			miR_related_mature=work_dir+'/demo/annotation_files/'+"ath_mature.fa"
		elif take == 'n':
			miR_related_mature='none'
		else:
			miR_related_mature=take
		
		print "Please provide a fasta file with known precursor miRNAs from "+ specie_name
		print "if demo press (Enter/Return) "
		print "if file not available press 'n' "
		print
		
		take=raw_input()
		
		if take == '':
			miR_known_precursor=work_dir+'/demo/annotation_files/'+"lja_stem_loop.fa"
		elif take == 'n':
			miR_known_precursor='none'
		else:
			miR_known_precursor=take
		
		
		
		
##########################################################
###################### Module 4.2 ########################
### this module maps the pedicted miRNAs to conserved  ###
### and know miRNAs 								   ###
##########################################################

if(MiRNA_dataset_mapping==True):
	######################################################
	### if True provide following parameters       	   ###
	### mapping_infile='' if miRNA prediction file is  ###
	### generated from module 4.1 else provide a 	   ###
	### miRDeep-P format outfile 					   ###
	######################################################
	print "## Configuring Module 4.2 ## "	
	print "Provide a fasta file to be mapped on miRNA database"
	print "Press (Enter/Return) for Default"
	print "Default is file generated from module-4.1"
	print
	take=raw_input()
	if(take==''):
		mapping_infile=''
	else:
		mapping_infile=raw_input()
	######################################################
	### miRNA_database is set of conserved miRNAs      ###
	### miRBase release 18 is attached with demo       ###
	######################################################
	
	print "Provide a fasta file of miRNA database"
	print "Press (Enter/Return) for Default"
	print "Default is "+work_dir+'/demo/annotation_files/'+'miRBase.fa'
	print
	take=raw_input()
	if(take==''):
		miRNA_database=work_dir+'/demo/annotation_files/'+'miRBase.fa'
	else:
		miRNA_database=take

##########################################################
###################### Module 5 ##########################
### this module predicts the ta-siRNAs                 ###
##########################################################


if(TasiRNA_prediction==True):						   ###
	######################################################
	### if True provide following parameters       	   ###
	### genome_file is reference file as used before   ###
	### for example genome_file=$HOME/reference.fa  ###
	######################################################
	print "## Configuring Module 5 ## "	
	print "Provide Genome Reference file"
	print "Press (Enter/Return) for Default"
	print "Default for demo is"+ work_dir+"/demo/annotation_files/lj_r25_demo.fa"
	print
	take=raw_input()
	if(take==''):
		genome_file=work_dir+'/demo/annotation_files/'+"lj_r25_demo.fa"
	else:
		genome_file=take
	
	######################################################
	### tasi_mapping_infile='' if mapped file is       ###
	### generated in module-2                          ###
	######################################################
	
	print "Provide a standard mapped file,i.e., .bed or .igv format"
	print "Press (Enter/Return) for Default"
	print "Default is the file generated by Module 2"
	print "For Demo, please provide file :"+work_dir+"/demo/annotation_files/sample.igv"
	print
	take=raw_input()
	if(take==''):
		tasi_mapping_infile=''
	else:
		tasi_mapping_infile=take
	######################################################
	### tasi_cluster_infile='' if the cluster file is  ###
	### same as generated in module -3 				   ###
	######################################################
	
	print "Provide a cluster file"
	print "Press (Enter/Return) for Default"
	print "Default is generated from module -3.1"
	print "For Demo, please provide file :"+work_dir+"/demo/annotation_files/sample.igv_cluster.bed"
	print
	take=raw_input()
	if(take==''):
		tasi_cluster_infile=''
	else:
		tasi_cluster_infile=take


##########################################################
####################### Module 6 #########################
### module adds annotations to small RNA profiles      ###
### store in the MySQL table						   ###
##########################################################


if(Add_annotation_to_mysql_database==True):			   ###	
	######################################################
	### provide annotation files, separated by commas  ###
	### these file could host genome, homologous genome###
	### repeat dataset, miRNA dataset or any other     ###
	### fatsa files, annotation coloumn is added and   ###
	### header represent name of the file			   ###
	######################################################
	print "## Configuring Module 6 ## "	
	print "Would like to keep default parameters mapping ?"
	print "Press (enter/return) for default [yes]"
	print "else press any other key for [No]"
	print 
	
	take=raw_input()
	if(take==''): 
		mismatches,max_maq_err,seedlen,no_of_cores = default_mapping_parameter()
	else:
		mismatches,max_maq_err,seedlen,no_of_cores = mapping_parameter()
	
	if (server==''):
		server,user,password,database=mysql()
	print "Provide annotation files seperated by commas(,)"
	print "Press (Enter/Return) for Default"
	print "Default is "+work_dir+"/demo/annotation_files/sample_reference.fa, "+work_dir+"/demo/annotation_files/miRBase.fa"
	print
	take=raw_input()
	if(take==''):
		annotation_file=work_dir+'/demo/annotation_files/'+"sample_reference.fa",work_dir+'/demo/annotation_files/'+'miRBase.fa'
	else:
		annotation_file=take
		
	infile=''
	######################################################
	### this module adds the annotation as exonic,     ###
	### intronic or intergenic depending upon the      ###
	### mapping position 							   ###
	######################################################
	
	print "Would you like to add genomic region annotation for each small RNA"
	print "Press (Enter/Return) for Default else anyother key"
	print "Default is False"
	print
	take=raw_input()
	if(take==''):
		add_genomic_region=False
	else:	
		add_genomic_region=True
	
	if(add_genomic_region==True):						   ###
		######################################################
		### if True provide following parameters       	   ###
		### gtf file is gene model and must be given       ###
		######################################################
		
		print "Provide a annotation (i.e.,.gtf,.gff) file"
		print "press (Enter/Return) for Demo"
		print "Default is "+work_dir+'/demo/annotation_files/'+"sample.gtf"
		print
		take=raw_input()
		if(take==''):
			gtf_file=work_dir+'/demo/annotation_files/'+"sample.gtf" 
		else:
			gtf_file=take


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
 
 
if(Pattern_plot==True):

		
	##########################################################
	### pattern_plot_file='' if same as produced in        ###
	### module-1, else provide a profile table             ###
	##########################################################
	print "## Configuring Module 7 ## "	
	print "Provide a profile table"
	print "Press (Enter/Return) for using output from Module 3"
	print "else provide full path of abundance table"
	print
	take=raw_input()
	if(take==''):
		pattern_plot_file=''
	else:
		pattern_plot_file=take
		
		
	print "have you already modified files 'libraries_reference' and 'compare_libraries' ? "
	print "libraries_reference contains reference for resepctive libraries, default is no refenrence"
	print "compare_libraries file contains information on which plot should be made, by default all vs all comparisons are plotted"
	
	print "press (Enter/Return) for default or first modify files and then press (Enter/Return)"
	print
	take=raw_input()
	
infile=''


##############################################################
######################## Module 8 ############################
### Pattern search 						                   ###
##############################################################

 
##############################################################
######################## Module 8.1 ##########################
### this module uses hidden morkov model to calcuate       ###
### frequency of next nucleotides given previous one       ###
### a morkov probability output file is generated          ###
##############################################################



if(Find_pattern==True):
	##########################################################
	### if True provide following parameters       	       ###
	### infile_for_pattern='' if same as produced in       ###
	### module-1, else provide a profile table             ###
	##########################################################
	print "## Configuring Module 8.1 ## "	
	print "Provide a profile table"
	print "Press (Enter/Return) for using output from Module 3"
	print "else provide full path of abundance table"
	print
	take=raw_input()
	if(take==''):
		infile_for_pattern=''
	else:
		infile_for_pattern=take
	##########################################################
	## ref_file_for_pattern is reference file as used before #
	## for example ref_file_for_pattern=$HOME/reference.fa  ##
	##########################################################
	
	print "Provide Genome Reference file"
	print "Press (Enter/Return) for Default"
	print "Default for demo is"+ work_dir+"/demo/annotation_files/lj_r25_demo.fa"
	print
	take=raw_input()
	if(take==''):
		ref_file_for_pattern=work_dir+'/demo/annotation_files/'+"lj_r25_demo.fa"
	else:
		ref_file_for_pattern=take


##############################################################
######################## Module 8.2 ##########################
### map all the sequences with a starting pattern          ###
##############################################################



if(Map_pattern==True):									   ###
	##########################################################
	### if True provide following parameters       	       ###
	### genome_file is reference file as used before       ###
	### for example genome_file=$HOME/reference.fa         ###
	##########################################################
	print "## Configuring Module 8.2 ## "	
	print "Provide Genome Reference file"
	print "Press (Enter/Return) for Default"
	print "Default for demo is"+ work_dir+"/demo/annotation_files/lj_r25_demo.fa"
	print
	take=raw_input()
	if(take==''):
		genome_file=work_dir+'/demo/annotation_files/'+"lj_r25_demo.fa"
	else:
		genome_file=take
	
	print "Would you like to keep default parameters mapping ?"
	print "Press (enter/return) for default [yes]"
	print "else press any other key for [No]"
	print 
	
	take=raw_input()
	if(take==''): 
		mismatches,max_maq_err,seedlen,no_of_cores = default_mapping_parameter()
	else:
		mismatches,max_maq_err,seedlen,no_of_cores = mapping_parameter()
	
	##########################################################
	### infile_for_pattern='' if same as produced in       ###
	### module-3, else provide a profile table             ###
	##########################################################
	
	print "Provide a profile table"
	print "Press (Enter/Return) for using output from Module 3"
	print "else provide full path of abundance table"
	print
	take=raw_input()
	if(take==''):
		infile_for_pattern_map=''
	else:
		infile_for_pattern_map=take
	
	###########################################################
	### pattern is the set of nucleotide which are present  ###
	### in the start of sequence                            ###
	###########################################################
	
	print "Provide a 5' nucleotide sequence"
	print "press (Enter/Return) for default"
	print "default is GAA"
	print
	take=raw_input()
	if(take==''):
		pattern="GAA"
	else:
		pattern=take


##############################################################
######################## Module 9  ###########################
### making mySQL queries 						           ###
##############################################################


if(MySQL_query==True):									   ###
	##########################################################
	### if True provide following parameters       	       ###
	### make a file with name 'query_infile' and write     ###
	### MySQL syntax query, it will be processed as output ###
	### will be send to text file                          ###
	##########################################################
	print "## Configuring Module 9 ## "	
	if(server==''):
		server,user,password,database=mysql()
	print "Please update the file named 'query_infile' to make a query"
	query_infile='query_infile'
	

