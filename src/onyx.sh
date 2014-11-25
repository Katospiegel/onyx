#!/bin/bash

##########################################################################################

#Onyx is a tools for know  miRNA measure expression levels in HTS data.
#
#Author: Carlos Vivar
#
#Parameters:
#<config.txt file> - The file where all the path and parameters are defined.

#-----------------------------PARAMETERS  AND PATHS------------------------------------#

STATIME=$(date +%s)

if [ $# != 1 ];
then
	echo "Insert the <config.txt file> parameter correctly"
else
	Config=$1
		
	#Path and parameters.
	Pro=$( cat $Config | grep -w '^PROJECT_PATH' | cut -d '=' -f2)
	Raw=$( cat $Config | grep -w '^RAW_DATA' | cut -d '=' -f2)
	Ref=$( cat $Config | grep -w '^REF_PATH' | cut -d '=' -f2)

	#Referencies
	Genome=$( cat $Config | grep -w '^GENOME_BOWTIE_NAME' | cut -d '=' -f2)
	Genome_path=$Ref/$Genome
	

	miRBase=$( cat $Config | grep -w '^miRBase_gff3_NAME' | cut -d '=' -f2)
	miRBase_path=$Ref/$miRBase

	#Folders creation.
	mkdir -p $Pro/result/

	mkdir -p $Pro/result/adtrim/
	mkdir -p $Pro/result/count/
	mkdir -p $Pro/result/map/
	mkdir -p $Pro/result/expression/

	mkdir -p $Pro/doc/
	mkdir -p $Pro/doc/prefastqc/
	mkdir -p $Pro/doc/postfastqc/
	mkdir -p $Pro/doc/log

	PRQC=$Pro/doc/prefastqc/
	ADTRIM=$Pro/result/adtrim/
	PSQC=$Pro/doc/postfastqc/
	MAP=$Pro/result/map/
	COUNT=$Pro/result/count/
	EXPRESSION=$Pro/result/expression/
	LOG=$Pro/doc/log

	#Tools Pathways
	FASTQC=$( cat $Config | grep -w '^FASTQC' | cut -d '=' -f2)
	CUTADAPT=$( cat $Config | grep -w '^CUTADAPT' | cut -d '=' -f2)
	BOWTIE=$( cat $Config | grep -w '^BOWTIE' | cut -d '=' -f2)
	HTseq=$( cat $Config | grep -w '^HTseq' | cut -d '=' -f2)
	R=$( cat $Config | grep -w '^R' | cut -d '=' -f2)

	#Tool Parameters.
	CUTADAPT_PAR=$( cat $Config | grep -w '^CUTADAPT_PAR' | cut -d '=' -f2)
	BOWTIE_PAR=$( cat $Config | grep -w '^BOWTIE_PAR' | cut -d '=' -f2)
	HTseq_PAR=$( cat $Config | grep -w '^HTseq_PAR' | cut -d '>' -f2)
	R_PAR=$( cat $Config | grep -w '^R_PAR' | cut -d '=' -f2)
	
	
	#General Information.
	PRO_NAME=$( cat $Config | grep -w '^PROJECT_NAME' | cut -d '=' -f2)
	SAMPLES=$( cat $Config | grep -w '^SAMPLES_NAME' | cut -d '=' -f2)
	GROUP=$( cat $Config | grep -w '^GROUPS' | cut -d '=' -f2)
	PAIRS=$( cat $Config | grep -w '^PAIRS' | cut -d '=' -f2)
	PAIR_SWITCH=$( cat $Config | grep -w '^PAIR_SWITCH' | cut -d '=' -f2)
	NUM_SAMPLES=$( echo $SAMPLES | awk -F\; '{print NF}' )
	
	#Optimization
	OPT=$( cat $Config | grep -w '^OPTIMIZATION' | cut -d '=' -f2)
	NPROC=$( nproc )
	WPROC=$( expr $NPROC / 2 )
	PROC_SAM=$( expr $WPROC / $NUM_SAMPLES )

#--------------------------------------PIPELINE-----------------------------------#

	echo 'Onyx miRNA expression Analysis running with '$Genome' mapping and '$miRBase' miRBase annotation for '$PRO_NAME' data'

#------------> LOG_FILE_CREATION
	for sample in $Raw/*fastq
	do
	
		NAME=$(basename "$sample" | cut -d '.' -f1)
		echo 'Onyx miRNA expression Analysis running with '$Genome' mapping and '$miRBase' miRBase annotation for ;
'$PRO_NAME' data' >> $LOG/$NAME.txt 2>&1

	done

#------------> PRE_FASTQC

	echo 'FASTQ QUALITY CONTROL OF THE RAW DATA'
	
	
	if [ $OPT == "YES" ];
	then

		opt_pre_fastqc(){
		local run=$1
		NAME=$(basename "$sample" | cut -d '.' -f1)	
		$FASTQC -o $PRQC -f fastq $sample 
		}
		for sample in $Raw/*fastq; do opt_pre_fastqc $sample & done;

		wait	

	else

		NAME=$(basename "$sample" | cut -d '.' -f1)
		$FASTQC -o $PRQC -f fastq $Raw/*fastq 
	
	fi	

#------------> TRIMMING
	
	echo 'TRIMMING' 
	
	if [ $OPT == "YES" ];
	then

		opt_trimming(){
		local run=$1
		NAME=$(basename "$file" | cut -d '.' -f1)
		echo 'TRIMMING' >> $LOG/$NAME.txt 2>&1		
		$CUTADAPT $CUTADAPT_PAR -o $ADTRIM/$NAME.fastq -f fastq $file  >> $LOG/$NAME.txt 2>&1
		}
		for file in $Raw/*fastq; do opt_trimming $file & done;
	
		wait

	else

		for file in $Raw/*fastq
		do
	
		NAME=$(basename "$file" | cut -d '.' -f1)
		echo 'TRIMMING' >> $LOG/$NAME.txt 2>&1
		$CUTADAPT $CUTADAPT_PAR -o $ADTRIM/$NAME.fastq -f fastq $file >> $LOG/$NAME.txt 2>&1

		done
	
	fi

#------------> POST_FASTQC

	echo 'FASTQ QUALITY CONTROL OF THE TRIMMED READS'

	if [ $OPT == "YES" ];
	then

		opt_post_fastqc(){
		local run=$1
		NAME=$(basename "$sample" | cut -d '.' -f1)
		$FASTQC -o $PSQC -f fastq $sample  
		}
		
		for sample in $ADTRIM/*fastq; do opt_post_fastqc $sample & done;
	
		wait	

	else
	
		NAME=$(basename "$sample" | cut -d '.' -f1)
		$FASTQC -o $PSQC -f fastq $ADTRIM/*fastq 

	fi

#------------> MAPPING

	echo 'MAPPING AGAINST '$Genome'' 

	for file in $ADTRIM/*fastq
	do
	
		NAME=$(basename "$file" | cut -d '.' -f1)
		echo 'MAPPING AGAINST '$Genome'' >> $LOG/$NAME.txt 2>&1
		#$BOWTIE $BOWTIE_PAR -x $Genome_path -U $file -S $MAP/$NAME.sam >> $LOG/$NAME.txt 2>&1
		$BOWTIE $BOWTIE_PAR $Genome_path $file $MAP/$NAME.sam >> $LOG/$NAME.txt 2>&1
	
	done

#-----------> ANNOTATION

	echo 'ANNOTATION AGAINST '$miRBase'' >> $LOG/$NAME.txt 2>&1

	if [ $OPT == "YES" ];
	then

		opt_annotation(){
		local run=$1
		NAME=$(basename "$file" | cut -d '.' -f1)
		echo 'ANNOTATION AGAINST '$miRBase'' >> $LOG/$NAME.txt 2>&1	
		$HTseq $HTseq_PAR -o $COUNT/$NAME.sam $file $miRBase_path  > $COUNT/$NAME.csv
		tail -n 5 $COUNT/$NAME.csv >> $LOG/$NAME.txt 2>&1
		}

		for file in $MAP/*sam; do opt_annotation $file & done;
	
		wait

	else

		for file in $MAP/*sam
		do
	
			NAME=$(basename "$file" | cut -d '.' -f1)
			echo 'ANNOTATION AGAINST '$miRBase'' >> $LOG/$NAME.txt 2>&1
			$HTseq 	$HTseq_PAR -o $COUNT/$NAME.sam $file $miRBase_path > $COUNT/$NAME.csv
			tail -n 5 $COUNT/$NAME.csv >> $LOG/$NAME.txt 2>&1

		done

	fi

#-----------> RAW EXPRESSION TABLE CREATION
	
	TABLE_1=$(echo $SAMPLES | sed 's/;/ /g')
	echo $SAMPLES
	for sample in $TABLE_1
	do 
		file=$COUNT/$sample.csv
		echo $file >> $COUNT/samples_temp.txt
	done
	
	NUM_MAX=$(expr $NUM_SAMPLES + $NUM_SAMPLES )
	
	for num in $(seq 1 $NUM_MAX)
	do 
		[ $((num % 2)) -eq 0 ] && echo \$$num >> $Pro/temp
	done
	
	PAIR_NUM=$(nawk '$1=$1' RS= $Pro/temp | sed 's/ /,/g')

	TABLE_2=$(nawk '$1=$1' RS= $COUNT/samples_temp.txt)
	paste $TABLE_2 | awk '{print $1,'$PAIR_NUM'}' > $COUNT/raw_expression_temp.csv

	head -n -5 $COUNT/raw_expression_temp.csv > $EXPRESSION/raw_expression.csv
	
	rm $COUNT/samples_temp.txt
	rm $COUNT/raw_expression_temp.csv
	rm $Pro/temp
		
	
#-----------> DIFFERENTIAL EXPRESSION ANALYSIS

	echo 'DIFFERENTIAL EXPRESSION ANALYSIS' 

	if [ $PAIR_SWITCH == "YES" ];
	then	

		Rscript $Pro/bin/differential_ex_pair.R $EXPRESSION/raw_expression.csv $EXPRESSION $SAMPLES $GROUP $PAIRS $PRO_NAME

	else

		Rscript $Pro/bin/differential_ex_nonpair.R $EXPRESSION/raw_expression.csv $EXPRESSION $SAMPLES $GROUP  $PRO_NAME

	fi

#-----------> SUMMARY

	for file in $LOG/*txt
	do
		NAME=$(basename "$file" | cut -d '.' -f1)
		echo $NAME >> $LOG/$NAME.csv
		cat $file | grep -w 'Processed reads:' | awk '{ print $3 }' >> $LOG/$NAME.csv
		cat $file | grep -w 'Processed bases:' | awk '{ print $3, $4 }' >> $LOG/$NAME.csv
		cat $file | grep -w 'Trimmed reads:' | awk '{ print $3, $4 }' >> $LOG/$NAME.csv
		cat $file | grep -w 'Quality-trimmed:' | awk '{ print $2, $3, $4 }' >> $LOG/$NAME.csv
		cat $file | grep -w 'Too short reads:' | awk '{ print  $4, $5 }' >> $LOG/$NAME.csv

		cat $file | grep -w 'reads processed:' | awk '{ print $4 }' >> $LOG/$NAME.csv
		cat $file | grep -w 'reads with at least' | awk '{ print $9, $10 }' >> $LOG/$NAME.csv
		cat $file | grep -w 'failed to align' | awk '{ print $7, $8 }' >> $LOG/$NAME.csv
		cat $file | grep -w 'Reported' | awk '{ print $2 }' >> $LOG/$NAME.csv

		cat $file | grep -w 'no_feature' | awk '{ print $2 }' >> $LOG/$NAME.csv
		cat $file | grep -w 'ambiguous' | awk '{ print $2 }' >> $LOG/$NAME.csv
		cat $file | grep -w 'not_aligned' | awk '{ print $2 }' >> $LOG/$NAME.csv

	done

	#Final table
	
	echo 'Samples' >> $LOG/summary.csv
	echo 'Trim_Processed_reads:'  >> $LOG/summary.csv
	echo 'Trim_Processed_bases:'  >> $LOG/summary.csv
	echo 'Trim_Trimmed_reads:'  >> $LOG/summary.csv
	echo 'Trim_Quality-trimmed:'  >> $LOG/summary.csv
	echo 'Trim_Too short_reads:'  >> $LOG/summary.csv

	echo 'Map_reads_processed:'  >> $LOG/summary.csv
	echo 'Map_reads_with_at_least:'  >> $LOG/summary.csv
	echo 'Map_failed_to_align:'  >> $LOG/summary.csv
	echo 'Maap_Reported:'  >> $LOG/summary.csv

	echo 'Ann_no_feature:'  >> $LOG/summary.csv
	echo 'Ann_ambiguous:'  >> $LOG/summary.csv
	echo 'Ann_not_aligned:' >> $LOG/summary.csv
	
	
	for name in $TABLE_1
	do
		file=$LOG/$name.csv
		echo $file >> $LOG/samples_temp.txt
	done

	TABLE_3=$(nawk '$1=$1' RS= $LOG/samples_temp.txt)
	paste $LOG/summary.csv $TABLE_3 | awk '{print }' > $LOG/final_summary.csv
	
	rm $LOG/samples_temp.txt
	rm $LOG/summary.csv

	for name in $TABLE_1
	do
		rm $LOG/$name.csv
	done
		

#-----------> HTML CREATION

	cp $Pro/bin/MainDocument.html $Pro/MainDocument.html

	#perl $Pro/bin/html_report.pl $config $Pro/MainDocument.html $LOG/final_summary.csv $EXPRESSION/dif_exp_sum.xls
	#cp onyx_pipeline.jpg $Pro/onyx_pipeline.jpg

#-----------> ENDING

	ENDTIME=$(date +%s)

	echo 'Onyx EXPRESSION ANALYSIS COMPLETE'

	echo 'Task time: '$(($ENDTIME - $STATIME))' s'

	
fi
	


