#!/bin/bash

######################################################################
#Which step to run? true("1") or false ("0") 
format_converter_enable=0
fam_computing=0
	hmmer_enable=0
	portho_enable=1
gecko_enable=0
results_parser_enable=0
plot_enable=1
######################################################################
actual_path=$PWD
path_where_your_scrips_are=/vol/cluster-data/castillo/scripts
## Path to .gbk genome files || input genomes
genome_directory_gbks=/vol/cluster-data/castillo/pezizomycetes_reduced_uptodate/
######################################################################
######## Format converter, the program "create_project"
create_project_path=/vol/cluster-data/castillo/ffgc_v0.3/bin/create_project
#old:create_project_path=/vol/cluster-data/castillo/ffgc_v0.2/build_ffgc/exe.linux-x86_64-2.7/create_project
#Folder name for the genomes
gos_file_folder_name=proyecto_danis_fungi_reduced
post_create_pro_path=$path_where_your_scrips_are/gostofastab
######################################################################
######## HMM stuff
## Path to HMM profiles to hmmscan them  
###you have to have hmmscan installed!!!###
hmm_script_path=$path_where_your_scrips_are/hmmcrap1
path_to_hmm_profile_folder=/vol/cluster-data/castillo/sacNOG_hmm
name_for_the_whole_hmm_file=aspergi
e_val_hmm=1e-5
post_hmm_path=$path_where_your_scrips_are/prepforgeck_egg
######################################################################
######## Proteinortho stuff
portho_script_path=$path_where_your_scrips_are/porthomystuff1
portho_path=/vol/cluster-data/castillo/proteinortho_v5.15/proteinortho5.pl
e_val_portho=1e-5
post_portho_path=$path_where_your_scrips_are/prepforgeck_portho
######################################################################
######## Gecko3 stuff
gecko_path=/homes/castillo/Desktop/castillo/Gecko3.1/bin/Gecko3
#Gecko input file name
cog_file_name_for_gecko=Danisfungicool6
#Gecko parameters in arrays
quorum=(7)
min_cluster_size=(3)
distance=(6 7)
######################################################################
######## Post Gecko3 
#must have python 2.7 installed
post_gecko_path=$path_where_your_scrips_are/gckresultsanalysis.py
######################################################################
######## Plot 
###must have R installed with ggplot2
#plot_path=$path_where_your_scrips_are/plotstuffnew.R
#plot_path_2=$path_where_your_scrips_are/newplotforclusts.R
plot_path_3=/homes/castillo/Desktop/castillo/PycharmProjects/plot2dhisto/plot2dh.py
######################################################################
######## Reported clusters data 
reported_clust_path=/homes/castillo/Desktop/castillo/reportedclustersupdt_verified/
#reported_clust_path=/homes/castillo/Desktop/castillo/reportedclusters2017_fit/
######################################################################
######################################################################

########### Create project
if [ $format_converter_enable -eq 1 ]
then 
	#get the .gbk file paths from the folder genome_directory_gbks
	files=$(ls -lrt -d -1 $genome_directory_gbks/{*,.*} | awk '{print $9}' | rev | grep -oP 'kbg\..*' | rev | tr '\n' ' ')

	#generate the gos from the gbk in the folder (genome_directory_gbks)
	if [ ! -d "$gos_file_folder_name" ]
	then
		$create_project_path -f GBK -p ./$gos_file_folder_name $files 
	else
		echo "The project directory already exists"
	fi
	
	# remove the error of having ge instead of gb
	cd ./$gos_file_folder_name/genomes
	sed -i -- 's/^>ge|/>gb|/g' *.gos
	cd ..
	cd .. 	

	#convert gos header into tabular form  
	$post_create_pro_path ./$gos_file_folder_name/genomes
fi

#########################################################
########### Family annotation 

if [ $fam_computing -eq 1 ]
then 

	if [ $hmmer_enable  -eq 1 ] && [ $portho_enable  -eq 1 ]
	then 
		echo "Can't do both families at the same time"		
		exit 1
	fi

	if [ $hmmer_enable -eq 0 ] && [  $portho_enable -eq 0 ]
	then
		echo "You are not computing families with hmmer nor portho"
		exit 1
	fi

	########### HMMER time! get it? 

	if [ $hmmer_enable -eq 1 ]
	then 
		$hmm_script_path $actual_path/$gos_file_folder_name/genomes $path_to_hmm_profile_folder $name_for_the_whole_hmm_file $e_val_hmm
		$post_hmm_path $actual_path/$gos_file_folder_name/genomes $cog_file_name_for_gecko 
	fi

	########### Proteinortho

	if [ $portho_enable -eq 1 ]
	then 
		$portho_script_path $actual_path/$gos_file_folder_name/genomes $portho_path $e_val_portho
		$post_portho_path $actual_path/$gos_file_folder_name/genomes $cog_file_name_for_gecko "$actual_path/$gos_file_folder_name/genomes/familiesportho"

	fi
fi
#########################################################
########### Gecko

if [ $gecko_enable  -eq 1 ]
then 

echo "Quorum values:"
	echo ${quorum[*]}
	echo "Minimum cluster size:"
	echo ${min_cluster_size[*]}
	echo "Distance/Indels:"
	echo ${distance[*]}


	qu=0
	mcs=0
	di=0
	for qu in "${!quorum[@]}"
	do	
		for mcs in "${!min_cluster_size[@]}"
		do	
			for di in "${!distance[@]}"
			do
				qu_a=${quorum[qu]}
				mcs_a=${min_cluster_size[mcs]}
				di_a=${distance[di]}

				if [ -f "stats-q$qu_a-mc$mcs_a-d$di_a.gckstats" ]
				then 
					echo "stats-q$qu_a-mc$mcs_a-d$di_a.gckstats already exists" 
				else
			
					echo "submitting:     quorum:$qu_a / min clusters:$mcs_a / distance:$di_a" 
					qsub -cwd -o "GeckoRslt$cog_file_name_for_geck-q$qu_a-mc$mcs_a-d$di_a.log" -V -j y -pe multislot 4 -l vf=12G -l idle=1 -P fair_share -b y $gecko_path -in "./$gos_file_folder_name/genomes/$cog_file_name_for_gecko.cog" -out "GeckoRslt$cog_file_name_for_geck-q$qu_a-mc$mcs_a-d$di_a.gck" -q $qu_a -s $mcs_a -d $di_a --resultOutput "clusterStatistics" showFiltered "stats-q$qu_a-mc$mcs_a-d$di_a.gckstats"

				fi
			done
		done
	done


fi


if [ $hmmer_enable  -eq 1 ] && [ $portho_enable  -eq 1 ]
then
		echo "cannot plot eggnog and proteinortho in the same project"
		exit 1
fi

if [ $hmmer_enable -eq 0 ] && [  $portho_enable -eq 0 ]
then
		echo "You are not plotting anything, enable hmmer or portho"
		exit 1
fi

if  [ $hmmer_enable  -eq 1 ]
then
		if [ $results_parser_enable -eq 1 ]
		then
			if [ -f "gckoeggout_summary" ]
			then 
				echo "gckoeggout_summary already exists"
			else
				python $post_gecko_path "$PWD" $reported_clust_path gckoeggout 
			fi
		fi

		if  [ $plot_enable  -eq 1 ]
		then
			#$plot_path gckoeggout_summary
			python $plot_path_3 gckoeggout2_summary plotty
		fi
fi
if [ $portho_enable  -eq 1 ]
then
		if [ $results_parser_enable -eq 1 ]
		then
			if [ -f "gckoporthoout_summary" ]
			then 
				echo "gckoporthoout_summary already exists"
			else
				python $post_gecko_path "$PWD" $reported_clust_path gckoporthoout 
			fi
		fi

		if  [ $plot_enable  -eq 1 ]
		then
			#Rscript $plot_path gckoporthoout_summary
			#Rscript $plot_path_2 gckoporthoout_summary
			python $plot_path_3 gckoeggout2_summary plotty
		fi
fi








