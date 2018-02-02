#!/bin/bash
######################################################################
#Which step to run? true("1") or false ("0") 
format_converter_enable=0 
domain_computing=0 
post_domain_comp=0
promiscous_removal=0
gecko_enable=1
#Don't use the non-parallel if you want to finish this year
results_parser_enable=0
#Use this one for faster
results_parser_parallel=0
after_parallel=0
plot_enable=0
######################################################################
actual_path=$(pwd)
path_where_your_scrips_are=/vol/cluster-data/castillo/scripts
## Path to .gbk genome files || input genomes
genome_directory_gbks=/vol/cluster-data/castillo/pezizomycetes_reduced_uptodate/
######################################################################
######## Format converter, the program "create_project"
create_project_path=/homes/castillo/Desktop/castillo/ffgc_v0.3.src/create_project.py
#old:create_project_path=/vol/cluster-data/castillo/ffgc_v0.2/build_ffgc/exe.linux-x86_64-2.7/create_project
#Folder name for the genomes
gos_file_folder_name=proyecto_danis_fungi_reduced
post_create_pro_path=$path_where_your_scrips_are/gostofastab
######################################################################
######################################################################
######## Domain stuff
interpro_scripto_path=$path_where_your_scrips_are/preinterscan 
interproscan_path=/vol/cluster-data/castillo/interproscan-5.23-62.0/interproscan.sh
######## Pre gecko domains
pregeckodomains_path=/homes/castillo/Desktop/castillo/PycharmProjects/domcogfilemaker/domcogmkr.py
###annotation type: pfam or ipr
annotation_type="ipr"
####################################
#######After, promiscuous protein domain families removal 
path_to_script_to_remove_promis=/homes/castillo/Desktop/castillo/PycharmProjects/promgenerem/prmgeneremvl.py
cog_file_name_for_gecko="Danisfungicool6"
##Promiscous removal suffix
promis_suffix="_prom_remvd"
percentil=99
######################################################################
######## Gecko3 stuff
#####normal gecko#####
#gecko_path=/vol/cluster-data/castillo/Gecko3.1/bin/Gecko3
#####a lot of memory gecko####
gecko_path=/vol/cluster-data/castillo/Gecko3.1/bin/Gecko3_32g
#Gecko input file name
cog_file_name_for_gecko_nopromis="Danisfungicool6$promis_suffix"
#Gecko parameters in arrays
quorum=(3 4 5 6 7 8 9 11)
min_cluster_size=(3)
distance=(0 1 2 3 4 5 6 7 8)
######################################################################
######## Post Gecko3 
#must have python 2.7 installed
post_gecko_path=$path_where_your_scrips_are/gckdomtogene.py
######################################################################
######## Post Gecko3 parallelized
#must have python 2.7 installed
post_gecko_parallel_submittion=/homes/castillo/Desktop/castillo/PycharmProjects/geckdomainstogenesubmitable/submitpostgeckdomains
post_gecko_parallel=/homes/castillo/Desktop/castillo/PycharmProjects/geckdomainstogenesubmitable/geckdomtogenesubmi.py
after_parallel_concat=/homes/castillo/Desktop/castillo/scripts/concatdomtogenegeckoout

######################################################################
######## Plot 
#must have R installed with ggplot2 
plot_path=$path_where_your_scrips_are/plotstuffnew.R
plot_path_2=$path_where_your_scrips_are/newplotforclusts.R
plot_path_3=/homes/castillo/Desktop/castillo/PycharmProjects/plot2dhisto/plot2dh.py
######################################################################
######## Reported clusters data 
#old clusts reported_clust_path=/homes/castillo/Desktop/castillo/reportedclustersupdt_verified/
#new clusts reported_clust_path=/homes/castillo/Desktop/castillo/reportedclusters2017_fit/
reported_clust_path=/homes/castillo/Desktop/castillo/reportedclustersupdt_verified/
######################################################################
######################################################################

########### Create project
if [ $format_converter_enable -eq 1 ]
then 
	#get the .gbk file paths from the folder genome_directory_gbks
	files=$(ls -lrt -d -1 $genome_directory_gbks/{*,.*} | awk '{print $9}' | rev | grep -oP 'kbg\..*' | rev | tr '\n' ' ')
	for elemento in $files
	do 
		echo "$elemento" >> genome_key_names
	done
	#generate the gos from the gbk in the folder (genome_directory_gbks)
	if [ ! -d "$gos_file_folder_name" ]
	then
		python $create_project_path -f GBK -p ./$gos_file_folder_name $files 
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
########### domain annotation 

if [ $domain_computing -eq 1 ]
then 
	cd ./$gos_file_folder_name/genomes
	$interpro_scripto_path ./$gos_file_folder_name/genomes $interproscan_path

fi

if [ $post_domain_comp -eq 1 ]
then 	

	cd ./$gos_file_folder_name/genomes
	python $pregeckodomains_path "$actual_path/$gos_file_folder_name/genomes" $cog_file_name_for_gecko $annotation_type
	
fi


if [ $promiscous_removal -eq 1 ]
then
##remove promiscous
	cd ./$gos_file_folder_name/genomes
	python $path_to_script_to_remove_promis "$actual_path/$gos_file_folder_name/genomes/$cog_file_name_for_gecko.cog" $percentil $promis_suffix
fi
 
#########################################################
########### Gecko



if [ $gecko_enable  -eq 1 ]
then 

cd $actual_path
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
					qsub -cwd -o "GeckoRslt$cog_file_name_for_geck-q$qu_a-mc$mcs_a-d$di_a.log" -N "geckodom$annotation_type-Q$qu_a-D$di_a" -V -j y -pe multislot 4 -l vf=16G -l idle=1 -P fair_share -b y $gecko_path -in "./$gos_file_folder_name/genomes/$cog_file_name_for_gecko_nopromis.cog" -out "GeckoRslt$cog_file_name_for_geck-q$qu_a-mc$mcs_a-d$di_a.gck" -q $qu_a -s $mcs_a -d $di_a --resultOutput "clusterStatistics" showFiltered "stats-q$qu_a-mc$mcs_a-d$di_a.gckstats"

				fi
			done
		done
	done


fi


if [ $results_parser_enable -eq 1 ]
			then
				cd $actual_path
				if [ -f "gckoeggout_summary" ]
				then 
					echo "gckoeggout_summary already exists"
				else
					python $post_gecko_path "$PWD" $reported_clust_path "./$gos_file_folder_name/genomes/" gckoeggout2 
				fi
			fi

if [ $results_parser_parallel -eq 1 ]
			then
				cd $actual_path
				if [ -f "gckoeggout2_summary" ]
				then 
					echo "gckoeggout2_summary already exists"
				else

					$post_gecko_parallel_submittion "./$gos_file_folder_name/genomes/" $reported_clust_path $post_gecko_parallel gckoeggout2
					
				fi
			fi

if [ $after_parallel -eq 1 ]
then
	cd $actual_path
	$after_parallel_concat

fi

if [ $plot_enable -eq 1  ] 
then 
	cd $actual_path

	#Rscript $plot_path gckoeggout2_summary
	#Rscript $plot_path_2 gckoeggout2_summary
	python $plot_path_3 gckoeggout2_summary plotty

fi





