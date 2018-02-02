#!/bin/bash

######################################################################
#Which step to run? true("1") or false ("0") 
format_converter_enable=0
dom_fam_computing=0
domains_to_sim_file=0
ffgc_macsi_enable=1
results_parser_enable=0
plot_enable=0
######################################################################
actual_path=$(pwd)
#/vol/cluster-data/castillo/scripts/newplotforclusts.R
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
######## Domain families stuff 
pre_interpro_step_path="$path_where_your_scrips_are/preinterscan"
interproscan_path=/vol/cluster-data/castillo/interproscan-5.23-62.0/interproscan.sh
######################################################################
######## Domains to sim stuff 
path_dom_to_sim_conv=/homes/castillo/Desktop/castillo/PycharmProjects/domtosimnopromis/domtosimnopromis.py
#path_dom_to_sim_conv=/vol/cluster-data/castillo/PycharmProjects/pdomtosim/pdomtosim.py
percentil_val=98
pfam_ipr=ipr
######################################################################
######## ffgc macsi stuff
macsi_path=/homes/castillo/Desktop/castillo/ffgc_v0.3.src/macsi.py
#parameters in arrays
quorum=(3 4 5 6)
min_cluster_size=(3)
distance=(0 1 2 3 4 5)
#for more than 5 distance, increase the memory in the submition

######################################################################
######## Post ffgc
weak_common_intervals_folder=/homes/castillo/Desktop/castillo/A_ffgcdom_nopromis/proyecto_danis_fungi_reduced/genomes/weak_common_intervals/
path_ffcg_out=/homes/castillo/Desktop/castillo/PycharmProjects/ffgcoutanalysiswerror/ffgcoutwerr.py
#path_ffcg_out=/vol/cluster-data/castillo/PycharmProjects/ffgcoutanalysis/ffgcout.py
output_summary_name="ffgcdomout_summary"
######################################################################
######## Plot 
#must have R installed with ggplot2 
#plot_path=$path_where_your_scrips_are/plotstuffnew.R
#plot_path_2=$path_where_your_scrips_are/newplotforclusts.R
plot_path_3=/homes/castillo/Desktop/castillo/PycharmProjects/plot2dhisto/plot2dh.py
######################################################################
######## Reported clusters data 
reported_clust_path=/homes/castillo/Desktop/castillo/reportedclusters2017_fit/
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
		

	#convert gos header into tabular form  
	$post_create_pro_path $actual_path/$gos_file_folder_name/genomes
fi

#########################################################
########### domain family annotation 

if [ $dom_fam_computing -eq 1 ]
then 
	cd "$actual_path/$gos_file_folder_name/genomes"

	"$path_where_your_scrips_are/preinterscan" "$actual_path/$gos_file_folder_name/genomes" $interproscan_path

fi



#########################################################
########### domains to sim file conversion

if [ $domains_to_sim_file -eq 1 ]
then 
	cd $actual_path/$gos_file_folder_name/genomes
	mkdir sim_files
	cd ./sim_files
 
	python $path_dom_to_sim_conv $actual_path/$gos_file_folder_name/genomes $percentil_val $pfam_ipr

	echo "This sim files were made with the $percentil_val percentile of the $pfam_ipr domain families computed with $interproscan_path" > details.log

fi


#########################################################
########### FFGC stuff

if [ $ffgc_macsi_enable  -eq 1 ]
then 


	cd $actual_path/$gos_file_folder_name/genomes
	mkdir weak_common_intervals
	cd ./weak_common_intervals

files_sim=$(ls -lrt -d -1 $actual_path/$gos_file_folder_name/genomes/sim_files/{*,.*} | awk '{print $9}' | rev | grep -oP 'mis\..*' | rev | tr '\n' ' ')

echo "Quorum values:"
	echo ${quorum[*]}
	echo "Minimum cluster size:"
	echo ${min_cluster_size[*]}
	echo "Distance/Indels:"
	echo ${distance[*]}
	
	#safecopy of yaml file, contains the date and a random number between 1 and 100
	
	#cp config.yaml "config_backupcopy$(date +"%Y%m%dT%H%M")$((1 + RANDOM % 100)).yaml" 

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

				
						
				qsub -cwd -o "ffgcres-q$qu_a-mc$mcs_a-d$di_a.log" -N "ffgcdomipr" -V -j y -pe multislot 8 -l vf=20G -l idle=1 -P fair_share -b y python $macsi_path -q $qu_a -d $di_a -m $mcs_a $files_sim

			done
		done
	done

fi




if [ $results_parser_enable -eq 1 ]
then

	python $path_ffcg_out $weak_common_intervals_folder $reported_clust_path $output_summary_name
		
fi
		

if  [ $plot_enable  -eq 1 ]
then
	#Rscript $plot_path $output_summary_name
	#Rscript $plot_path_2 $output_summary_name
	python $plot_path_3 gckoeggout2_summary plotty
fi








