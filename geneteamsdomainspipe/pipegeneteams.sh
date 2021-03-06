#!/bin/bash
######################################################################
#Which step to run? true("1") or false ("0") 
making_graph_files=0
gene_teams_qsub_jobarray=0
result_summary=0
plot_enable=0
######################################################################
actual_path=$(pwd)
######################################################################
######## making graph files
#use the file computed for gecko, domains or raw
cog_to_graph_gml=/homes/castillo/Desktop/castillo/PycharmProjects/cogfiletograph/cogtograf.py
cog_file_for_gecko=/homes/castillo/Desktop/castillo/A_geckoupdgenomes_porto_newffgc/proyecto_danis_fungi_reduced/genomes/Danisfungicool6.cog
######################################################################
######## Gene teams
gene_teams_script=/homes/castillo/Desktop/castillo/GraphTeams/scripts/graph_teams_mod01.py
#job array
script_with_gt_job_array=/homes/castillo/Desktop/castillo/scripts/gtjobarray
######################################################################
#####create parameter file
> parametersgeneteams
echo "quorum_combinations=$actual_path/filepathlistsquorum" >> parametersgeneteams
echo "distance=0 1 2 3 5 10" >> parametersgeneteams
## quorums already ran: 3,
echo "quorum=4" >> parametersgeneteams

######################################################################
#####result analysis summary maker
summary_output=gtoutput
summ_maker=/homes/castillo/Desktop/castillo/PycharmProjects/aftergeneteamsnodom/aftergtnd.py

######## Plot 
#must have R installed with ggplot2 
#plot_path=$path_where_your_scrips_are/plotstuffnew.R
#plot_path_2=$path_where_your_scrips_are/newplotforclusts.R
plot_path_3=/homes/castillo/Desktop/castillo/PycharmProjects/plot2dhisto/plot2dh.py
######################################################################
######## Reported clusters data 
#you have to update the clusters for this shit
#reported_clust_path=/homes/castillo/Desktop/castillo/reportedclustersupdt_verified/
reported_clust_path=/homes/castillo/Desktop/castillo/reportedclusters2017_fit/
#reported_clust_path=/homes/castillo/Desktop/castillo/reportedclusts2017_fit_fordatasetlike2017/
######################################################################
######################################################################
######################################################################
########### make graph files
if [ $making_graph_files -eq 1 ]
then 
	python $cog_to_graph_gml $cog_file_for_gecko "raw" $actual_path
fi

#########################################################
########### geneteams submition as job array 

mkdir geneteamoutputs

if [ $gene_teams_qsub_jobarray  -eq 1 ]
then 
	quorum_combinations=$(grep -oP 'quorum_combinations=\K.*' parametersgeneteams)
	distance_str=$(grep -oP 'distance=\K.*' parametersgeneteams)
	distance=( $distance_str )
	quorum_str=$(grep -oP 'quorum=\K.*' parametersgeneteams)
	quorum=( $quorum_str )
	echo "Distance/Indels:"
	echo ${distance[*]}
	echo "Quorum:"
	echo ${quorum[*]}

	files=$(ls -lrt -d -1 $actual_path/{*,.*} | awk '{print $9}' | rev | grep -oP 'lmg\..*' | rev | tr '\n' ' ')

	qu=0
	di=0
	for qu in "${!quorum[@]}"
	do
		for di in "${!distance[@]}"
		do
			qu_a=${quorum[qu]}
			di_a=${distance[di]}
			IFS=$'\n'
			quorum_combi_entries=( $(grep -oP "^$qu_a\t\K.*" $quorum_combinations) ) 
			for combi in ${!quorum_combi_entries[@]}
			do
				qcombi=${quorum_combi_entries[combi]}
				submittions+=("python $gene_teams_script -d $di_a $qcombi ")
			done	
		done
	done

	IFS=$' \t\n'

	amount_o_runs=${#submittions[@]}

	echo "el amount de runs es: $amount_o_runs"

	qsub -t 1-$amount_o_runs -o "./outputshido.log" -cwd -V -j y -pe multislot 2 -l vf=8G -l idle=1 -P fair_share -b y $script_with_gt_job_array

fi

if [ $result_summary  -eq 1 ]
then 

	<script path>	<graph teams results folder path> <cluster path> <outputname>
	$summ_maker	$actual_path/geneteamoutputs	$reported_clust_path	$summary_output

fi


if  [ $plot_enable  -eq 1 ]
then
	#Rscript $plot_path $summary_output
	#Rscript $plot_path_2 $summary_output
	python $plot_path_3 gckoeggout2_summary plotty
fi

