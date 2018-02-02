import os
import sys
from copy import deepcopy
import numpy

def deleteContent(fName):
    with open(fName, "w"):
        pass

def make_sequence(a,b,increment):
    list=[]
    for i in range(a,b+1,increment):
        list.append(i)
    return list

def get_genome_key(file_name_whole):
    genome_key_dict = {}
    with open(file_name_whole[0], 'r') as file_chido:
        contador=0
        flag1 = False
        for line in file_chido:
            if flag1 == True:
                genome_key_dict[int(line.strip("G"))]=contador
                genome_key_file_2.write("%s is the genome number %s\n" % (str(line), str(contador)))

                contador = contador + 1
                flag1 = False
            if line == "<genome>\n":
                flag1 = True
            if line == "</genomes>\n":
                break

    #print genome_key_dict
    return genome_key_dict



def jaccardindex(list_a,list_b):
    set_a = set(list_a)
    set_b = set(list_b)
    union_chida = set_a.union(set_b)
    interseccion_chida = set_a.intersection(set_b)
    jaccard_index = float(len(interseccion_chida)) / float(len(union_chida))
    return jaccard_index

def found_clust_func(list_a,list_b,p):
    set_a = set(list_a)
    set_b = set(list_b)
    interseccion_chida = set_a.intersection(set_b)
    set_a_length = float(len(set_a))
    intersec_len = float(len(interseccion_chida))
    inter_over_length = intersec_len/set_a_length
    length_percentage = float(p/100.0)

    if inter_over_length >= length_percentage:
        return 1
    else:
        return 0

if __name__ == '__main__':

    if len(sys.argv) == 4:
        target_directory_path = sys.argv[1]
        print "The folder containing the .gck files is: " + target_directory_path
        target_directory_path_2 = sys.argv[2]
        print "The folder containing the .cluster files is: " + target_directory_path_2
        output_file_name = sys.argv[3]
        print "The output file name is: " + output_file_name
    else:
        print "Argument 1: You must provide the directory path which contains the gck files"
        print "Argument 2: You must provide the directory path which contains the cluster files"
        print "Argument 3: You must provide a name for the output file"
    
        sys.exit()
    #gck files
    file_names_list = os.listdir(target_directory_path)  #
    file_name_whole_list = []
    for file_name in file_names_list:
        if file_name.endswith('.gck'):
            file_name_whole = target_directory_path + "/" + file_name
            try:
                open(file_name_whole, "r")
            except IOError:
                print "Error: Couldnt open the file named: " + file_name_whole
                sys.exit()
            file_name_whole_list.append(file_name_whole)
            print ".gck file: " + file_name_whole + " " + "found"
    ####cluster files
    file_names_list_2 = os.listdir(target_directory_path_2)  #
    file_name_whole_list_2 = []

    genome_key_file = open("clusterxgenome_key", 'a+')
    genome_key_file_2 = open("genomexgenome_key", 'a+')

    for numerodefile,file_name in enumerate(file_names_list_2):
        if file_name.endswith('.cluster'):
            file_name_whole_2 = target_directory_path_2 + "/" + file_name
            try:
                open(file_name_whole, "r")
            except IOError:
                print "Error: Couldnt open the file named: " + file_name_whole_2
                sys.exit()
            file_name_whole_list_2.append(file_name_whole_2)
            print ".cluster file: " + file_name_whole_2 + " " + "found"
	    genome_key_file.write("%s\t%s\n" % (str(numerodefile),str(file_name_whole_2)))

    genome_key=get_genome_key(file_name_whole_list)

    findings  = []
    flag = 0
    flag_quiebre = 0
    best_of_the_best_ever_dict={}
    best_of_the_best_ever_entry=[]
    all_stats={}
    all_stats_entry=[]
    stats_increment=1
    all_stats[0]=["reported_cluster","quorum","minclust","indels","mean","std","n","hit_50%","hit_70%","hit_100%"]

    cluster_local_id_list=[]
    cluster_ranges_list=[]

    summary_output_name = output_file_name + "_summary"
    deleteContent(summary_output_name)
    output_file_pointer_2 = open(summary_output_name, 'a+')
    output_file_pointer_2.write("quorum\tminclust\tindels\tref_clust_id\tcluster_size_ref\tfound_clust_id\tcluster_size_found\tjaccard\thit_50%\thit_70%\thit_100%\treference_clust\tcomputed_clust\tother_genomes\tother_clusters_included\tcluster_members\n")
    loggylog = output_file_name + "_log"

    for ref_clust_id,file_path in enumerate(file_name_whole_list_2):

        file_name_short = file_path.split("/")
        file_name_cluster_crap=deepcopy(file_path)
        file_name_short = file_name_short[len(file_name_short) - 1]
        gn_file = file_name_short.split(".")[1]
        gn_file_pre_number = gn_file.strip("G")

        gn_file_number=genome_key[int(gn_file_pre_number)]
        #print gn_file_number

        with open(file_path, 'r') as file_chido:
            cluster_local_id_list = []
            for num_entries,line in enumerate(file_chido):

                            cluster_data = line.split('\t')
                            if cluster_data != "\n" and len(cluster_data) > 8:
                                #print cluster_data
                                cluster_local_id = cluster_data[7]
                                #print cluster_local_id
                                if cluster_local_id != "0":
                                    cluster_local_id_list.append(cluster_local_id)

            range_a_ref = int(cluster_local_id_list[0])
            range_b_ref = int(cluster_local_id_list[len(cluster_local_id_list)-1])

        stuff_in_between = []
        cluster_stuff  = []

        for file_path in file_name_whole_list:
            file_name_short = file_path.split("/")
            file_name_short = file_name_short[len(file_name_short) - 1]

            quorum = file_name_short.split("-")[1]
            quorum = int(quorum.strip("q"))
            minclus = file_name_short.split("-")[2]
            minclus = int(minclus.strip("mc"))
            indels = file_name_short.split("-")[3]
            indels = indels.strip("d")
            indels = int(indels.strip(".gck"))


            print str(quorum) + "/" + str(minclus) + "/" + str(indels)

            output_file_pointer_log = open(loggylog, 'a+')
            output_file_pointer_log.write("quorum: %s\tminclust: %s\tindels:  %s\t%s\t%s\t%s\tel id interno de gecko del genoma actual es: %s\n" % (str(quorum),str(minclus),str(indels),file_name_short,file_name_cluster_crap, file_name_short,str(gn_file_number)))

            flag1 = False
            flag2 = False
            flag3 = False
            flaga= False
            flagb= False
            with open(file_path,'r') as file_chido:
                for line in file_chido:
                    if flaga == True:
                        cluster_stuff.append(line)
                    if line == "<cluster>\n":
                        flaga= True
                    if line == "<occ>\n":
                        flag1 = True
                        flaga = False
                    if line == "</occ>\n":
                        flag3 = True
                        flag1 = False
                    if flag1 == True:
                        stuff_in_between.append(line)
                    if flag3 == True:
                        found_clust_id = cluster_stuff[0].split('\t')[0]

                        clusti_list_o_list=[]
                        for j, item in enumerate(stuff_in_between):
                            if j > 1:
                                item_plus=item.split('\t')
                                if item_plus[0] != "</occ>\n":
                                    #print item_plus
                                    clusti_list_o_list.append(item_plus)

                                    #print clusti_list_o_list

                        other_genomes_w_findings=[]
                        other_genomes=set()
                        for item in clusti_list_o_list:
                            #print item
                            # de aqu[i hay que ir a trav[es de esta lista a ver cual item concuerda con el clus referencia, de ser as[i se guarda
                            ##t odo, incluyendo todos los n[umero de lso genomas

                            if int(item[0]) == int(gn_file_number):
                                range_a_found = int(item[3])
                                range_b_found = int(item[4])

                                computed_clust = make_sequence(range_a_found,range_b_found,1)
                                reference_clust = make_sequence(range_a_ref, range_b_ref,1)
                                clust_len_found = len(computed_clust)
                                clust_len_ref = len(reference_clust)
                                jacindx=jaccardindex(computed_clust,reference_clust)
                                if jacindx > 0.0:
                                    hit50 = found_clust_func(reference_clust,computed_clust, 50.0)
                                    hit70 = found_clust_func(reference_clust,computed_clust, 70.0)
                                    hit100 = found_clust_func(reference_clust,computed_clust, 100.0)
                                    # print "el jaccard"
                                    # print jacindx
                                    # print "computed_clust"
                                    # print computed_clust
                                    # print reference_clust
                                    #      ("quorum\tminclust\tindels\tref_clust_id\tcluster_size_ref\tfound_clust_id\tcluster_size_found\tjaccard\thit_50%\thit_70%\thit_100%\tfamilies_ref\tfamilies_found\n" )

                                    for thing in clusti_list_o_list:
                                        clust_entry = thing[0] + "," + thing[2] + "," + thing[3] + "," + thing[4]
                                        other_genomes_w_findings.append(clust_entry)
                                        other_genomes.add(thing[0])

                                    otha_genm_list = list(other_genomes)
                                    output_file_pointer_2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(quorum),str(minclus),str(indels),str(ref_clust_id), str(clust_len_ref), str(found_clust_id), str(clust_len_found),str(jacindx), str(hit50), str(hit70), str(hit100), str(reference_clust), str(computed_clust),str(otha_genm_list),str(other_genomes_w_findings),str(cluster_stuff[1].strip("\n"))))

                        cluster_stuff = []                   #
                        stuff_in_between = []
                        flag1 = False
                        flag2 = False

                    flag3=False

    # summary_output_name1 = output_file_name + "50_summary"
    # output_file_pointer_2_1 = open(summary_output_name1, 'a+')
    # output_file_pointer_2_1.write(
    #     "quorum\tminclust\tindels\tref_clust_id\tcluster_size_ref\tfound_clust_id\tcluster_size_found\tjaccard\thit_50%\thit_70%\thit_100%\treference_clust\tcomputed_clust\tother_genomes\tother_clusters_included\n")
    #
    # summary_output_name2 = output_file_name + "70_summary"
    # output_file_pointer_2_2 = open(summary_output_name2, 'a+')
    # output_file_pointer_2_2.write(
    #     "quorum\tminclust\tindels\tref_clust_id\tcluster_size_ref\tfound_clust_id\tcluster_size_found\tjaccard\thit_50%\thit_70%\thit_100%\treference_clust\tcomputed_clust\tother_genomes\tother_clusters_included\n")
    #
    # summary_output_name3 = output_file_name + "100_summary"
    # output_file_pointer_2_3 = open(summary_output_name3, 'a+')
    # output_file_pointer_2_3.write(
    #     "quorum\tminclust\tindels\tref_clust_id\tcluster_size_ref\tfound_clust_id\tcluster_size_found\tjaccard\thit_50%\thit_70%\thit_100%\treference_clust\tcomputed_clust\tother_genomes\tother_clusters_included\n")
    #
    # with open(summary_output_name, 'r') as file_chido:
    #     for line_number,entry in enumerate(file_chido):
    #         entradas=entry.split("\t")
    #         #print entradas
    #         if line_number != 0:
    #             if int(entradas[8]) == 1:
    #                 output_file_pointer_2_1.write("%s\n" % str(entry))
    #             if int(entradas[9]) == 1:
    #                 output_file_pointer_2_2.write("%s\n" % str(entry))
    #             if int(entradas[10]) == 1:
    #                 output_file_pointer_2_3.write("%s\n" % str(entry))
    #
    #






