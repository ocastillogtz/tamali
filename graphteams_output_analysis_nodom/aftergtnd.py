import os
import sys
import re
from copy import deepcopy

def deleteContent(fName):
    with open(fName, "w"):
        pass
def make_sequence(a, b, increment):
    list = []
    for i in range(a, b + 1, increment):
        list.append(i)
    return list

def regex_test(rule, target):
    regexp = re.compile(rule)
    if regexp.search(target):
        return True
    else:
        return False

# def get_genome_key(file_name_whole):
#     genome_key_dict = {}
#     with open(file_name_whole[0], 'r') as file_chido:
#         contador = 0
#         flag1 = False
#         for line in file_chido:
#             if flag1 == True:
#                 genome_key_dict[int(line.strip("G"))] = contador
#                 genome_key_file_2.write("%s is the genome number %s\n" % (str(line), str(contador)))
#
#                 contador = contador + 1
#                 flag1 = False
#             if line == "<genome>\n":
#                 flag1 = True
#             if line == "</genomes>\n":
#                 break
#     # print genome_key_dict
#     return genome_key_dict

def jaccardindex(list_a, list_b):
    set_a = set(list_a)
    set_b = set(list_b)
    union_chida = set_a.union(set_b)
    interseccion_chida = set_a.intersection(set_b)
    jaccard_index = float(len(interseccion_chida)) / float(len(union_chida))
    return jaccard_index

def found_clust_func(list_a, list_b, p):
    set_a = set(list_a)
    set_b = set(list_b)
    interseccion_chida = set_a.intersection(set_b)
    set_a_length = float(len(set_a))
    intersec_len = float(len(interseccion_chida))
    inter_over_length = intersec_len / set_a_length
    length_percentage = float(p / 100.0)

    if inter_over_length >= length_percentage:
        return 1
    else:
        return 0


if __name__ == '__main__':

    if len(sys.argv) == 4:

        target_directory_path = sys.argv[1]
        print "The folder containing the graph_team files are: " + target_directory_path
        target_directory_path_2 = sys.argv[2]
        print "The folder containing the cluster files are: " + target_directory_path
        output_file_name = sys.argv[3]
        print "The output file name is: " + output_file_name

    else:
        print "Argument 1: You must provide the directory path which contains the graph_team files are"
        print "Argument 2: You must provide the directory path which contains the .cluster files are"
        print "Argument 3: Output name "

        sys.exit()
######################################################################### ALL MODE




    file_names_list = os.listdir(target_directory_path)  #
    file_name_whole_list = []
    print "archivos encontrados en el folder de la output de graph teams"
    print str(len(file_names_list))
    for file_name in file_names_list:
        if regex_test('graph_teams',file_name):

            file_name_whole = target_directory_path + "/" + file_name
            try:
                open(file_name_whole, "r")
            except IOError:
                print "Error: Couldnt open the file named: " + file_name_whole
                sys.exit()
            file_name_whole_list.append(file_name_whole)
            #print "graph_teams file: " + file_name_whole + " " + "found"


    ####cluster files
    file_names_list_2 = os.listdir(target_directory_path_2)  #
    file_name_whole_list_2 = []

    genome_key_file = open("clusterxgenome_key", 'a+')
    #genome_key_file_2 = open("genomexgenome_key", 'a+')

    summary_output_name = output_file_name + "_summary"
    deleteContent(summary_output_name)
    output_file_pointer_2 = open(summary_output_name, 'a+')
    output_file_pointer_2.write(
        "quorum\tminclust\tindels\tref_clust_id\tcluster_size_ref\tfound_clust_id\tcluster_size_found\tjaccard\thit_50%\thit_70%\thit_100%\treference_clust\tcomputed_clust\tother_genomes\tother_clusters_included\tcluster_members\n")

    for numerodefile, file_name in enumerate(file_names_list_2):
        if file_name.endswith('.cluster'):
            file_name_whole_2 = target_directory_path_2 + "/" + file_name
            try:
                open(file_name_whole_2, "r")
            except IOError:
                print "Error: Couldnt open the file named: " + file_name_whole_2
                sys.exit()
            file_name_whole_list_2.append(file_name_whole_2)
            print ".cluster file: " + file_name_whole_2 + " " + "found"
            genome_key_file.write("%s\t%s\n" % (str(numerodefile), str(file_name_whole_2)))


    ###make big summary file with all the clusters found by graph_teams



        # graph_teams - d1 - q7 - uniq10111415158
        # graph_teams - d5 - q7 - uniq111324579

    print str(len(file_name_whole_list))
#####################################################################################################################
    for ref_clust_id, file_path in enumerate(file_name_whole_list_2):
        file_name_short = file_path.split("/")
        file_name_cluster_crap = deepcopy(file_path)
        file_name_short = file_name_short[len(file_name_short) - 1]
        gn_file = file_name_short.split(".")[1]
        gn_file_pre_number = gn_file.strip("G")
        #print "The actual genome file is: " + gn_file

        with open(file_path, 'r') as file_chido:
            #print file_path
            cluster_local_id_list = []
            for num_entries, line in enumerate(file_chido):
                cluster_data = line.split('\t')
                if cluster_data != "\n" and len(cluster_data) > 8:
                    # print cluster_data
                    cluster_local_id = cluster_data[7]
                    # print cluster_local_id
                    if cluster_local_id != "0":
                        cluster_local_id_list.append(cluster_local_id)
         #   print "el contenido del cluster referencia es "
          #  print cluster_local_id_list
            range_a_ref = int(cluster_local_id_list[0])
            range_b_ref = int(cluster_local_id_list[len(cluster_local_id_list) - 1])
            clust_len_ref = len(cluster_local_id_list)
            reference_clust = cluster_local_id_list
            ######looping en los archivos de los resultados de graph teams
            genomes_in_this = []
            cluster_members = []
            cluster_entries = 0
            computed_clust = []
            other_genomes=[]


            for file_path_1 in file_name_whole_list:

                # print "my graphfile is: "
                # print file_path_1
                file_name_short = file_path_1.split("/")
                file_name_short_1 = file_name_short[len(file_name_short) - 1]
                indels = file_name_short_1.split("-")[1].strip("d")
                quorum = file_name_short_1.split("-")[2].strip("q")
                minclus = 3

                with open(file_path_1, 'r') as graph_teams_file_pointer:
                    counter_1 = 0
                    for line_1 in graph_teams_file_pointer:
                       # print line_1
                        # regex_text(rule, text)
                        if regex_test("^No\s", line_1):
                          #  print "no clusters found"
                            break
                        else:

                            # get first line details, the genomes involved in the quorum of this thing
                            if regex_test("^common", line_1):
                                genomes_in_this = []
                                genomes_involved_paths = line_1.split("\t")[1:]
                                for genome_path in genomes_involved_paths:
                                    genomes_in_this.append(genome_path.split("/")[-1].strip("\n").strip(".gml"))
                               # print genomes_in_this


                            # check further lines for clusters with members 3 or more
                            else:
                                cluster_members=[]
                                if len(genomes_in_this) > 1:
                                    if gn_file in genomes_in_this:
                                        for j, genome_file_name in enumerate(genomes_in_this):
                                            if genome_file_name == gn_file:
                                                ref_cluster_genome_field = j
                                            other_genomes.append(genome_file_name.strip("G"))
                                            #          print "the column where the content of the genome " + str(gn_file) + "  is the field/column " + str(ref_cluster_genome_field)
                                        cluster_data_fields = line_1.split("\t")
                                        for item in cluster_data_fields[0].split(";"):
                                            cluster_members.append(item)
                                            #         print cluster_members
                                        if len(cluster_members) > 2:

                                         #   print "el indice de la columna es "
                                          #  print ref_cluster_genome_field
                                            computed_clust=[]
                                            cluster_genome_entry = cluster_data_fields[ref_cluster_genome_field].split(
                                                ";")
                                            for things_1 in cluster_genome_entry:
                                                computed_clust.append(things_1)
                                            clust_len_found = len(computed_clust)
                                            found_clust_id = counter_1
                                            otha_genm_list = genomes_in_this
                                            cluster_stuff = cluster_members
                                            other_genomes_w_findings = "NA"
                                            clust_len_found = len(computed_clust)
                                            clust_len_ref = len(reference_clust)
                                            jacindx = jaccardindex(computed_clust, reference_clust)

                                            if jacindx > 0.2:
                                                hit50 = found_clust_func(reference_clust, computed_clust, 50.0)
                                                hit70 = found_clust_func(reference_clust, computed_clust, 70.0)
                                                hit100 = found_clust_func(reference_clust, computed_clust, 100.0)

                                                output_file_pointer_2.write(
                                                    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                                                        str(quorum), str(minclus), str(indels), str(ref_clust_id),
                                                        str(clust_len_ref),
                                                        str(found_clust_id), str(clust_len_found), str(jacindx), str(hit50),
                                                        str(hit70),
                                                        str(hit100), str(reference_clust), str(computed_clust),
                                                        str(otha_genm_list),
                                                        str(other_genomes_w_findings), str(cluster_stuff)))

                                else:
                                    print "Error you can't start reading the cluster data without having the info on which genomes are involved"
                                    raw_input()
                                    break

                            counter_1 = counter_1 + 1








