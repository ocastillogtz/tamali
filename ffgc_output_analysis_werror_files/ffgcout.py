


#/vol/cluster-data/castillo/A_ffgcv03test_updt_gen/bin/proyecto_danis_fungi_reduced_updated/weak_common_intervals/
#macsi_evalue=0.00001_NS_n2_l10_s0.5_d1_m4_q2

#macsi.mwci
#/vol/cluster-data/castillo/A_ffgcv03test_updt_gen/bin/proyecto_danis_fungi_reduced_updated/weak_common_intervals/macsi_evalue=0.00001_NS_n2_l10_s0.5_d1_m4_q2/macsi.mwci

import os
import sys
import numpy
import datetime



def deleteContent(fName):
    with open(fName, "w"):
        pass

def make_sequence(a, b, increment):
    list = []
    for i in range(int(a), int(b)+1, increment):
        list.append(i)
    return list

def get_genome_key(file_name_whole):
    genome_key_dict = {}
    with open(file_name_whole[0], 'r') as file_chido:
        contador = 0
        flag1 = False
        for line in file_chido:
            if flag1 == True:
                genome_key_dict[int(line.strip("G"))] = contador
                contador = contador + 1
                flag1 = False
            if line == "<genome>\n":
                flag1 = True
            if line == "</genomes>\n":
                break
    # print genome_key_dict
    return genome_key_dict

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

def get_gene_real_position(target_dir_3,G_number,protein_dom_pos):
    file_full_path = target_dir_3 + "/" + "G" + str(G_number) + ".geckable_dom"
    with open(file_full_path, 'r') as file_chido:
        for num_entries, line in enumerate(file_chido):
            if num_entries == (protein_dom_pos + 2):
                #print line

                real_local_id = line.split("\t")[6]
                return real_local_id
        return "error"


def file_len(fname):
    with open(fname) as f:
        global cuenta_1
        cuenta_1 = 0
        for cuenta_1, l in enumerate(f):
            pass
    return cuenta_1 + 1


def get_file_paths_from_dir(target_directory_path,extension):
    file_names_list = os.listdir(target_directory_path)
    files_ext = [i for i in file_names_list if i.endswith(extension)]
    file_name_whole_list = []
    for file_name in files_ext:
        file_name_whole = target_directory_path + "/" + file_name
        try:
            open(file_name_whole, "r")
        except IOError:
            print "Error: Couldnt open the file named: " + file_name_whole
            sys.exit()


        file_name_whole_list.append(file_name_whole)
        print extension + "File: " + file_name_whole + " " + "found"
    return file_name_whole_list


#get_file_paths_from_ffgc_wci_folders

def get_file_paths_from_ffgc_wci_folders(target_directory_path):
    file_names_list = os.listdir(target_directory_path)
    file_name_whole_list = []
    for file_name in file_names_list:
        record = True
        file_name_whole = target_directory_path + "/" + file_name + "/macsi.mwci"
        try:
            open(file_name_whole, "r")
        except IOError:
            print "Error: Couldnt open the file named: " + file_name_whole
            errand_log_pointer = open("errorlog", 'a+')
            #t = datetime.time(1, 2, 3)
            errand_log_pointer.write("error: Couldn't open the file: %s\n" % file_name_whole)
            record = False
        if record:
            file_name_whole_list.append(file_name_whole)
            print "File: " + file_name_whole + " " + "found"
    return file_name_whole_list

def indel_calc(members_A,members_B):
    return len(members_A.symmetric_difference(members_B))


class cluster:
    def __init__(self, genome, start, end, members,genome_ref,start_ref,end_ref,members_ref):
        self.genome = genome
        self.start = start
        self.end = end
        self.members = members
        self.genome_ref = genome_ref
        self.start_ref = start_ref
        self.end_ref = end_ref
        self.members_ref = members_ref




if __name__ == '__main__':

    print "Start after FFGC using blast cutoff relationships analysis"
    if len(sys.argv) == 4:
        target_directory_path = sys.argv[1]
        print "The path for the weak_common_intervals folder: " + target_directory_path
        target_directory_path_2 = sys.argv[2]
        print "The folder containing the .cluster files is: " + target_directory_path_2
        output_file_name = sys.argv[3]
        print "The output file name is: " + output_file_name
    else:
        print "Argument 1: You must provide the path for the weak_common_intervals folder"
        print "Argument 2: You must provide the directory path which contains the cluster files"
        print "Argument 3: You must provide a name for the output file"
        sys.exit()
    # gck files

    # /vol/cluster-data/castillo/A_ffgcv03test_updt_gen/bin/proyecto_danis_fungi_reduced_updated/weak_common_intervals/
    # macsi_evalue=0.00001_NS_n2_l10_s0.5_d1_m4_q2
    # macsi.mwci
    # /vol/cluster-data/castillo/A_ffgcv03test_updt_gen/bin/proyecto_danis_fungi_reduced_updated/weak_common_intervals/macsi_evalue=0.00001_NS_n2_l10_s0.5_d1_m4_q2/macsi.mwci

    #### mwci files
    file_name_whole_list = get_file_paths_from_ffgc_wci_folders(target_directory_path)
    ####cluster files
    file_name_whole_list_2 = get_file_paths_from_dir(target_directory_path_2, ".cluster")

    findings = []
    flag = 0
    flag_quiebre = 0
    best_of_the_best_ever_dict = {}
    best_of_the_best_ever_entry = []
    all_stats = {}
    all_stats_entry = []
    stats_increment = 1
    all_stats[0] = ["reported_cluster", "quorum", "minclust", "indels", "mean", "std", "n", "hit_50%", "hit_70%",
                    "hit_100%"]

    cluster_local_id_list = []
    cluster_ranges_list = []

    summary_output_name = output_file_name + "_summary"
    deleteContent(summary_output_name)
    output_file_pointer_2 = open(summary_output_name, 'a+')
    output_file_pointer_2.write(
        "quorum\tminclust\tindels\tref_clust_id\tcluster_size_ref\tfound_clust_id\tcluster_size_found\tjaccard\thit_50%\thit_70%\thit_100%\treference_clust\tcomputed_clust\tother_genomes\tother_clusters_included\n")

    gene_location_set = set()

    list_of_clusters=[]
    flag_1 = False
    flag_2 = False
    for ref_clust_id, file_path in enumerate(file_name_whole_list_2):

        file_name_short = file_path.split("/")
        file_name_short = file_name_short[len(file_name_short) - 1]
        file_name_short_clust = file_name_short
        gn_file = file_name_short.split(".")[1]
        print "Looking for the cluster:" + file_name_short
        print gn_file
        #gn_file_pre_number = gn_file.strip("G")

        #gn_file_number = genome_key[int(gn_file_pre_number)]
        # print gn_file_number

        with open(file_path, 'r') as file_chido:
            cluster_local_id_list = []
            for num_entries, line in enumerate(file_chido):

                cluster_data = line.split('\t')
                if cluster_data != "\n" and len(cluster_data) > 8:
                    # print cluster_data
                    cluster_local_id = cluster_data[7]
                    # print cluster_local_id
                    if cluster_local_id != "0":
                        cluster_local_id_list.append(cluster_local_id)

            range_a_ref = int(cluster_local_id_list[0])
            range_b_ref = int(cluster_local_id_list[len(cluster_local_id_list) - 1])



        stuff_in_between = []
        cluster_stuff = []
        #gn_file is the full G# name
        #now loop on the mwci files
        # /vol/cluster-data/castillo/A_ffgcv03test_updt_gen/bin/proyecto_danis_fungi_reduced_updated/weak_common_intervals/macsi_evalue=0.00001_NS_n2_l10_s0.5_d1_m4_q2/macsi.mwci

        genomes_other_found = []
        range_a_other_found = []
        range_b_other_found = []
        genoms_w_clust_was_found = []
        for file_path_1 in file_name_whole_list:
            file_name_short = file_path_1.split("/")
            file_name_short = file_name_short[len(file_name_short) - 2]
            file_name_shorter = file_name_short.split("_")
            quorum = file_name_shorter[- 1]
            quorum = int(quorum.strip("q"))
            minclus = file_name_shorter[- 2]
            minclus = int(minclus.strip("m"))
            indels = file_name_shorter[- 3]
            indels = indels.strip("d")

            #print file_name_short
            #print "Quorum: " + str(quorum) + "/ " + "minclustsize: " + str(minclus) + "/ " + "indels: " + str(indels)
            actual_GA_range_a=int(0)
            actual_GA_range_b=int(0)
            past_GA_range_a=int(999999999)
            past_GA_range_b=int(999999999)
            contador_1=int(0)
            genes_inclust_ref=[]

            max_lines_ffgc_out=file_len(file_path_1)
            file_chido = open(file_path_1, 'r')



            for line_num,line in enumerate(file_chido):
                line_seped = line.split("\t")
                genome_A = line_seped[0]
                genome_B = line_seped[4]
                GA_range_a = line_seped[2]
                GA_range_b = line_seped[3]
                GB_range_a = line_seped[6]
                GB_range_b = line_seped[7]
# G1	BA000049	1	3	G10	NT_165939	164	165
# G1	BA000049	1	3	G14	NC_007199	7965	7966
# G1	BA000049	2	5	G14	NC_007199	7965	7966


                if genome_A == gn_file and int(GA_range_a) >= int(range_a_ref) and int(GA_range_a) <= int(range_b_ref) or genome_A == gn_file and int(GA_range_b) >= int(range_a_ref) and int(GA_range_b) <= int(range_b_ref):

                    actual_GA_range_a = int(GA_range_a)
                    actual_GA_range_b = int(GA_range_b)

                    if actual_GA_range_a == past_GA_range_a and actual_GA_range_b == past_GA_range_b or past_GA_range_a == int(999999999) and past_GA_range_b == int(999999999):
                        past_GA_range_a = actual_GA_range_a
                        past_GA_range_b = actual_GA_range_b


                        if int(GB_range_a) < int(GB_range_b):
                            genes_inclust = make_sequence(int(GB_range_a), int(GB_range_b), 1)
                            if int(GA_range_a) < int(GA_range_b):
                                genes_inclust_ref = make_sequence(GA_range_a, GA_range_b, 1)
                                # list_of_clusters.append(cluster(genome_B,GB_range_a,GB_range_b, set(genes_inclust),genome_A,GA_range_a,GA_range_b,genes_inclust_ref))
                                list_of_clusters.append(
                                    [genome_B, GB_range_a, GB_range_b.strip("\n"), set(genes_inclust), genome_A,
                                     GA_range_a, GA_range_b, genes_inclust_ref])
                        else:
                            genes_inclust = make_sequence(GB_range_b, GB_range_a, 1)
                            if int(GA_range_a) > int(GA_range_b):
                                genes_inclust_ref = make_sequence(GA_range_b, GA_range_a, 1)
                                # list_of_clusters.append(cluster(genome_B, GB_range_a, GB_range_b, set(genes_inclust), genome_A, GA_range_b,GA_range_a, genes_inclust_ref))
                                list_of_clusters.append(
                                    [genome_B, GB_range_a, GB_range_b.strip("\n"), set(genes_inclust), genome_A,
                                     GA_range_b, GA_range_a, genes_inclust_ref])

                            else:
                                genes_inclust_ref = make_sequence(GA_range_b, GA_range_a, 1)
                                list_of_clusters.append(
                                    [genome_B, GB_range_a, GB_range_b.strip("\n"), set(genes_inclust), genome_A,
                                     GA_range_a, GA_range_b, genes_inclust_ref])
                    else:

                        other_genomes_found=set()
                        other_clusters=[]
                        #indels_inside
                        for clust_item in list_of_clusters:
                            #print clust_item
                            other_genomes_found.add(clust_item[0])
                            indels_inside=indel_calc(clust_item[3], clust_item[7])
                            the_other_cluster=[clust_item[0],indels_inside,clust_item[1],clust_item[2]]
                            other_clusters.append(the_other_cluster)

                        clust_len_ref=len(cluster_local_id_list)
                        found_clust_id=contador_1
                        clust_len_found=len(genes_inclust_ref)

                        jacindx=jaccardindex(cluster_local_id_list, map(str,genes_inclust_ref))
                        #print jacindx
                        hit25 = found_clust_func(cluster_local_id_list, map(str,genes_inclust_ref), 25.0)
                        hit50 = found_clust_func(cluster_local_id_list, map(str,genes_inclust_ref), 50.0)
                        hit70 = found_clust_func(cluster_local_id_list, map(str,genes_inclust_ref), 70.0)
                        hit100 = found_clust_func(cluster_local_id_list, map(str,genes_inclust_ref), 100.0)

                        reference_clust=cluster_local_id_list
                        computed_clust=list(genes_inclust_ref)
                        otha_genm_list=list(other_genomes_found)
                        other_genomes_w_findings=list(other_clusters)
                        if hit25 == 1:
                            output_file_pointer_2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(quorum), str(minclus), str(indels), str(ref_clust_id), str(clust_len_ref),str(found_clust_id), str(clust_len_found), str(jacindx), str(hit50), str(hit70),str(hit100), str(reference_clust), str(computed_clust), str(otha_genm_list),str(other_genomes_w_findings)))
                        list_of_clusters = []
                        contador_1 = contador_1 + 1
                        flag_2 == False
                        flag_1 == False
                        past_GA_range_a = actual_GA_range_a
                        past_GA_range_b = actual_GA_range_b
                        if int(GB_range_a) < int(GB_range_b):
                            genes_inclust = make_sequence(int(GB_range_a), int(GB_range_b), 1)
                            if int(GA_range_a) < int(GA_range_b):
                                genes_inclust_ref = make_sequence(GA_range_a, GA_range_b, 1)
                                #list_of_clusters.append(cluster(genome_B,GB_range_a,GB_range_b, set(genes_inclust),genome_A,GA_range_a,GA_range_b,genes_inclust_ref))
                                list_of_clusters.append([genome_B,GB_range_a,GB_range_b.strip("\n"), set(genes_inclust),genome_A,GA_range_a,GA_range_b,genes_inclust_ref])
                        else:
                            genes_inclust = make_sequence(GB_range_b, GB_range_a, 1)
                            if int(GA_range_a) > int(GA_range_b):
                                genes_inclust_ref = make_sequence(GA_range_b, GA_range_a, 1)
                                #list_of_clusters.append(cluster(genome_B, GB_range_a, GB_range_b, set(genes_inclust), genome_A, GA_range_b,GA_range_a, genes_inclust_ref))
                                list_of_clusters.append([genome_B, GB_range_a, GB_range_b.strip("\n"), set(genes_inclust), genome_A, GA_range_b,GA_range_a, genes_inclust_ref])

                            else:
                                genes_inclust_ref = make_sequence(GA_range_b, GA_range_a, 1)
                                list_of_clusters.append([genome_B, GB_range_a, GB_range_b.strip("\n"), set(genes_inclust), genome_A, GA_range_a,GA_range_b, genes_inclust_ref])


