import os
import sys

import numpy

def deleteContent(fName):
    with open(fName, "w"):
        pass

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
    length_percentage = float(p)
    if inter_over_length >= length_percentage:
        return 1
    else:
        return 0

if __name__ == '__main__':

    if len(sys.argv) == 5:
        target_directory_path = sys.argv[1]
        print "The folder containing the .gck files is: " + target_directory_path
        target_directory_path_2 = sys.argv[2]
        print "The folder containing the .cluster files is: " + target_directory_path_2
        output_file_name = sys.argv[3]
        print "The output file name is: " + output_file_name
        famtype = sys.argv[4]
        print "Family type: " + famtype
    else:
        print "Argument 1: You must provide the directory path which contains the gck files"
        print "Argument 2: You must provide the directory path which contains the cluster files"
        print "Argument 3: You must provide a name for the output file"
        print "Argument 4: <egg|portho> according to wanting eggnog families or proteinortho"
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
    #cluster files
    file_names_list_2 = os.listdir(target_directory_path_2)  #
    file_name_whole_list_2 = []
    for file_name in file_names_list_2:
        if file_name.endswith('.cluster'):
            file_name_whole_2 = target_directory_path_2 + "/" + file_name
            try:
                open(file_name_whole, "r")
            except IOError:
                print "Error: Couldnt open the file named: " + file_name_whole_2
                sys.exit()
            file_name_whole_list_2.append(file_name_whole_2)
            print ".cluster file: " + file_name_whole_2 + " " + "found"

    deleteContent(output_file_name)

    #Read the .cluster files
    reported_clust_egg = {}
    reported_clust_egg_len = {}
    reported_clust_portho = {}
    reported_clust_portho_len = {}
    egg_cluster=[]
    num_entries = 0
    for i,file_path in enumerate(file_name_whole_list_2):
        file_name_short = file_path.split("/")
        file_name_short = file_name_short[len(file_name_short) - 1]
        egg_cluster = []
        if famtype == "egg":
            print "family type is eggnog"
            with open(file_path, 'r') as file_chido:
                for num_entries,line in enumerate(file_chido):
                    cluster_data = line.split("\t")
                    eggnog_fam = cluster_data[0]
                    if eggnog_fam != "\n":
                        if eggnog_fam != "0":
                            egg_cluster.append(eggnog_fam)
        if famtype == "portho":
            print "family type is proteinortho"
            with open(file_path, 'r') as file_chido:
                for num_entries,line in enumerate(file_chido):
                    line_1=line.strip()
                    cluster_data = line_1.split("\t")
                    if len(cluster_data) > 1:
                        eggnog_fam = cluster_data[8]
                        if eggnog_fam != "\n":
                            if eggnog_fam != "0":
                                egg_cluster.append(eggnog_fam.strip(" "))
        if famtype != "portho" and famtype != "egg":
            print "Family type is not valid"
            sys.exit()

        reported_clust_egg[i]=egg_cluster # minus one because of the newline item
        reported_clust_egg_len[i] = num_entries
    #Read the .gck files
    findings  = []
    flag = 0
    flag_quiebre = 0
    best_of_the_best_ever_dict={}
    best_of_the_best_ever_entry=[]
    all_stats={}
    all_stats_entry=[]
    stats_increment=1
    all_stats[0]=["reported_cluster","quorum","minclust","indels","mean","std","n","hit_50%","hit_70%","hit_100%"]
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

        #print str(quorum) + "/" + str(minclus) + "/" + str(indels)



        findings = []
        dict_count = 1
        clust_dict = {}
        clust_dict_len = {}
        cluster_1 = str()
        cluster_2 = []
        with open(file_path,'r') as file_chido:
            #print findings
            flag = 0
            counter = 0
            for line in file_chido:
                if line == "<occ>\n":
                    flag = 0
                    counter=0
                    cluster_1 = str()
                    cluster_2 = []
                if flag == 1 and counter == 1:
                    cluster_1 = cluster_1 + line
                    #print cluster_1
                    fam_list=list(cluster_1.split(","))
                    cluster_length=len(fam_list)

                    for index1,item in enumerate(fam_list):
                        if index1 == 0:
                            cluster_2.append(fam_list[index1][1:].strip(" "))
                        elif index1 == len(fam_list)-1:
                            cluster_2.append(fam_list[len(fam_list)-1][:-2].strip(" "))
                        else:
                            cluster_2.append(fam_list[index1].strip(" "))
                    clust_dict[dict_count] = cluster_2
                    clust_dict_len[dict_count] = cluster_length
                    dict_count = dict_count + 1

                if flag == 1:
                    counter=1
                if line == "<cluster>\n":
                    flag = 1


        #print file_name_short
        #print clust_dict
        jacij=float()
        dict_best_clust={}
        current_best_clust=[]
        dict_best_clust[0]=["genome_a_id", "size_a", "genome_b_id", "size_b", "jaccard","hit_50%","hit_70%","hit_100%", "fams_a", "fams_b"]
        shared_clust_count=1

        p1=float(0.5)
        p1_lista=[]
        p2=float(0.7)
        p2_lista = []
        p3=float(1.0)
        p3_lista = []
        for keyi, valuei in reported_clust_egg.iteritems():
            #the best entry per reported cluster
            best_of_the_best_ever_entry = [0,0,0,0,float(0.0),0,0,0,"x","x"]
            jaccardlist = []
            all_stats_entry = []
            p1_lista = []
            p2_lista = []
            p3_lista = []
            for keyj, valuej in clust_dict.iteritems():
                 jacij=jaccardindex(valuei, valuej) #input are lists that are going to be converted into sets, output is an integer
                 found_or_not_5 = found_clust_func(valuei, valuej, p1)
                 found_or_not_7 = found_clust_func(valuei, valuej, p2)
                 found_or_not_10 = found_clust_func(valuei, valuej, p3)


                 if jacij > 0.0:

                    current_best_clust = [keyi,reported_clust_egg_len[keyi],keyj,len(valuej), jacij,found_or_not_5,found_or_not_7,found_or_not_10,valuei,valuej]
                    dict_best_clust[shared_clust_count] = current_best_clust
                    jaccardlist.append(jacij)
                    p1_lista.append(found_or_not_5)
                    p2_lista.append(found_or_not_7)
                    p3_lista.append(found_or_not_10)
                    shared_clust_count = shared_clust_count + 1


                    # keep the best of the best
                    if jacij > best_of_the_best_ever_entry[4]:
                        best_of_the_best_ever_entry = current_best_clust


            single_file_name_2 = str(file_name_whole_list_2[keyi])
            list_of_single_file_name_2 = single_file_name_2.split("/")
            single_file_name_2 = list_of_single_file_name_2[len(list_of_single_file_name_2) - 1]
            list_of_single_file_name_2 = single_file_name_2.split(".")
            single_file_name_2 = list_of_single_file_name_2[0]

            whole_output_name = output_file_name + "-" + famtype + "-" + single_file_name_2 + str(quorum) + "-" + str(minclus) + "-" + str(indels)
            deleteContent(whole_output_name)
            output_file_pointer = open(whole_output_name, 'a+')
            i = 0
            while i < len(dict_best_clust):
                output_file_pointer.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(dict_best_clust[i][0]),str(dict_best_clust[i][1]),str(dict_best_clust[i][2]),str(dict_best_clust[i][3]),str(dict_best_clust[i][4]),str(dict_best_clust[i][5]),str(dict_best_clust[i][6]),str(dict_best_clust[i][7]),str(dict_best_clust[i][8]),str(dict_best_clust[i][9])))
                i += 1


            #save stats
            if len(p1_lista) == 0:
                hit_50 = 0
            else:
                hit_50 = max(p1_lista)
            if len(p2_lista) == 0:
                hit_70 = 0
            else:
                hit_70 = max(p2_lista)
            if len(p3_lista) == 0:
                hit_100 = 0
            else:
                hit_100 = max(p3_lista)

            all_stats_entry = [keyi,quorum,minclus,indels,numpy.mean(jaccardlist),numpy.std(jaccardlist),len(jaccardlist), hit_50,hit_70,hit_100]

            #print all_stats_entry

            all_stats[stats_increment] = all_stats_entry
            stats_increment = stats_increment + 1
        # print "Finished"



            # save the best
            entry_name_best_of_best= single_file_name_2 + "-" + str(quorum) + "-" + str(minclus)  + "-" + str(indels)
            best_of_the_best_ever_dict[entry_name_best_of_best] = best_of_the_best_ever_entry

    summary_output_name = output_file_name + "_summary_" + famtype
    deleteContent(summary_output_name)
    output_file_pointer_2 = open(summary_output_name, 'a+')

    output_file_pointer_2.write("quorum\tminclust\tindels\tref_clust_id\tcluster_size_ref\tfound_clust_id\tcluster_size_found\tjaccard\thit_50%\thit_70%\thit_100%\tfamilies_ref\tfamilies_found\n" )

    for key, item in best_of_the_best_ever_dict.iteritems():
        furtherinfo=str(key)
        quorum=furtherinfo.split("-")[1]
        minclus=furtherinfo.split("-")[2]
        indels=furtherinfo.split("-")[3]


        output_file_pointer_2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (quorum,minclus,indels,str(item[0]), str(item[1]), str(item[2]), str(item[3]),str(item[4]), str(item[5]), str(item[6]), str(item[7]), str(item[8]), str(item[9])))

    stats_output_name = output_file_name + "_all_stats_stuff_" + famtype
    deleteContent(stats_output_name)
    output_file_pointer_3 = open(stats_output_name, 'a+')

    for key, item in all_stats.iteritems():
        #all_stats[0] = ["reported_cluster", "quorum", "minclust", "indels", "mean", "std", "n"]

        rep_clust = item[0]
        quorum = item[1]
        minclus = item[2]
        indels = item[3]
        mean = item[4]
        std = item[5]
        size_n = item[6]
        hit_50 = item[7]
        hit_70 = item[8]
        hit_100 = item[9]
        output_file_pointer_3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(rep_clust), str(quorum), str(minclus), str(indels), str(mean), str(std), str(size_n),str(hit_50),str(hit_70),str(hit_100)))
