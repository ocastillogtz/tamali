

import os
import sys
import re
import numpy
import networkx as nx
import itertools
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

def get_file_paths_from_ffgc_logs_error(target_directory_path):
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

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)

def regex_test(rule, target):
    regexp = re.compile(rule)
    if regexp.search(target):
        return True
    else:
        return False

def get_file_paths_from_dir_fastab(target_directory_path):
    file_names_list = os.listdir(target_directory_path)
    files_pfam = [i for i in file_names_list if i.endswith('.fastab')]
    file_name_whole_list = []
    for file_name in files_pfam:
        file_name_whole = target_directory_path + "/" + file_name
        try:
            open(file_name_whole, "r")
        except IOError:
            print "Error: Couldnt open the file named: " + file_name_whole
            sys.exit()
        file_name_whole_list.append(file_name_whole)
        print "Fastab file: " + file_name_whole + " " + "found"
    return file_name_whole_list

def get_all_the_domainids_from_cog(file):
    domains=set()
    domain_dict = {}

    with open(file) as file_pointer:  # open file
        for line in file_pointer:  # read file by line
            linedata = line.split('\t')  # tokenize by tab

            if len(linedata) > 4:
                if not regex_test('^G\d+\s|^G\d+\n', line):
                    if not regex_test('^\d+\s+proteins', line):
                        if not regex_test('^\n|^\s+', line):
                            domain_id = linedata[0]
                            domains.add(domain_id)
    for i,item in enumerate(domains):
        domain_dict[item] = i
    return domain_dict

## make a dictionary with all the domains found in all the genomes and gives an integer to each individual one
## domain_dict[<Domain name>]= <assigned integer>




if __name__ == '__main__':

    #### check arguments
    print "########################################################"
    print "##  Making graph file our of the cog file for gecko   ##"
    print "########################################################"

    print "Arg1: file path for the cog file for gecko"
    print "Arg2: \"domains\" or \"raw\""
    print "Arg3: actual path"
    if len(sys.argv) == 4:
        target_directory_path = sys.argv[1]
        print "The folder containing the files is: " + target_directory_path
    else:
        print "You must provide the file path for the cog file for gecko"
        sys.exit()
    if not sys.argv[2] == "raw" or sys.argv[2] == "domains":
        print "Arg2: \"domains\" or \"raw\""
        sys.exit()
    
    current_dir=sys.argv[3]
    previous_linea=99999999
    my_graph = nx.Graph()
    cuenta=0
    starting_cuenta=[]
    overall_data_filename="graphdataall"
    deleteContent(overall_data_filename)
    oFile_1 = open(overall_data_filename, "w")
    file_name_list=[]
    if sys.argv[2] == "domains":
        domain_dict = get_all_the_domainids_from_cog(target_directory_path)
        domain_file_pointer = open("domaindictionary_cog_to_gt", 'w+')
        for key, value in domain_dict.iteritems():
            domain_file_pointer.write("%s\t%s\n" % (str(key), str(value)))
    with open(target_directory_path, 'r') as cog_file_pointer:
        for num_entries, line in enumerate(cog_file_pointer):
            if regex_test('^G\d+\s|^G\d+\n', line):
                starting_cuenta = []
                genome_name=re.findall('^G\d+', line)
                #print "El nombre del genoma es: " + str(genome_name[0])
                output_file_name=str(genome_name[0]) + ".gml"
                file_name_list.append(output_file_name)
                deleteContent(output_file_name)
                oFile=open(output_file_name,"w")
                oFile.write("graph [\n    directed 0\n")


            if not regex_test('^\n|^\s+', line.split("\t")[0]):
                if len(line.split("\t")) > 3:
                    linea_separa = line.split("\t")

                    #IPR014031	+	?	BAE66272.1	unknown	AO090010000404	11713
                    if linea_separa[0] != 0:
                        starting_cuenta.append(cuenta)
                        family = linea_separa[0]
                        if sys.argv[2] == "domains":
                            family_number= domain_dict[linea_separa[0]]
                        else:
                            family_number=family
                        protein_id=linea_separa[3]
                        gene_id = linea_separa[5]
                        local_id = linea_separa[6].strip("\n")
                        #,protein=protein_id,gene=gene_id,local_id=local_id,genome=genome_name
                        #my_graph.add_node(num_entries, Class=linea_separa[0],label=int(local_id),id=int(cuenta),'protein_id'=str(protein_id))
                        if family == 0 or family == "0":
                            family_number = str(cuenta) + "NA"
                        inside_data= str(cuenta) + "\t" + str(family) + "\t" + str(family_number) + "\t" + str(protein_id) + "\t" + str(gene_id) + "\t" + str(local_id) + "\n"
                        stringo="    node [\n        id " + str(cuenta) + "\n        label " + str(cuenta) +  "\n        class \"" + str(family_number) + "\"\n    ]\n"
                        oFile.write("%s\n" % stringo)
                        oFile_1.write("%s" % inside_data)



                        #my_graph.node[cuenta]['protein_id'] = str(protein_id)
                        #my_graph.node[cuenta]['gene_id'] = str(gene_id)
                        #my_graph.node[cuenta]['local_id'] = str(local_id)
                        if previous_linea != 99999999:
                            my_graph.add_edge(num_entries-1, num_entries)
                    previous_linea=line
            else:
                previous_linea = 99999999
                #print int(starting_cuenta[0])
                #print cuenta
                for i in range(int(starting_cuenta[0]),cuenta - 1):
                    orig=i
                    dest=i+1
                    stringo = "edge [\n        source " + str(orig) + "\n        target " + str(dest) + "\n        weight 1\n    ]\n"
                    oFile.write("%s" % stringo)
                oFile.write("]")
                #graph_file_name = str(genome_name[0]) + ".gml"
                #nx.write_gml(my_graph, graph_file_name)
                starting_cuenta=[]
            cuenta=cuenta+1
    print file_name_list

    quorum_paths_filename="filepathlistsquorum"
    deleteContent(quorum_paths_filename)
    quorumpaths_pointer = open(quorum_paths_filename, "w")
    file_name_list_plus=[]
    for item in file_name_list:

        full_path= str(current_dir) + "/" + str(item)
        file_name_list_plus.append(full_path)



    for index_2 in range(1,len(file_name_list_plus)):
        iteracion = list(itertools.combinations(file_name_list_plus, index_2))
        for items_comb in iteracion:
            quorum = len(items_comb)
            quorumpaths_pointer.write("%s\t" % str(quorum))
            for items_inside_comb in items_comb:
                quorumpaths_pointer.write("%s " % str(items_inside_comb))
            quorumpaths_pointer.write("\n")




