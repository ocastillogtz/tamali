import os
import sys

#For IPR = 11, for PF equals 4
annota_column = int(11)


##gets the paths of the individual files in the given directory
def get_file_paths_from_dir(target_directory_path):
    file_names_list = os.listdir(target_directory_path)
    files_pfam = [i for i in file_names_list if i.endswith('.pfam.tsv')]
    file_name_whole_list = []
    for file_name in files_pfam:
        file_name_whole = target_directory_path + "/" + file_name
        try:
            open(file_name_whole, "r")
        except IOError:
            print "Error: Couldnt open the file named: " + file_name_whole
            sys.exit()
        file_name_whole_list.append(file_name_whole)
        print "pfam file: " + file_name_whole + " " + "found"
    return file_name_whole_list

def deleteContent(fName):
    with open(fName, "w"):
        pass

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




## make a dictionary with all the domains found in all the genomes and gives an integer to each individual one
## domain_dict[<Domain name>]= <assigned integer>
def get_all_the_domainids(file_name_whole_list):
    domains=set()
    domain_dict = {}
    for file in file_name_whole_list:
        with open(file) as file_pointer:  # open file
            for line in file_pointer:  # read file by line
                linedata = line.split('\t')  # tokenize by tab

                if len(linedata) > annota_column:
                    domain_id = linedata[annota_column]
                    domains.add(domain_id)
    for i,item in enumerate(domains):
        domain_dict[item] = i
    return domain_dict

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if __name__ == '__main__':

#### check arguments
    if len(sys.argv) == 2:
        target_directory_path = sys.argv[1]
        print "The folder containing the files is: " + target_directory_path
    else:
        print "You must provide the directory path which contains the .pfam.tsv files"
        sys.exit()
    ##get the file paths from the folder
    file_name_whole_list = get_file_paths_from_dir(target_directory_path)
    ##make the dictionary of domains
    domain_dict = get_all_the_domainids(file_name_whole_list)
    domain_file_pointer = open("domaindictionary", 'w+')
    for key, value in domain_dict.iteritems() :
        domain_file_pointer.write("%s\t%s\n" % (str(key),str(value)))
    #print domain_dict

    file_name_fastab_whole_list = get_file_paths_from_dir_fastab(target_directory_path)

    #Making a list of dictionaries containing the genomes
    dicts_of_dicts_of_genomes = {}
##### data sample
##gb5161|2116216|locus|2_#_5056_#_75529|9270545:9271573|chromosome|chromosome_2|strand|+	d6ad767136b1a32060f1716642967f36	301	Pfam	PF06747	CHCH domain	188	224	1.9E-6	T	17-04-2017	IPR010625	CHCH
    print "Domain dictionary done"




    for i, file_a in enumerate(file_name_whole_list):
        current_genome_number = file_a.split("/")[-1].strip(".pfam.tsv").strip("G")
        print "Current genome: " + current_genome_number

        dicts_of_dicts_of_genomes[int(current_genome_number)]={}
        domains=set()
        #ridicoulous number to make finding the first entry easy
        previous_gene_id=int(99999999)
        ### open genome file
        with open(file_a) as file_pointer_a:  # open file
            for line1 in file_pointer_a:  # read file by line
                linedata = line1.split('\t')  # tokenize by tab
                #print linedata
                local_id_genome_a = linedata[0].split("|")[0].strip("gb")
                #print local_id_genome_a
                #print previous_gene_id
                if ((int(local_id_genome_a) == int(previous_gene_id)) or (int(previous_gene_id) == int(99999999))):
                    #print "entro"
                    if len(linedata) > annota_column:
                        domain_id_a = linedata[annota_column]
                        local_domain_id_a = domain_dict[domain_id_a]
                        domains.add(local_domain_id_a)
                        #print domains
                        previous_gene_id = local_id_genome_a
                    else:
                        previous_gene_id = local_id_genome_a
                else:
                    #print "found domains"
                    #print domains
                    dicts_of_dicts_of_genomes[int(current_genome_number)][int(previous_gene_id)] = list(domains)
                    #print dicts_of_dicts_of_genomes[int(current_genome_number)]
                    domains = set()
                    if len(linedata) > annota_column:
                        domain_id_a = linedata[annota_column]
                        local_domain_id_a = domain_dict[domain_id_a]
                        domains.add(local_domain_id_a)
                        #print domains
                        previous_gene_id = local_id_genome_a
                    else:
                        previous_gene_id = local_id_genome_a
                    #raw_input()
                #print dicts_of_dicts_of_genomes
    print "Dictionary is done"
    #print dicts_of_dicts_of_genomes[3]
    #print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    #print dicts_of_dicts_of_genomes[2]
    # raw_input()


    lista_dom_a=[]
    lista_dom_b=[]
    domains_in_a = set()
    domains_in_b = set()
    foundy_found=False
    for i,file_a in enumerate(file_name_fastab_whole_list):
        file_pointer_a = open(file_a, 'r')
        short_fila_a = file_a.split("/")[-1].strip(".fastab")
        shorter_file_a = short_fila_a.strip("G")
        for iii, line1 in enumerate(file_pointer_a):  # read file by line
            linedata_a = line1.split('\t')  # tokenize by tab
            local_id_genome_a = iii + 1
            tamali_a = dicts_of_dicts_of_genomes[int(shorter_file_a)].get(int(local_id_genome_a), "notfound")
            if tamali_a == "notfound":
                foundy_found = False
            else:
                foundy_found = True

            if foundy_found == True:

                for j,file_b in enumerate(file_name_fastab_whole_list):
                    #print file_b
                    file_pointer_b = open(file_b, 'r')
                    #print str(i) + " " + str(j) + " " + file_a + " " + file_b
                    if i != j:
                        short_fila_b = file_b.split("/")[-1].strip(".fastab")
                        shorter_file_b = short_fila_b.strip("G")
                        output_file_name = str(short_fila_a) + "_" + str(short_fila_b) + ".sim"
                        sim_file_pointer = open(output_file_name, 'a')


                        #print output_file_name

                        for jjj,line2 in enumerate(file_pointer_b):  # read file by line
                            linedata_b = line2.split('\t')  # tokenize by tab
                            local_id_genome_b = jjj + 1
                            tamali_b = dicts_of_dicts_of_genomes[int(shorter_file_b)].get(int(local_id_genome_b),"notfound")

                            #domains_in_a = set(dicts_of_dicts_of_genomes[int(shorter_file_a)][local_id_genome_a])
                            if tamali_b != "notfound":


                                lista_dom_a=list(dicts_of_dicts_of_genomes[int(shorter_file_a)][int(local_id_genome_a)])
                                lista_dom_b=list(dicts_of_dicts_of_genomes[int(shorter_file_b)][int(local_id_genome_b)])
                                domains_in_a=set(lista_dom_a)
                                domains_in_b=set(lista_dom_b)

                                ##print "gene in A:" + str(domains_in_a)
                                ##print "gene in B:" + str(domains_in_b)
                                intersec_domains=domains_in_a.intersection(domains_in_b)

                                if len(intersec_domains) > 0:
                                    # print "el genoma A: " + str(short_fila_a)
                                    # print "El gen con entradas en el diccionario encontrado en el genoma A es:"
                                    # print "la llave es: " + str(local_id_genome_a)
                                    # print dicts_of_dicts_of_genomes[int(shorter_file_a)][int(local_id_genome_a)]
                                    # print "el genoma B:" + str(shorter_file_b)
                                    # print "El gen con entradas en el diccionario encontrado en el genoma B es:"
                                    # print "la llave es: " + str(local_id_genome_b)
                                    # print  dicts_of_dicts_of_genomes[int(shorter_file_b)][int(local_id_genome_b)]
                                    # print "la cardinalidad de la intersecci[on: " + str(len(intersec_domains))
                                    
                                    sim_file_pointer.write("%s\t%s\t%s\t%s\t%d\t%s\n" % (str(short_fila_a), str(local_id_genome_a),str(short_fila_b),str(local_id_genome_b),int(1),str("0.5")))
                foundy_found=False
