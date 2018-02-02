import os
import sys
import re


##gets the paths of the individual files in the given directory
def get_dom_file_paths_from_dir(target_directory_path):
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



## make a dictionary with all the domains found in all the genomes and gives an integer to each individual one
## domain_dict[<Domain name>]= <assigned integer>
def get_all_the_domainids(dom_file_name_whole_list):
    domains=set()
    domain_dict = {}
    for file in dom_file_name_whole_list:
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
    i=0
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

####################################
############################################
############################################################
############################################################################
####################################################################################
if __name__ == '__main__':

#### check arguments
    print "###########################################################"
    print "##  Make cog file for Gecko input with protein domains   ##"
    print "###########################################################"
    if len(sys.argv) == 4:
        target_directory_path = sys.argv[1]
        print "The folder containing the files is: " + target_directory_path
    else:
        print "Arg1: You must provide the directory path which contains the .pfam.tsv files"
        print "Arg2: You must provide an output file name"
        print "Arg3: pfam or ipr, whether you want pfam or interpro annotations"
        sys.exit()

    output_file_name = sys.argv[2]
    cog_domain_file_name = target_directory_path + "/" + output_file_name + ".cog"

    print "The output file name is: " + cog_domain_file_name

    annota = sys.argv[3]
    if regex_test('pfam|ipr', annota):
        print "Annotation is: " + str(annota)
    else:
        print "Error, must provide a string which is \"pfam\" or \"ipr\""
        sys.exit()

    if annota == "pfam":
        annota_column = int(4)
    else:
        annota_column = int(11)

    ##get the file paths from the folder
    dom_file_name_whole_list = get_dom_file_paths_from_dir(target_directory_path)
    ##make the dictionary of domains
    domain_dict = get_all_the_domainids(dom_file_name_whole_list)
    domain_file_pointer = open("domaindictionary", 'w+')
    for key, value in domain_dict.iteritems() :
        domain_file_pointer.write("%s\t%s\n" % (str(key),str(value)))
    #print domain_dict



    file_name_fastab_whole_list = get_file_paths_from_dir_fastab(target_directory_path)

    if len(file_name_fastab_whole_list) != len(dom_file_name_whole_list):
        print "number of files must match"
        sys.exit()

    #Making a list of dictionaries containing the genomes
    dicts_of_dicts_of_genomes = {}


    sim_file_pointer = open(cog_domain_file_name, 'w')
    sim_file_pointer.close()
##### data sample
##gb5161|2116216|locus|2_#_5056_#_75529|9270545:9271573|chromosome|chromosome_2|strand|+	d6ad767136b1a32060f1716642967f36	301	Pfam	PF06747	CHCH domain	188	224	1.9E-6	T	17-04-2017	IPR010625	CHCH
    print "making domain dictionary"
    print dom_file_name_whole_list



    for i, file_a in enumerate(dom_file_name_whole_list):

        current_genome_number = file_a.split("/")[-1].strip(".pfam.tsv").strip("G")
        print "Current genome: " + current_genome_number
    
        dicts_of_dicts_of_genomes[int(current_genome_number)] = {}
        domains = set()
        # ridicoulous number to make finding the first entry easy
        previous_gene_id = int(99999999)
        ### open genome file
        with open(file_a) as file_pointer_a:  # open file
            for line1 in file_pointer_a:  # read file by line
                linedata = line1.split('\t')  # tokenize by tab
                #print linedata
                local_id_genome_a = linedata[0].split("|")[0].strip("gb")
                # print local_id_genome_a
                # print previous_gene_id
                if ((int(local_id_genome_a) == int(previous_gene_id)) or (int(previous_gene_id) == int(99999999))):
                    # print "entro"
                    if len(linedata) > annota_column:
                        domain_id_a = linedata[annota_column]
                        #print domain_id_a
                        #local_domain_id_a = domain_dict[domain_id_a]
                        domains.add(domain_id_a)
                        #print domains
                        previous_gene_id = local_id_genome_a
                    else:
                        previous_gene_id = local_id_genome_a
                else:
                    #print "found domains"
                    #print domains

                    dicts_of_dicts_of_genomes[int(current_genome_number)][int(previous_gene_id)] = list(domains)
                    # print dicts_of_dicts_of_genomes[int(current_genome_number)]

                    domains = set()
                    if len(linedata) > annota_column:
                        domain_id_a = linedata[annota_column]
                        #local_domain_id_a = domain_dict[domain_id_a]
                        domains.add(domain_id_a)
                        # print domains
                        previous_gene_id = local_id_genome_a
                    else:
                        previous_gene_id = local_id_genome_a
                        # raw_input()
                        # print dicts_of_dicts_of_genomes
    print "Dictionary is done"
    #print dicts_of_dicts_of_genomes[4]
    #print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    #print dicts_of_dicts_of_genomes[2]

##G10
#10401	proteins
#0	-	?	XP_001208347.1	unknown	ATEG_00982	unknown
#0	-	?	XP_001208348.1	unknown	ATEG_00983	unknown
#0	+	?	XP_001208349.1	unknown	ATEG_00984	unknown
#0	+	?	XP_001208350.1	unknown	ATEG_00985	unknown
#ENOG4119FWE	-	?	XP_001208351.1	unknown	ATEG_00986	unknown
#0	+	?	XP_001208352.1	unknown	ATEG_00987	unknown
#0	-	?	XP_001208353.1	unknown	ATEG_00988	unknown
#0	-	?	XP_001208354.1	unknown	ATEG_00989	unknown
#ENOG4119G3N	+	?	XP_001208355.1	unknown	ATEG_00990	unknown
#0	-	?	XP_001208356.1	unknown	ATEG_00991	unknown
#0	+	?	XP_001208357.1	unknown	ATEG_00992	unknown



    lista_dom_a = []
    lista_dom_b = []
    domains_in_a = set()
    domains_in_b = set()
    foundy_found = False
    for i in range(1,len(file_name_fastab_whole_list)+1):

        file_a = target_directory_path + "/G" + str(i) + ".fastab"
        file_pointer_a = open(file_a, 'r')
        proteinas = file_len(file_a)
        #COG_like file header

        sim_file_pointer = open(cog_domain_file_name, 'a')
        sim_file_pointer.write("G%s\n" % (str(i)))
        sim_file_pointer.write("%s\tproteins\n" % (str(proteinas)))

	individual_cog_file_path = target_directory_path + "/G" + str(i) +  ".geckable_dom"
  	individual_cog_pointer = open(individual_cog_file_path, 'w')
    	individual_cog_pointer.close()
	individual_cog_pointer = open(individual_cog_file_path, 'a')
        individual_cog_pointer.write("G%s\n" % (str(i)))
        individual_cog_pointer.write("%s\tproteins\n" % (str(proteinas)))

        for iii, line1 in enumerate(file_pointer_a):  # read file by line
            linedata_a = line1.split('\t')  # tokenize by tab
            protein_id = linedata_a[0]
            gene_id = linedata_a[1]
            strand = linedata_a[6].strip("\n")
            local_id_genome_a = iii + 1

            #check if the local id can be found in the dictionary
            tamali_a = dicts_of_dicts_of_genomes[int(i)].get(int(local_id_genome_a), "notfound")
            if tamali_a == "notfound":
                fam = "0"
                sim_file_pointer.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    str(fam), str(strand), str("?"), str(protein_id), str("unknown"), str(gene_id),
                    str(local_id_genome_a)))
		individual_cog_pointer.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    str(fam), str(strand), str("?"), str(protein_id), str("unknown"), str(gene_id),
                    str(local_id_genome_a)))	
            else:
                fam = tamali_a
                for dominio in fam:
                    sim_file_pointer.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( str(dominio), str(strand),str("?"),str(protein_id), str("unknown"),str(gene_id), str(local_id_genome_a) ))
		    individual_cog_pointer.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( str(dominio), str(strand),str("?"),str(protein_id), str("unknown"),str(gene_id), str(local_id_genome_a) ))
                    foundy_found = False

        sim_file_pointer.write("\n")
