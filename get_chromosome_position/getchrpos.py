from xml.dom.minidom import parse, parseString
import sys
import wget
from Bio import Entrez
### check arguments
####################################
input_file_path_batch = "genome_list.garnacha"
output_file_path = ""
species_name = ""
words = []
input_genomes=[]
species_name_search = ""
if len(sys.argv) == 1:
    try:
        open(input_file_path_batch, "r")
    except IOError:
        print('Error: Couldnt find the file named: genome_list.garnacha')
        sys.exit()
    print "Default path. Input file path: " + input_file_path_batch
    input_file_pointer = open(input_file_path_batch, "r")
    for line in input_file_pointer:
        line = line[:-1] + ".tabdom"
        input_genomes.append(line.replace(' ', '_'))
    for i in range(0, len(input_genomes)):
        print ("input genome %i: %s" % (i + 1, input_genomes[i]))

elif len(sys.argv) > 1:
    for i in range(1, len(sys.argv)):
        words.append(sys.argv[i])
    for i in range(0, len(words)):
        species_name = species_name + " " + words[i]
        species_name_search = species_name_search + "+" + words[i]
        if i < 1:
            output_file_path = output_file_path + " " + words[i]
        elif i == 1:
            output_file_path = output_file_path + " " + words[i] + ".txt"
    print ("Search for the species: %s" % (species_name))

print "Output file path: " + output_file_path
for i in range(0, len(input_genomes)):
    with open(input_genomes[i]) as file_pointer_input_file:  # open file
        for j,line in enumerate(file_pointer_input_file):  # read file by line
            if j != 0:

                linedata = line.split('\t')  # tokenize by tab
                print linedata
                ###Entry	Entry name	Gene names	Cross-reference (Pfam)	Cross-reference (InterPro)	Cross-reference (GeneID)	Organism ID	Organism	Cross-reference (RefSeq)	Gene ontology IDs
                search_GI=linedata[5].replace(";","")
                ###https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=3533391

                print search_GI
                Entrez.email = "castillo@cebitec.uni-bielefeld.de"

                handle = Entrez.esummary(db="gene", id=search_GI)
                record = Entrez.read(handle)
                #print record
                #print record.values()
                #print record["DocumentSummarySet"]
                type(record["DocumentSummarySet"]["DocumentSummary"])

               # print record["DocumentSummarySet"]["DocumentSummary"]

                input("presione para continuar...")