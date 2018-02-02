import os
import sys
import ast
from Bio import Entrez
import time


Entrez.email = "castillo@cebitec.uni-bielefeld.de"



def retrieve_annotation(id_list):

    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene",id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print "Not found"
        annotations = "NA"
        return annotations


    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key =
            queryKey)
    annotations = Entrez.read(data)

    #print "Retrieved %d annotations for %d genes" % (len(annotations),
    #        len(id_list))

    return annotations

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if __name__ == '__main__':

#### check arguments
    if len(sys.argv) > 3:
        target_directory_path = sys.argv[1]
        print "The folder containing the fastab files is: " + target_directory_path
    else:
        print "You must provide the folder where the .fastab files are"
        sys.exit()


    concat_other_clusts=""
    for i in range(2,len(sys.argv)):
        concat_other_clusts = concat_other_clusts + sys.argv[i]
        concat_other_clusts = concat_other_clusts + "\t"

    gene_list = ast.literal_eval(concat_other_clusts)


    print "genes to look for:"
    print gene_list

    print gene_list[0]
    print gene_list[1]




    for x1,sub_gene in enumerate(gene_list):
        print "cluster number " + str(x1)
        sub_sub_gene=sub_gene.split(",")
        g_file = "G" + sub_sub_gene[0] + ".fastab"
        g_file_path = target_directory_path + "/" + g_file
        lista_local_id =range(int(sub_sub_gene[2]), int(sub_sub_gene[3]) + 1, 1)



        for item in lista_local_id:
            with open(g_file_path) as file_pointer_a:  # open file
                for l_num,line1 in enumerate(file_pointer_a):  # read file by line
                    #print "linea numero:" + str(l_num)
                    #print line1

                    if (l_num + 1) == int(item):
                        linedata = line1.split('\t')  # tokenize by tab
                        name=linedata[1]
                        #print linedata
                        time.sleep(0.1)

                        handle = Entrez.esearch(db="gene", term=name)
                        record = Entrez.read(handle)
                        #print record
                        migeb=retrieve_annotation(record["IdList"])
                       # print migeb
                        if migeb != "NA":
                            name_1 = migeb["DocumentSummarySet"]["DocumentSummary"][0]["Name"]
                            symbol_1 = migeb["DocumentSummarySet"]["DocumentSummary"][0]["NomenclatureSymbol"]
                            #start_1 = migeb["DocumentSummarySet"]["DocumentSummary"][0]["GenomicInfo"][0]["ChrStart"]
                            #stop_1 = migeb["DocumentSummarySet"]["DocumentSummary"][0]["GenomicInfo"][0]["ChrStop"]
                            description_1 = migeb["DocumentSummarySet"]["DocumentSummary"][0]["Description"]
                            #print name_1 + "\t" +  symbol_1 + "\t" + start_1 + "\t" + stop_1 + "\t" + description_1
                            print name_1 + "\t" +  symbol_1 + "\t" + description_1


                        
                        #print migeb["Chromosome"]
                        #print migeb["NomenclatureSymbol"]

                        #print migeb["ChrStart"]
                        #print migeb["ChrStop"]
                        #print migeb["Description"]