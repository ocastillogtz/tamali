from collections import Counter
import os
import sys
import numpy
import datetime
import re
def deleteContent(fName):
    with open(fName, "w"):
        pass


def file_len(fname):
    with open(fname) as f:
        global cuenta_1
        cuenta_1 = 0
        for cuenta_1, l in enumerate(f):
            pass
    return cuenta_1 + 1

def regex_test(rule, target):
    regexp = re.compile(rule)
    if regexp.search(target):
        return True
    else:
        return False

if __name__ == '__main__':

        #### check arguments
    print "##############################################"
    print "##  Removing promiscuous genes in cog file  ##"
    print "##############################################"

    print "Arg1: cog file path"
    print "Arg2: integer value the percentil of the protein domain families frequencies, recommended 95"
    print "Arg3: the suffix of the output file suffix"
   
    if len(sys.argv) == 4:
        target_directory_path = sys.argv[1]
        print "The cog file is here: " + target_directory_path

    else:
        print "Argument 1: You must provide the path for the cog file"
        sys.exit()

    try:
        open(target_directory_path, "r")
    except IOError:
        print "Error: Couldnt open the cog file in: " + target_directory_path
        sys.exit()

    percentil_value  = sys.argv[2]
    print "the Protein domain family frequency percentil cut out is: " + percentil_value
    outputsufix = sys.argv[3]
    promis_fams_report_filename= "promis_fams_report"

#G10
#10401	proteins
#0	-	?	XP_001218754.1	unknown	ATEG_10406	unknown
#0	+	?	XP_001218753.1	unknown	ATEG_10405	unknown
#0	-	?	XP_001218575.1	unknown	ATEG_10227	unknown
#15155	-	?	XP_001218576.1	unknown	ATEG_10228	unknown
#0	-	?	XP_001218577.1	unknown	ATEG_10229	unknown
#0	-	?	XP_001218578.1	unknown	ATEG_10230	unknown
    summary_output_name = target_directory_path.split("/")[-1].strip(".cog") + str(outputsufix) +".cog"
    print "Output file is: " + summary_output_name
    deleteContent(summary_output_name)
    output_file_pointer = open(summary_output_name, 'a+')


    families_list=[]
    with open(target_directory_path, 'r') as file_chido:
        for num_entries, line in enumerate(file_chido):
            if regex_test("^G\d+",line) == False and regex_test("^\d+\s+proteins",line) == False and str(line) != "" and regex_test("^\s+$",line) == False:
                #print line
                if str(line.split("\t")[0]) != "0":
                    families_list.append(line.split("\t")[0])


    la_cuenta = Counter(families_list)
#find the most promiscuos fams
   
    los_mas_entrometidos=[]

    lista_frecuencia=[]
    for item_1 in la_cuenta:
        lista_frecuencia.append(la_cuenta[item_1])
    cut_off = numpy.percentile(lista_frecuencia, int(percentil_value))



    left_out_because_their_promiscuity=[]
    for item_2 in la_cuenta:
        if int(la_cuenta[item_2]) > int(cut_off):
            left_out_because_their_promiscuity.append(item_2)


    promis = set(left_out_because_their_promiscuity)
    deleteContent(promis_fams_report_filename)
    promis_fam_rep_pointer=open(promis_fams_report_filename, 'a+')
    promis_fam_rep_pointer.write("Promiscous gene families\n")
    promis_fam_rep_pointer.write("Cutting out the %s percentile\n" % str(percentil_value))
    for member in promis:
        promis_fam_rep_pointer.write("%s\n" % str(member))

    with open(target_directory_path, 'r') as file_chido:
        for num_entries, line in enumerate(file_chido):
            if str(line.split("\t")[0]) in promis:
                entries_max_len=len(line.split("\t"))
                for j in range(0,entries_max_len):
                    if j == 0:
                        output_file_pointer.write("0\t")
                    else:
                        output_file_pointer.write("%s\t" % str(line.split("\t")[j].strip('\n')))
                output_file_pointer.write("\n")
            else:
                output_file_pointer.write("%s\n" % str(line.strip('\n')))




