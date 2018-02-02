import os
import sys

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

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

if __name__ == '__main__':

#### check arguments
    if len(sys.argv) == 3:
        target_directory_path = sys.argv[1]
        print "The folder containing the cluster files is: " + target_directory_path
        target_directory_path_1 = sys.argv[2]
        print "The folder containing the geckable_porto files is: " + target_directory_path_1
    else:
        print "You must provide the directory path which contains the .cluster files (1st arg)"
        print "You must provide the directory path which contains the .fastab files (2st arg)"
        sys.exit()

    cluster_file_list = get_file_paths_from_dir(target_directory_path, ".cluster")
    fastab_file_list = get_file_paths_from_dir(target_directory_path_1, ".geckable_portho")

    output_folder_path = "updatedclusts"
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
    # column in which you find the gene annotation in cluster file
    column_with_gene_annot = 5
    # column in which you find the gene annotation in geckable_portho file
    column_with_gene_annot_1 = 5
    contador=0
    for i, file_path in enumerate(cluster_file_list):
        if i > 0:
            if contador == 0:
                print "Couldn't generate cluster file " + str(file_name_short) + " it had no matches with their genenames"
        file_name_short = file_path.split("/")
        file_name_short = file_name_short[len(file_name_short) - 1]
        gn_file = file_name_short.split(".")[1]
        #print gn_file

        contador=0
        with open(file_path, 'r') as file_ptr:
            for num_entries, line in enumerate(file_ptr):
                gene_name = line.split("\t")[column_with_gene_annot]

                for j, file_path_1 in enumerate(fastab_file_list):
                    file_name_short_1 = file_path_1.split("/")
                    file_name_short_1 = file_name_short_1[len(file_name_short_1) - 1]
                    gn_file_1 = file_name_short_1.split(".")[0]

                    with open(file_path_1, 'r') as file_ptr_1:
                        for num_entries_1, line_1 in enumerate(file_ptr_1):
                            if len(line_1.split("\t")) >= 5:

                                gene_name_1 = line_1.split("\t")[column_with_gene_annot_1]

                                if gene_name == gene_name_1 or gene_name.split(".")[0] == gene_name_1.split(".")[0]:
                                    contador=contador+1
                                    summary_output_name = "./" + output_folder_path + "/" + file_name_short.split(".")[0] + "." + gn_file_1 + ".cluster"
                                    print summary_output_name
                                    output_file_pointer = open(summary_output_name, 'a+')
                                    #print line_1
                                    #print "in cluster: " + gene_name
                                    #print gn_file
                                    #print "in geckable_portho: " + gene_name_1
                                    #print gn_file_1
                                    # plus 1 because it starts in line zero instead of 1 and minus 2 because the cog format header that uses two lines
                                    local_id = num_entries_1 - 1
                                    #print "in geckable_portho line number: " + str(local_id)
                                    item_1 = line_1.split("\t")[0]
                                    item_2 = line_1.split("\t")[1]
                                    item_3 = line_1.split("\t")[2]
                                    item_4 = line_1.split("\t")[3]
                                    item_5 = line_1.split("\t")[4]
                                    item_6 = line_1.split("\t")[5]


                                    output_file_pointer.write("%s\t%s\t%s\t%s\t%s\t%s\tunknown\t%s\t0\n" % (item_1,item_2,item_3,item_4,item_5,item_6,local_id))




