import pandas as pd
import sys
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np


def plotty_by_hitpercentage(hit_column_name,output_filename):

    newthing = df.groupby(["quorum", "indels","ref_clust_id"])[hit_column_name].sum()
    summary1_filename = "sumario_de_" + hit_column_name.strip("%")
    newthing.to_csv(summary1_filename, sep='\t')
    df_50_1 = pd.read_csv(summary1_filename, sep='\t',  names = ["quorum", "indels", "ref_clust_id", hit_column_name], index_col=False)
    mask = df_50_1[hit_column_name] > 1
    column_name = hit_column_name
    df_50_1.loc[mask, column_name] = 1
    newthing_1 = df_50_1.groupby(["quorum", "indels"])[hit_column_name].sum()
    summary2_filename = "sumario_de_" + hit_column_name.strip("%") + "_a"
    newthing_1.to_csv(summary2_filename, sep='\t')
    df_50_2 = pd.read_csv(summary2_filename, sep='\t', names=["quorum", "indels", hit_column_name],index_col=False)
    max_quorum = df_50_1['quorum'].max()
    min_quorum = df_50_1['quorum'].min()
    max_indels = df_50_1['indels'].max()
    min_indels = df_50_1['indels'].min()
    quorums_list = []

    for value in range(min_quorum, max_quorum+1):
        if any(df_50_2.quorum == value):
            quorums_list.append(value)



    for i in quorums_list:
        filtered = df_50_2[(df_50_2['quorum'] == i)]
        for indice_1 in range(min_indels, max_indels):
            if not any(filtered.indels == indice_1):
                new_row = [i, indice_1, 0]
                df_50_2.loc[len(df_50_2)] = new_row

    fig, axs = plt.subplots(len(quorums_list), 1 , figsize=(9, 14), facecolor='w')
    fig.subplots_adjust(hspace=.5, wspace=.001)

    axs = axs.ravel()



    for i,j in enumerate(quorums_list[::-1]):

        filtered = df_50_2[(df_50_2['quorum'] == j)]

        axs[i].set_axisbelow(True)
        axs[i].yaxis.grid(color="#b2b2b2", linestyle='-', linewidth=1)
        axs[i].bar(filtered['indels'].tolist(), filtered[hit_column_name].tolist(), align='center',color="#165b42", linewidth = 0)
        axs[i].set_title("Quorum " + str(j), rotation='vertical',x=-0.05,y=0.7, size = "smaller")
        axs[i].set_xlim([min_indels-1, max_indels+1])
        axs[i].set_ylim([0, 30])

        rects = axs[i].patches

        # Now make some labels
        #labels = filtered['hit_50%'].tolist()
        #labels = ["%d Clts" % i for i in filtered['hit_50%'].tolist()]
        labels=[]
        for item in filtered[hit_column_name].tolist():
            labels.append(str(item) + " clusters")
        #print labels
        for rect, label in zip(rects, labels):
            height = rect.get_height()
            axs[i].text(rect.get_x() + rect.get_width() / 2, height + 2, label, ha='center', va='bottom', size = "smaller")

    #plt.show()
    fig_file_name = hit_column_name.strip("%") + output_filename +".png"
    plt.savefig(fig_file_name, dpi = 300)


if __name__ == '__main__':

    if len(sys.argv) == 3:
        datafile_path = sys.argv[1]
        print "Data file for plotting: " + datafile_path
        output_filename = sys.argv[2]


    else:
        print "Argument 1: You must provide the the path of your datafile"
        print "Argument 2: Output file name"
        sys.exit()

    df = pd.read_csv(datafile_path, sep='\t' ,usecols=range(0,11), index_col=False)

    print "plotting 50% hits"
    plotty_by_hitpercentage("hit_50%",output_filename)
    print "plotting 70% hits"
    plotty_by_hitpercentage("hit_70%",output_filename)
    print "plotting 100% hits"
    plotty_by_hitpercentage("hit_100%",output_filename)

