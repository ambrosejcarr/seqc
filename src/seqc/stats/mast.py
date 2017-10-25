import math
import subprocess
import imp
import os
import pandas as pd
import numpy as np


def run_mast(counts_filtered, clustering_communities, output_prefix):
    # Differentially Expression Analysis using MAST
    log_counts = (counts_filtered + 1.0).applymap(math.log2)
    de_results = []  # array containing the differentially expression analysis for each cluster
    for c in range(np.max(clustering_communities) + 1):
        tmp_input_file = output_prefix + "_cluster_" + str(c) + "_mast_input.csv"
        tmp_output_file = output_prefix + "_cluster_" + str(c) + "_mast_results.csv"
        reduced_tdf1 = log_counts.iloc[np.where(clustering_communities == c)[0]]
        reduced_tdf2 = log_counts.iloc[np.where(clustering_communities != c)[0]]
        reduced_df = pd.concat([reduced_tdf1, reduced_tdf2])
        reduced_df.index = pd.Index([1 if i < len(reduced_tdf1.index) else 0 for i in range(len(reduced_tdf1.index) + len(reduced_tdf2.index))])
        reduced_df.to_csv(tmp_input_file)

        path_to_run_mast = imp.find_module('seqc')[1]
        args = 'Rscript {p} {i} {o}'.format(p=os.path.join(path_to_run_mast, 'run_mast.R'), i=tmp_input_file, o=tmp_output_file)
        with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True) as p:
            out, err = p.communicate()
            if os.path.isfile(tmp_output_file):
                de_gene_df = pd.read_csv(tmp_output_file)
                if len(de_gene_df.index) > 0:
                    de_results.append(de_gene_df)
                else:  # if no differentially expressed genes
                    de_results.append(None)
            else:
                de_results.append(None)

    de_gene_list_file = output_prefix + "_de_gene_list.txt"
    with open(de_gene_list_file, "w") as f:
        f.write("Differential Expression Analysis Using MAST\n\n")
        c = 1
        for de_result in de_results:
            if de_result is not None:
                f.write("Differentially expressed genes for cluster %d:\n" % (c))
                f.write("%-10s  %-10s  %-10s  %-10s\n" % ("Gene", "p", "p.fdr", "logFC"))

                for i in range(len(de_result)):
                    p_v = "%.2e" % de_result.loc[i][1]
                    p_fdr = "%.2e" % de_result.loc[i][2]
                    logFC = "%.2f" % de_result.loc[i][3]
                    f.write("%-10s  %-10s  %-10s  %-10s\n" % (de_result.loc[i][0], p_v, p_fdr, logFC))
            else:
                f.write("No differentially expressed genes has been found for cluster %d.\n" % (c))
            c += 1
            f.write("\n")
        f.close()
    return de_gene_list_file