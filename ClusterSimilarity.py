# !/usr/bin/python3
"""This script is to be used in conjunction with the Camoco software.
It will compare a set of genes within a cluster across networks and
prioritize the biggest connestion differences for a set of genes
"""

# even though only a few packages are called in this script
# all of them are needed for Camoco software built in functions
import sys
import os
import pandas as pd
import numpy as np
import scipy as sp
import camoco as co
import collections
import itertools
from itertools import chain
import getopt


def usage():
    print('''
          This is the usage function:
          python3 ClusterEnrichment.py -c <COB> -s <COB2>

          -c or --cob The name of your co-expression network
          -s or --secondnetwork The name of your comparision co-expression network''')

try:
    opts, args = getopt.getopt(sys.argv[1:],
                               "c:s:h",
                               ["cob=", "secondnetwork=", "help"])

except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-c", "--cob"):
        cobname = arg
    elif opt in ("-s", "--secondnetwork"):
        cob2name = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"


# load the network object
cob = co.COB(cobname)
cob_compare = co.COB(cob2name)

cob.set_sig_edge_zscore(2.5)
cob_compare.set_sig_edge_zscore(2.5)

Ortho = pd.read_csv("/home/mich0391/Downloads/Updated_Brachy_Wheat_Ortho.csv")
#Ortho = pd.read_csv("/home/mich0391/Downloads/CamocoWork/Vahid/Brachy_Wheat_orthologs.csv")

# you can input a list
# Ortho.loc[Ortho["W2691"].isin(blort)]
# Get the index of the pandas column that matches gene
# test = Ortho[Ortho.W2691 == "TRAESCS7D01G316200"].index[-1]
# Use that index to get the Ortholog
# Ortho.get_value(Ortho.index[test[-1]], "Bd21_NGFF")
# change from np.ndarray to pd dataframe
Clusters = pd.DataFrame(cob.clusters)

# Make a ordered dictionary to each key
# is a cluster and each value is the
# genes in that cluster
ClustDict = collections.OrderedDict()
for index, row in Clusters.iterrows():
    if row[0] in ClustDict.keys():
        ClustDict[row[0]].append(index)
    else:
        ClustDict[row[0]] = [index]

# Output header information
print("ClusterID", "\t",
      cobname, "_Density", "\t",
      cob2name, "_Gene_Density", "\t",
      "Density_Difference", "\t",
      "Percent_Density_Difference",
      )

cob_densities = list()
cob_compare_densities = list()
Differences = list()
Percentages = list()
# Iterate through the cluster list numbers and measure their density
# go through each cluster and set of genes
for Cluster, Genes in ClustDict.items():
    ClusterID = Cluster
    # convert gene list to IDs in network
    Gene_List = (cob.refgen.from_ids(Genes))
    if len(Genes) >= 9:
        cob_density = cob.density(Gene_List, min_distance=2.5)
        # Acess ortholog dataframe and get rows with matching gene names
        OrthologsIndex = Ortho.loc[Ortho[cobname].isin(Genes)]
        # Index the second network column and convert it to a list
        Orthologs = OrthologsIndex[cob2name].tolist()
        Ortho_Gene_List = (cob_compare.refgen.from_ids(Orthologs))
        cob_compare_density = cob_compare.density(Ortho_Gene_List)
        cob_compare_density = cob_compare.density(Ortho_Gene_List,
                                                  min_distance=2.5)
        cob_compare_densities.append(cob_compare_density)
        Difference = (cob_density - cob_compare_density)
        Differences.append(Difference)
        Percent_Difference = 100 * (Difference/cob_density)
        Percentages.append(Percent_Difference)
        print(Cluster,
              cob_density,
              cob_compare_density,
              Difference,
              Percent_Difference,
              sep="\t")

print("Averages: ",
      "\n",
      sum(cob_densities)/len(cob_densities),
      sum(cob_compare_densities)/len(cob_compare_densities),
      sum(Differences)/len(Differences),
      sum(Percentages)/len(Percentages),
      sep="\t")
