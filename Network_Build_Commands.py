#!/usr/bin/python3

#packages used to build the three networks
import camoco as co
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


#Brachipodium Ref
Bd21 = co.RefGen("Bd21")
co.RefGen.from_gff("BdistachyonBd21-3v1.1.gene.gff3" ,"Bd21", "BdistachyonBd21-3v1.1", "BdistachyonBd21-3v1.1", "Brachipodium Ref")

#Brachipodium GO propogation
co.GOnt.from_obo('go.obo',
    'Brach.gene.GO.txt',
    'Bd21GO',
    'Brachi Gene Ontology',
     Bd21)

#Build Brachy network
co.COB.from_table('Bd21_Meesh_V3.csv',
'Bd21_Treated',
'Bd21_Treated Samples with mock on day two',
Bd21,
rawtype='RNASEQ',
max_gene_missing_data=0.4,
max_accession_missing_data=0.4,
min_single_sample_expr=1,
min_expr=0.001,
quantile=False,
max_val=300,
sep=','
)



#Wheat refernce generation
WheatRef = co.RefGen("WheatRef")
co.RefGen.from_gff("iwgsc_refseqv1.0_HighConf_2017Mar13-transposed.gff3" ,"WheatRef", "Chinese_Spring_v1.0", "TChinese_Spring_v1.0", "Wheat Ref")

#Wheat GO propogation
co.GOnt.from_obo('go.obo',
    'WheatGOnt_CSW.txt',
    'WheatGO',
    'Wheat Gene Ontology',
     WheatRef)


#Build Wheat networks
co.COB.from_table('Wheat_W2691_V2.csv',
'W2691',
'Wheat_W2691 Samples with mock on day two',
WheatRef,
rawtype='RNASEQ',
max_gene_missing_data=0.4,
max_accession_missing_data=0.4,
min_single_sample_expr=1,
min_expr=0.001,
quantile=False,
max_val=300,
sep=','
)

co.COB.from_table('Wheat_Sr9b_V2.csv',
'Sr9b',
'Wheat Sr9b Samples with mock on day two',
WheatRef,
rawtype='RNASEQ',
max_gene_missing_data=0.4,
max_accession_missing_data=0.4,
min_single_sample_expr=1,
min_expr=0.001,
quantile=False,
max_val=300,
sep=','
)
