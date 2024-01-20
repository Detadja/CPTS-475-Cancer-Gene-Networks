# CPTS-475-Cancer-Gene-Networks

* This is a school project created in Python and R, by me and my partner, Pei-Yau Weng, in the doctorate program (Cpts 575). It aims to take in csv files with a series of genes related to specific types of cancers (given to us by our Project Mentor, Shruti Sunil Patil), `donor_art_seq_468.csv`, `icgc_468_gene_tissue.csv`, `icgc_data_468.csv`, and `USA_T10_2021_final.csv`. Then clean, sort, identify and analyze key genes that have high likelyhoods of being directly involved with that specific cancer.

* It does this via using a Sequential Similarity Network (SSN) to find the genes which have high relation to each other, and constructing a Protein-Protein Interaction Network (PPIN). To do this we used the Directed Weighted All Nearest Neighbours (DiWANN) Network<sup>1</sup> to construct the SSN, and once we got the results, we directly compared our results with the information on the Human Protein Atlas (HPA)<sup>2</sup>. As such, it outputs the networks made using the information gained from the algorithm. (Unfortunately, the code which contained the cleaning of the original dataset and created the SSN is inaccessible and I have no further contact with my previous partner).

* Last updated on: 1st Dec 2022, as the Term Project for the course Cpts 475 - Data Science.

* This program was written in R. To run, import files to local system and run on R Studio. The Python files are inaccessible and hence could not be uploaded.

* [Reference 1](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2453-2#Abs1)
  [Reference 2](www.proteinatlas.org)
