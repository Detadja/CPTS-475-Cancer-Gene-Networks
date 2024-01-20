#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from igraph import Graph, mean
from igraph import VertexSeq
from igraph import EdgeSeq
from igraph import summary
from igraph import plot
from igraph import GraphBase
import igraph
from igraph import VertexClustering
from igraph import clustering
import time
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as sch
#import louvain



#This function returns a sanitized protein with spaces and digits removed, and case set to lower		
def sanitize(string):
	string=re.sub(r'[^\w]','',string)	
	string=re.sub(r'[\d_]','',string)		
	string=string.lower()
	return string

#returns a list as a string with sep separators.
def listtostring(list, sep):
	string=""
	if len(list)>0:
		string+=list[0]#.decode('utf-8')
	for item in list[1:]:
		string+=sep+item#.decode('utf-8')
	return string

#function to read fasta files
def readFasta(file, keepFile="NA"):
	import re
	sequences=[]
	kp=set()
	with open(file,'r') as input:
		if keepFile!="NA":
			with open(keepFile,'r') as keep:
				for line in keep:
					kp.add(line[:-1])
		#print kp
		for line in input:
			if line[0]!=">":
				if len(kp)>0:
					if sequences[-1][0] in kp:
						#print "keep"
						sequences[-1][1]+=sanitize(line)
					else:
						#print sequences[-1], "not on keep list"
						continue
				else:
					sequences[-1][1]+=sanitize(line)
			else:
				if len(sequences)>0 and sequences[-1][1]=="":
					del sequences[-1]
				sequences.append( [re.sub('[^0-9]', '',line[1:]),""] )
	print(len(sequences), sequences[0])
	return sequences

#Edit Distance
def levenshteinDistance(s1,s2):
		
    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2,char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            #scoremod=matrix[(char1,char2)]-4 #Value between 0 and -8
            if char1 == char2:
                newDistances.append(distances[index1])
                #newDistances.append(distances[index1]-scoremod)
            else:
                newDistances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             newDistances[-1])))
                #newDistances.append(1 + min((distances[index1],
                #                             distances[index1+1],
                #                             newDistances[-1])) - scoremod)
        distances = newDistances
#   print distances[-1]
    return distances[-1]
#END

#bounded edit distance    
def levenshteinDistanceThresh(s1,s2,k):
#    print "LD:", k
#    print s1
#    print s2

    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62

    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2,char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            if abs(index1-index2) > k:
                newDistances.append(sys.maxsize)
            elif char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],distances[index1+1],newDistances[-1])))	
        distances = newDistances

        if min(distances)>k:
            return sys.maxsize		
    return distances[-1]

#DiWANN graph construction
def makeGraphEfficiently(species, repeatmatrix, filename, sparse=False):
	
	G = Graph(directed=True)
	G.add_vertices(len(species))
	
	#summary (G)
	
	names=[]
	for r in labelss:
		names.append('; '.join(r))
	#G.vs["label"] = names
	
###
	for row in	 range(len(species)):
		#populate graph
		#print len(repeatmatrix[row])
		if sparse==True:
			import scipy
			m=min(list(scipy.sparse.find(repeatmatrix.tocsc()[:,row])[2]) + list(scipy.sparse.find(repeatmatrix.tocsc()[row,:])[2]))
			edgeTo = getSparseMinED(repeatmatrix, row, m)
			#return
		else:
			m=min(l for l in repeatmatrix[row] if l > 0)
			edgeTo = getMinED(repeatmatrix, row, m)
		#print edgeTo
		for e in edgeTo:
			G.add_edge(e[0], e[1], weight=m)
		
	#summary(G)
	G.simplify(combine_edges=min)
	summary(G)
	G.write(filename+"DWRNgraphEff.gml","gml")
	#summary(G)
	return G

#threshold based network construction
def makeGraph(species, RepeatMatrix, filename):

	G_new = Graph(directed=False)
	G_new.add_vertices(len(species))
	
	G = Graph()
	G.add_vertices(len(species))
	
	names=[]
	for r in labelss:
		names.append('; '.join(r))

	#G_new.vs["label"] = names
	#G.vs["label"] = names
	
	its=len(species)*len(species)
	it=0
	for i in range(len(species)):
		# print RepeatMatrix[i]
		minED = min(l for l in RepeatMatrix[i] if l > 0)
		
		for j in range(len(species)):
			#print "adding"
			if (RepeatMatrix[i][j] == minED):
				G_new.add_edge(i,j,weight=RepeatMatrix[i][j])
				#weights.append(RepeatMatrix[i][j])
			if i>j and RepeatMatrix[i][j]<2000000:
				G.add_edge(i,j,weight=RepeatMatrix[i][j])
			#print "added"
			it+=1
		print(round((it*1.0)/its,3))
	summary(G_new)
	summary(G)
	
	
	G_new.write(filename+"DWN-base.gml","gml")
	G.write(filename+"Thresh.gml","gml")
 
	return G

#finds lowest in a vector
def getMinED(repeatmatrix, row, minVal):
#	minVal=sys.maxint
	
	#print repeatmatrix
	ret=[]
	
	for col in range(len(repeatmatrix[row])):
#		#print "a"
#		if repeatmatrix[row][col] < minVal:
#			minVal = repeatmatrix[row][col]
#			ret = [(row,col)]
#		elif repeatmatrix[row][col] == minVal:
		if repeatmatrix[row][col] == minVal:
			ret.append((row,col))
		if repeatmatrix[row][col] < minVal and repeatmatrix[row][col] < 0:
			print("Warning, min passed was wrong")
			ret = [(row,col)]
			minVal=repeatmatrix[row][col] 
	
	#print minVal
	#print ret
	return ret
	
def getSparseMinED(repeatmatrix, row, minVal):
#	minVal=sys.maxint
	
	#print repeatmatrix
	ret=[]
	
	for col in range(len(repeatmatrix[row])):
#		#print "a"
#		if repeatmatrix[row][col] < minVal:
#			minVal = repeatmatrix[row][col]
#			ret = [(row,col)]
#		elif repeatmatrix[row][col] == minVal:
		if repeatmatrix[row,col] == minVal:
			ret.append((row,col))
		if repeatmatrix[row,col] < minVal and repeatmatrix[row,col] < 0:
			print("Warning, min passed was wrong")
			ret = [(row,col)]
			minVal=repeatmatrix[row,col] 
	
	#print minVal
	#print ret
	return ret

#creats DiWANN similarity matrix
def createMinMatrix(species, useBoundedLD=False):
    #print("started creatMinMatrix")
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio import pairwise2
    matrix = matlist.blosum62
	#print(matrix)
    skippedCells=0	
    repeatmatrix=[]
    its=(len(species)*len(species))/2
    it=len(species)
    #print("its=",its)
    #print("it=",it)
    for row in range(len(species)):
        #Add new row to matrix
        repeatmatrix.append([])
        repeatmatrix[row].extend([sys.maxsize]*len(species))
    #print(repeatmatrix)
    for row in range(len(species)):
        #print "row:", row
        if row>0:
            print("Percent complete:", round((it*1.0)/its,3))
            #calculate ranges
            maxED=[]
            minED=[]
            if row>1:
                rowMin=min(repeatmatrix[row][0:row-1])
                #print repeatmatrix[row][0:row-1]
                #print rowMin
            else:
                rowMin=sys.maxsize
			
            for col in range(row+1,len(species)):
                minED.append(abs(repeatmatrix[0][col]-repeatmatrix[0][row]))
                maxED.append(repeatmatrix[0][col]+repeatmatrix[0][row])
			#then only compare possible mins
			#print row, len(minED)
			
            if len(maxED)>0:
                lowestMax = min(maxED)
            else:
                lowestMax = sys.maxsize
                #print lowestMax
			
            for col in range(row+1,len(species)):
                it+=1
                colMin = min(repeatmatrix[col])
                #print colMin
                #if col - (row+1) < 0:
                #print "ERROR"
                if (minED[col-(row+1)] > lowestMax or minED[col-(row+1)] > rowMin) and (minED[col-(row+1)] > colMin): 
                    repeatmatrix[row][col]=sys.maxsize
                    repeatmatrix[col][row]=sys.maxsize
					#print "skipping ", row, col
                    skippedCells+=1
                else:
                    if useBoundedLD == True:
                        repeatmatrix[row][col]=levenshteinDistanceThresh(species[row][0],species[col][0], max(colMin,rowMin))
					 #repeatmatrix[row][col]=pairwise2.align.globaldx(species.repeats[row].sequence.upper(), species.repeats[col].sequence.upper(), matrix, score_only=True, one_alignment_only=True)
                        if repeatmatrix[row][col] > max(colMin,rowMin):
                            repeatmatrix[row][col]=sys.maxsize
                        #print max(colMin,rowMin)
                    else:
                        repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
                        #repeatmatrix[row][col]=pairwise2.align.globaldx(species.repeats[row].sequence.upper(), species.repeats[col].sequence.upper(), matrix, score_only=True, one_alignment_only=True)
                    repeatmatrix[col][row]=repeatmatrix[row][col]
                    #if repeatmatrix[row][col] == 0:
                        #print "WARNING!!!"
                    if repeatmatrix[row][col] < rowMin:
                        rowMin=repeatmatrix[row][col]

        else:
            for col in range(len(species)):
                repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
                #repeatmatrix[row][col]=pairwise2.align.globaldx(species.repeats[row].sequence.upper(), species.repeats[col].sequence.upper(), matrix, score_only=True, one_alignment_only=True)
                #print repeatmatrix[col][row]
                repeatmatrix[col][row]=repeatmatrix[row][col]
            repeatmatrix[0][0] = sys.maxsize

    #print("Skipped computing", skippedCells, "of", ((len(species)*len(species))-len(species))/2, "cells")
    return repeatmatrix

#Creats similarity matric for threshold-based methods
def createRepeatMatrix(species, thresh=2):
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio import pairwise2
    matrix = matlist.blosum62

    repeatmatrix=[]
    for row in range(len(species)):
        #Add new row to matrix
        repeatmatrix.append([ listtostring(labelss[row],";") ])
        repeatmatrix[row].extend([0]*len(species))
    its=len(species)*len(species)/2
    i=0
    for row in range(len(species)):
        for col in range(row,len(species)):
            #repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
            repeatmatrix[row][col]=levenshteinDistanceThresh(species[row][0],species[col][0],thresh)
            #repeatmatrix[row][col]=pairwise2.align.globaldx(species[row][0].upper(), species[col][0].upper(), matrix, score_only=True, one_alignment_only=True)
            repeatmatrix[col][row]=repeatmatrix[row][col]
            
            i+=1
        print("Percent complete:", round((i*1.0)/its,3))

    return repeatmatrix


def MatrixtoCSV(mat,file):
	with open(file, 'wb') as csv:
		csv.write(u'\ufeff'.encode('utf-8'))
		for row in mat:
			for col in row:
				csv.write(str(col)+",")
			csv.write("\n")

def create_dendrogram(df_matrix):
    hc=AgglomerativeClustering(n_clusters=7,affinity='precomputed',linkage='complete')
    y_hc=hc.fit_predict(df_matrix)
    linkage_matrix=sch.linkage(df_matrix,'complete')
    #visualize clusters
    plt.figure(figsize=(35, 15))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    sch.dendrogram(linkage_matrix,leaf_rotation=90.,leaf_font_size=8.,labels=labelss)
    plt.show()
    plt.savefig("Dendrogram01.jpg")
    return y_hc

def makeNetwork(filename,thresh,bthresh,type="", species=None,BL=""):
    graph=""
    if type=="PB":
        print("Prune + bound*****************")
        print(filename)
        start = time.time()
        C = createMinMatrix(species, False)
        C = np.asarray(C)
        end = time.time()
       # y_hc=create_dendrogram(C)
        print("Pruning and Bounding:", (end-start))
        #MatrixtoCSV(C,"DWRN"+filename+".csv")
        graph = makeGraphEfficiently(species, C, filename+"minNet")
        summary(graph)
    elif type=="BFT":
        print("Brute force*****************")
        #thr=len(species[0][0])/3
        thr=thresh
        print(filename)#, "up to", thr
        start = time.time()
        C = createRepeatMatrix(species, thr)
        #C = createRepeatMatrix(species)
        end = time.time()
        print("Brute Force (my implementation):", (end-start))
        #MatrixtoCSV(C,"BF"+filename+".csv")
        graph = makeGraph(species, C, filename+"fullNet")
        graph.simplify(combine_edges=min)
        summary(graph)
    elif type=="Blast":
        print("Approximate distances*****************")
        print(filename, BL)
        import csv
        nodes={}
        graph=Graph()
        sum=0
        with open(BL, "r") as file:
            reader=csv.reader(file,delimiter="\t")
            for line in reader:
                if line[0] not in nodes:
                    nodes[line[0]]=len(nodes)
                    #print(line[0])
                    #print(nodes)
                    graph.add_vertex(name=line[0])
                if line[1] not in nodes:
                    nodes[line[1]]=len(nodes)
                    #print(line[1])
                    graph.add_vertex(name=line[1])
                #if float(line[10])<1.0e-10:
                #if float(line[2])>92.0:
                if float(line[11])>bthresh:
					
                    graph.add_edge(nodes[line[0]],nodes[line[1]],weight=float(line[11]))
                    sum+=1
            #summary(graph)
            graph=graph.simplify(combine_edges=min)
            summary(graph)
            graph.write(filename+"BlastGraph.gml","gml")
    return graph

style = {}
style["edge_curved"] = False
style["margin"] = 50
style["edge_arrow_size"]=0.1
style["vertex.label.cex"]=0.4
style["vertex.label.family"]="Helvetica"
style["vetrex.label.color"]='black'
style["edge_width"]=0.5
style["edge_color"]="black"
style["vertex_size"]=9
style["vertex.label.font"]=1.5

USA_unique_T10 = pd.read_csv('USA_T10_2021_final.csv',header=0)

uni_seq_T10 = USA_unique_T10.iloc[:,1]
uni_seq_T10 = uni_seq_T10.values.tolist()
uni_labels_T10 = USA_unique_T10.iloc[:,0]
uni_labels_T10 = uni_labels_T10.values.tolist()
uni_acc_T10 = USA_unique_T10.iloc[:,2]
uni_acc_T10 = uni_acc_T10.values.tolist()
uni_date_T10 = USA_unique_T10.iloc[:,3]
uni_date_T10 = uni_date_T10.values.tolist()
uni_count_T10 = USA_unique_T10.iloc[:,4]
uni_count_T10 = uni_count_T10.values.tolist()

time_seq = USA_unique_T10.iloc[:,[1]]
time_seq = time_seq.values.tolist()
df = pd.DataFrame()
df.insert(0, 'label', (range(len(uni_seq_T10))))
df = df.applymap(str)
x_labels = df.iloc[:,[0]]
labelss=x_labels.values.tolist()
rfile="USATime10"

#create DiWANN network
diwann_usa_time10=Graph()				
diwann_usa_time10=makeNetwork(rfile,thresh="",bthresh="",type="PB",species=time_seq)
plot(diwann_usa_time10,'diwann_USA_Time10.png',**style)

#Visualization example using a list color_mapT10 containing node colors according to the state label
'''
bstyle = {}
bstyle["edge_curved"] = False
bstyle["margin"] = 50
bstyle["edge_arrow_size"]=0.1
bstyle["vertex.label.cex"]=0.2
bstyle["vertex.label.family"]="Helvetica"
bstyle["vetrex.label.color"]='black'
bstyle["edge_width"]=0.2
#bstyle["edge_color"]=black
size=[]
#seq_label=[]
for i in range(len(uni_count_T10)):
    if uni_count_T10[i] == 1:
        size.append(7)
    if uni_count_T10[i] > 1 and uni_count_T10[i] <= 10:
        size.append(7)
    if uni_count_T10[i] > 10 and uni_count_T10[i] <= 60:
        size.append(8)
    if uni_count_T10[i] > 60 and uni_count_T10[i] <= 500:
        size.append(10)
    if uni_count_T10[i] > 500:
        size.append(14)
        #seq_label.append(uni_labels_AllT[i])
shape=[]
for a in aa10_95:
    if a == 'I':
        shape.append('triangle-down')
    if a == 'T':
        shape.append('circle')
    if a == 'X':
        shape.append('rectangle')
bstyle["vertex_shape"] = shape
bstyle["vertex_size"]=size
bstyle["vertex.label.font"]=0.9
bstyle["vertex_color"] = color_mapT10
'''








