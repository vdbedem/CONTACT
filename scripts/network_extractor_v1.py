import networkx as nx, re
import matplotlib.pyplot as plt
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout
layout = graphviz_layout

import sys
"""
usage:  python network_extractor_v1.py cypa_detailed_SC25.txt

ARGV[1] should be name of detailed pathways output

Default settings are:
*take all pathways (even length 2) - but filter on relief threshold 0.9 CHANGED FROM 0.95
*take networks length 3 or more

*output pymol instructions and pngs of subgraphs
"""

#EXTRACT PATHWAYS AS LISTS OF LISTS
def extract_pathways(filename, relief_threshold=0.9, filter_residues=0, filter_residues_list=[], min_pathway_length=1):
    """
outputs lists of:
    * list of lists: all_paths, each list has residue numbers connected in an individual pathway
    * list of lists: min_length_paths , which are only include pathways of length greater than min_pathway_length
    * dictionary: residue_counter keys are residues and value is number of pathways that residue is in, does not take into account min_pathway_length **or** relief threshold

any pathway with a single relief value below relief_threshold not included in output list of pathways

when filter_residues flag = 1, only pathways containing residues in filter_residues_list are output to all_paths and min_length_paths

if you define filter_residue_list to a network, then len(all_paths) gives you number of pathways in that network
    """
    f = open(filename)
    new_path = 0
    all_paths = []
    min_length_paths = []
    residues = []
    reliefs = []
    residue_counter = {}
    for line in f:
        if re.match("PATH",line):
            new_path = 1
        elif re.match("1", line):
            if new_path == 1:
                residue = int(line.split(",")[2][0:-2])
                clash_ratio = float(line.split(",")[4])
                relief_ratio = float(line.split(",")[5])
                if residue in residue_counter:
                    residue_counter[residue] +=1
                else:
                    residue_counter[residue] = 1
                if relief_ratio == 0.0:
                    relief_ratio = 1.0
                residues.append(residue)
                reliefs.append(relief_ratio)
        else:
            #CHECKS WHETHER ANY RELIEFS ARE BELOW THRESHOLD
            threshold = 0
            for relief in reliefs:
                if relief < relief_threshold:
                    threshold = 1
            if threshold:
                pass
            else:
                if filter_residues:
                    filtered = 0
                    for residue in residues:
                        if residue in filter_residues_list:
                            filtered =1
                    if filtered == 1:
                        all_paths.append(residues)
                    #ADDS PATHWAYS TO LIST OF PATHWAYS LENGTH REQUIREMENT
                        if len(residues) > min_pathway_length:
                            min_length_paths.append(residues)
                    else:
                        pass                      
                else:
                    all_paths.append(residues)
                #ADDS PATHWAYS TO LIST OF PATHWAYS LENGTH REQUIREMENT
                    if len(residues) > min_pathway_length:
                        min_length_paths.append(residues)
            new_path = 0
            residues = []
            reliefs = []
#   NOTE IF min_pathway_length=1 or 2, then min_length_paths and all_paths should be the same
    return min_length_paths, all_paths, residue_counter

def build_graph(paths, network_length_threshold = 3):
    """
    input is a list of lists from extract_pathway
    will connect networks from that list and throw out any networks with less than network_length_threshold nodes
    output is a graph with subgraphs representing different networks, edges are weighted by number of pathways connecting nodes
    """
    pathway_lengths = []
    G=nx.Graph()
    for pathway in paths:
        pathway_lengths.append(len(pathway))
        i = 0
        while len(pathway) > i+1:
            #CHECK IF IT EXISTS AND THEN ADD WEIGHTS 
            pathway_exists = 0
            for e in G.edges():
                if pathway[i] in e and pathway[i+1] in e:
                    #G.edge[pathway[i]][pathway[i+1]]['weight'] += 1
                    G.get_edge_data(pathway[i],pathway[i+1])['weight'] += 1
                    pathway_exists = 1
            if not pathway_exists:
                G.add_edge(pathway[i],pathway[i+1], weight=1)
            i +=1 
    connects = nx.connected_components(G)
    for c in sorted(nx.connected_components(G), key=len, reverse=True):
      if len(c) < network_length_threshold:
        c.clear()
    ###FILTER GRAPH TO ONLY NETWORKS > THRESHOLD
#   for network in connects:
#       if len(network) < network_length_threshold:
#           G.remove_node(network[0])
#           G.remove_node(network[1])
    #for e in G.edges():
#         print e, G.get_edge_data(e[0],e[1])
#     plt.hist(pathway_lengths)
#     plt.show()
#     print len(pathway_lengths) # total number of pathways
    return G


#HERE is where you actually run things
r = float(sys.argv[2])

#get just the network and how many paths

min_length_paths, all_paths, residue_counter = extract_pathways(sys.argv[1],relief_threshold=r,filter_residues=0,filter_residues_list=[99,111,113,115,55,57,122,61,63])

# from sets import Set
# 
# for path in all_paths:
#     print path
# print len(all_paths)
# # a = Set(all_paths)
### get pathways and SET RELIEF THRESHOLD HERE
# min_length_paths, all_paths, residue_counter = extract_pathways(sys.argv[1],relief_threshold=r)



#FILTER GRAPH TO ONLY NETWORKS > 3 from all pathways
G = build_graph(all_paths,3)
#make pictures and output networks
colors = ['firebrick','cyan','green','magenta','yellow','orange','purple','grey','blue']
i = 0
C=nx.connected_component_subgraphs(G)
for g in C:
    resis = []
    for resi in g:
        resis.append(str(resi))
    if len(resis) > 2:
        #pos=nx.spring_layout(g,weight=0.1)
        pos=layout(g,prog='twopi')
        nx.draw_networkx_nodes(g,pos,node_size=150,node_color=colors[i])
        nx.draw_networkx_labels(g,pos,font_size=7,font_family='sans-serif')
        for e in g.edges():
#             print e, G.get_edge_data(e[0],e[1])['weight']
            nx.draw_networkx_edges(G,pos,edgelist=[e],width=(G.get_edge_data(e[0],e[1])['weight'])/10000+1,color=colors[i] )#for cypa, remove **0.5, for dhfr include
        print "color %s, resi %s" %(colors[i],"+".join(resis))
        print "create %s_sec, resi %s" %(colors[i],"+".join(resis))
        print "show surface, %s_sec" %(colors[i])
        plt.savefig("%s.png" %colors[i], dpi=600)
        plt.clf()    
        i+=1
        if i > 8:
            i = 0

# OTHER THINGS YOU CAN DO:
# ###COUNTING HOW MANY TIMES EACH RESIDUE OCCURS IN ALL PATHWAYS (even length 2, ignoring relief thresholds)
# for k in residue_counter:
#     print k, residue_counter[k]
# #To get a specific network:    
# scpathways3, scpathways2, residue_counter = extract_pathways(sys.argv[1],0.95, filter_residues=3,filter_residues_list=[99,111,113,115,55,57,122,61,63])
# 
# print len(scpathways2)
