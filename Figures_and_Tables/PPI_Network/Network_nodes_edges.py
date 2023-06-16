import pandas as pd 

import networkx as nx 

from matplotlib import pylab as plt 

 

file_nodes = '/lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI/20230609genes.tsv' 

file_edges = '/lab02/Data_Raw/Xiaolong/HumanCommon/20210711ConsensusPathDB.StringDB.GIANTv2/20210824.ConsensusPathDB_v35.StringDB_v11.5.GIANTv2_combined.3DB.combined.tsv.gz' 

 

df_nodes = pd.read_csv(file_nodes, sep='\t', index_col=0)#2413 

 

# add column "ADHD_linkage" 

genelists = ['target', 'ADHD_evidence', 'other_Neuro', 'SFARI'] 

dc_genes = {} 

for keyword in genelists: 

    dc_genes[keyword] = set(df_nodes[df_nodes[keyword] > 0].index) 

    print(keyword,len(dc_genes[keyword])) 

 

# target 270 

# ADHD_evidence 17 

# other_Neuro 6 

# SFARI 45 

 

 

nodes_all = set(df_nodes.index) 

df_edges = pd.read_csv(file_edges, sep='\t')#656252 

df_edges = df_edges[df_edges['node1'].isin(nodes_all) & df_edges['node2'].isin(nodes_all)]#134 

df_edges = df_edges.groupby(['node1','node2']).apply(lambda x:x.iloc[0]).reset_index(drop=True)#134 

# add edge type 

def getEdgeType(r): 

    edgetype = [] 

    if r['STRING_combined_score']: 

        edgetype.append('S') 

    if r['consensus_score'] != 0: 

        edgetype.append('C') 

    if r['GIANT_score'] != 0: 

        edgetype.append('G') 

    edgetype = ''.join(edgetype) 

    if len(edgetype) == 1: 

        return edgetype 

    return 'M' 

df_edges['edgeType'] = df_edges.apply(getEdgeType,axis=1)#134 

 

target_genes = nodes_all 

tdf_edges = df_edges 

# create a Graph with target_genes, plus intermediate nodes supported by edgeType M 

G = nx.Graph() 

# create network 

G.add_nodes_from(target_genes) 

for row, r in tdf_edges.iterrows(): 

    G.add_edge(r['node1'], r['node2'], weight=r['DataBaseCount']) 

 

ls_networks = list(nx.connected_components(G))# 163 independent subgraph. one with 50 genes 

ls_networks = sorted(ls_networks, key=len, reverse=True) 

print(len(ls_networks), [len(i) for i in ls_networks])#163 [50, 22, 11, 8, 5, 4, 3, 3, 3, 2, ...] 

 

nodes_use = ls_networks[0] 

G0 = G.subgraph(nodes_use).copy() 

nodes_use = ls_networks[1] 

G1 = G.subgraph(nodes_use).copy() 

nodes_use = ls_networks[2] 

G2 = G.subgraph(nodes_use).copy() 

nodes_use = ls_networks[3] 

G3 = G.subgraph(nodes_use).copy() 

 

 

tG4 = G.subgraph([i for j in ls_networks[:4] for i in j]) 

print(sorted(list(G0.degree), key=lambda x:x[1], reverse=True)[:10]) 

#[('DLGAP2', 9), ('LRRK2', 7), ('RBFOX1', 7), ('CYP2A6', 6), ('SNTG1', 5), ('CLTC', 5), ('IQGAP1', 5), ('RIMS2', 4), ('SLC6A12', 3), ('ALDH1A1', 3)] 

 

# get node position 

pos0 = nx.kamada_kawai_layout(G0, center=[-2,2],scale=3) 

pos1 = nx.kamada_kawai_layout(G1, center=[2,2],scale=2) 

pos2 = nx.kamada_kawai_layout(G2, center=[-2,-2],scale=1) 

pos3 = nx.kamada_kawai_layout(G3, center=[2,-2],scale=1) 

 

def getNodesEdgesOfGraph(net_G, net_pos): 

    nodes_net = set(net_G.nodes) 

    tdf_nodes_net = df_nodes.loc[list(net_G.nodes)] 

    tdf_nodes_net['x'] = [net_pos[i][0] for i in tdf_nodes_net.index] 

    tdf_nodes_net['y'] = [net_pos[i][1] for i in tdf_nodes_net.index] 

    edges_net = set(net_G.edges) 

    tdf_edges_net = df_edges[df_edges.apply(lambda x:((x['node1'], x['node2']) in edges_net) or ((x['node2'], x['node1']) in edges_net), axis=1)].copy() 

    tdf_edges_net['x'] = tdf_edges_net['node1'].apply(lambda x:net_pos[x][0]) 

    tdf_edges_net['y'] = tdf_edges_net['node1'].apply(lambda x:net_pos[x][1]) 

    tdf_edges_net['xend'] = tdf_edges_net['node2'].apply(lambda x:net_pos[x][0]) 

    tdf_edges_net['yend'] = tdf_edges_net['node2'].apply(lambda x:net_pos[x][1]) 

    return tdf_nodes_net, tdf_edges_net 

 

ls_pos = [(G0,pos0)] 

results = [getNodesEdgesOfGraph(net_G, net_pos) for net_G, net_pos in ls_pos] 

tdf_nodes_final = pd.concat([i[0] for i in results]) 

tdf_edges_final = pd.concat([i[1] for i in results]) 

tdf_nodes_final = tdf_nodes_final.fillna(0) 


# save tdf_edges and tdf_nodes 

HEADER='20230609_truncated.networkx' 

tdf_nodes_final.to_csv(f'/lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI//{HEADER}.nodes.tsv',sep='\t') 

tdf_edges_final.to_csv(f'/lab01/Projects/Rohan_Projects/CNV_Project/2023/Submission_Pipeline/Results/Plots/Network_PPI//{HEADER}.edges.tsv',sep='\t',index=None) 

