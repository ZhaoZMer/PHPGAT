import os
import pandas as pd
import pickle as pkl
import networkx as nx
import numpy as np
from utils import check_folder

################################################################################
###############################  create folder   ###############################
################################################################################

phage_phage_ntw = "out/phage_phage.ntw"
phage_host_ntw = "out/phage_host.ntw"
host_host_ntw = "out/host_host.ntw"
phage_pred_ntw = "out/phage_pred.ntw"
crispr_pred_ntw = "out/crispr_pred.ntw"
RBP_pred_ntw = "out/RBP_pred.ntw"
rRNA_pred_ntw = "out/rRNA_pred.ntw"

check_folder("GCN_data")

################################################################################
############################  Edge construction  ###############################
################################################################################

G = nx.Graph()
with open(phage_phage_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0].split('.')[0]
        node2 = tmp[1].split('.')[0]
        G.add_edge(node1, node2, weight=1)

with open(host_host_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0].split('.')[0]
        node2 = tmp[1].split('.')[0]
        G.add_edge(node1, node2, weight=1)

with open(phage_host_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0].split('.')[0]
        node2 = tmp[1].split('.')[0]
        G.add_edge(node1, node2, weight=1)

with open(phage_pred_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0].split('.')[0]
        node2 = tmp[1].split('.')[0]
        G.add_edge(node1, node2, weight=1)

with open(crispr_pred_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0].split('.')[0]
        node2 = tmp[1].split('.')[0]
        G.add_edge(node1, node2, weight=1)

with open(RBP_pred_ntw) as file_in:
    for line in file_in.readlines():
        tmp = line[:-1].split(",")
        node1 = tmp[0].split('.')[0]
        node2 = tmp[1].split('.')[0]
        G.add_edge(node1, node2, weight=1)

host_df = pd.read_csv('dataset/prokaryote.csv')
phages_df = pd.read_csv('dataset/phages.csv')

host_list = os.listdir('prokaryote/')
host_list = [name.split('.')[0] for name in host_list]

train_phages_list = os.listdir('train_phage/')
train_phages_list = [name.split('.')[0] for name in train_phages_list]

#species2host = {host_df[host_df['Accession'] == item]['Species'].values[0]: item for item in host_list}
#species2host = {host_df[host_df['Accession'] == item]['Genus'].values[0]: item for item in host_list}
#species2host = {host_df[host_df['Accession'] == item]['Family'].values[0]: item for item in host_list}
#species2host = {host_df[host_df['Accession'] == item]['Order'].values[0]: item for item in host_list}
# crispr_pred = pkl.load(open('out/crispr_pred.dict', 'rb'))
# for phages, host in crispr_pred.items():
#     if host in species2host:
#         G.add_edge(phages, species2host[host])
#
# RBP_pred = pkl.load(open('out/RBP_pred.dict', 'rb'))
# for phages, host in RBP_pred.items():
#     if host in species2host:
#         G.add_edge(phages, species2host[host])

for host in host_list:
    #species = host_df[host_df['Accession'] == host]['Species'].values[0]
    #species = host_df[host_df['Accession'] == host]['Genus'].values[0]
    #species = host_df[host_df['Accession'] == host]['Family'].values[0]
    #species = host_df[host_df['Accession'] == host]['Order'].values[0]
    #phage_list = phages_df[phages_df['Species'] == species]['Accession'].values
    #phage_list = phages_df[phages_df['Genus'] == species]['Accession'].values
    #phage_list = phages_df[phages_df['Family'] == species]['Accession'].values
    #phage_list = phages_df[phages_df['Order'] == species]['Accession'].values
    for phage in phage_list:
        if phage in train_phages_list:
            G.add_edge(host, phage, weight=1)

################################################################################
############################  Nodes construction ###############################
################################################################################

phages2id = pkl.load(open("node_feature/phages.dict", 'rb'))
# phagesF = pkl.load(open("node_feature/phages.F", 'rb'))
prokaryote2id = pkl.load(open("node_feature/prokaryote.dict", 'rb'))
# prokaryoteF = pkl.load(open("node_feature/prokaryote.F", 'rb'))
test_phages2id = pkl.load(open("node_feature/test_phages.dict", 'rb'))
# test_phagesF = pkl.load(open("node_feature/test_phages.F", 'rb'))

# node_feature = []
# for node in G.nodes():
#     if node in prokaryote2id.keys():
#         node_feature.append(prokaryoteF[prokaryote2id[node]])
#     elif node in phages2id.keys():
#         node_feature.append(phagesF[phages2id[node]])
#     elif node in test_phages2id.keys():
#         node_feature.append(test_phagesF[test_phages2id[node]])
#     else:
#         print(f"node error {node}")
#         exit()
#
# node_feature = np.array(node_feature)

################################################################################
############################  Label construction ###############################
################################################################################

crispr_pred = pkl.load(open('out/crispr_pred.dict', 'rb'))
RBP_pred = pkl.load(open('out/RBP_pred.dict', 'rb'))
rRNA_pred = pkl.load(open('out/rRNA_pred.dict', 'rb'))
phages_pred = pkl.load(open('out/phage_pred.dict', 'rb'))
phages_df = pd.read_csv("dataset/phages.csv")
prokaryote_df = pd.read_csv("dataset/prokaryote.csv")

idx = 0
test_id = {}
node2label = {}
cnt = 0
for node in G.nodes():
    if "newphage" in node:
        neighbor_label = []
        for _, neighbor in G.edges(node):
            if neighbor in phages2id.keys():
                #phages_label = phages_df[phages_df['Accession'] == neighbor]['Species'].values[0]
                #phages_label = set(phages_df[phages_df['Accession'] == neighbor]['Species'])
                #phages_label = set(phages_df[phages_df['Accession'] == neighbor]['Genus'])
                #phages_label = set(phages_df[phages_df['Accession'] == neighbor]['Family'])
                phages_label = set(phages_df[phages_df['Accession'] == neighbor]['Order'])
                neighbor_label.extend(phages_label)
            elif neighbor in prokaryote2id.keys():
                #prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == neighbor]['Species'].values[0]
                #prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == neighbor]['Genus'].values[0]
                #prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == neighbor]['Family'].values[0]
                prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == neighbor]['Order'].values[0]
                neighbor_label.append(prokaryote_label)
        if len(set(neighbor_label)) == 1:
            node2label[node] = neighbor_label
            test_id[node] = 1
        elif node in crispr_pred:
            node2label[node] = []
            for value in crispr_pred[node]:
                #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == value]['Species'].values[0])
                #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == value]['Genus'].values[0])
                #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == value]['Family'].values[0])
                #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == value]['Order'].values[0])
            test_id[node] = 1
        elif node in RBP_pred:
            node2label[node] = []
            for value in RBP_pred[node]:
                #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == value]['Species'].values[0])
                #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == value]['Genus'].values[0])
                #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == value]['Family'].values[0])
                node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == value]['Order'].values[0])
            test_id[node] = 1
        elif node in phages_pred:
            node2label[node] = []
            for value in phages_pred[node]:
                #node2label[node] = set(phages_df[phages_df['Accession'] == value]['Species'])
                #node2label[node] = set(phages_df[phages_df['Accession'] == value]['Genus'])
                #node2label[node] = set(phages_df[phages_df['Accession'] == value]['Family'])
                phages_pred_label = set(phages_df[phages_df['Accession'] == value]['Order'])
                node2label[node].extend(phages_pred_label)
            test_id[node] = 1
        else:
            node2label[node] = []
            node2label[node].append('unknown')
            test_id[node] = 2
    elif node in prokaryote2id.keys():
        node2label[node] = []
        #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == node]['Species'].values[0])
        #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == node]['Genus'].values[0])
        #node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == node]['Family'].values[0])
        node2label[node].append(prokaryote_df[prokaryote_df['Accession'] == node]['Order'].values[0])
        test_id[node] = 0
    elif node in phages2id.keys():
        #node2label[node] = set(phages_df[phages_df['Accession'] == node]['Species'])
        #node2label[node] = set(phages_df[phages_df['Accession'] == node]['Genus'])
        #node2label[node] = set(phages_df[phages_df['Accession'] == node]['Family'])
        node2label[node] = set(phages_df[phages_df['Accession'] == node]['Order'])
        test_id[node] = 0
    else:
        print("Error: " + node)
    idx += 1

# for node1 in prokaryote2id.keys():
#     for node2 in prokaryote2id.keys():
#         if (node1 == node2):
#             continue
#         if (node1 not in node2label or node2 not in node2label):
#             continue
#         if (node2label[node1] == node2label[node2]):
#             G.add_edge(node1, node2)

# for sub in nx.connected_components(G):
#     flag = 0
#     for node in sub:
#         if "newphage" not in node:
#             flag = 1
#     if not flag:
#         CRISPR_label = ""
#         CRISPR_cnt = 0
#         RBP_label = ""
#         RBP_cnt = 0
#         rRNA_label = ""
#         rRNA_cnt = 0
#         for node in sub:
#             if node in crispr_pred:
#                 CRISPR_cnt += 1
#                 CRISPR_label = crispr_pred[node]
#             elif node in RBP_pred:
#                 RBP_cnt += 1
#                 RBP_label = RBP_pred[node]
#         if CRISPR_cnt == 1:
#             for node in sub:
#                 node2label[node] = CRISPR_label
#         elif RBP_cnt == 1:
#             for node in sub:
#                 node2label[node] = RBP_label

for sub in nx.connected_components(G):
    sub_label = []
    for node in sub:
        if node in phages2id.keys():
            #phages_label = set(phages_df[phages_df['Accession'] == node]['Species'])
            #phages_label = set(phages_df[phages_df['Accession'] == node]['Genus'])
            #phages_label = set(phages_df[phages_df['Accession'] == node]['Family'])
            phages_label = set(phages_df[phages_df['Accession'] == node]['Order'])
            sub_label.extend(phages_label)
        elif node in prokaryote2id.keys():
            #prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == node]['Species'].values[0]
            #prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == node]['Genus'].values[0]
            #prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == node]['Family'].values[0]
            prokaryote_label = prokaryote_df[prokaryote_df['Accession'] == node]['Order'].values[0]
            sub_label.append(prokaryote_label)
    if set(sub_label) == 1:
        for node in sub:
            node2label[node] = sub_label

for key in node2label:
    node2label[key] = list(set(node2label[key]))

id2node = {idx: node for idx, node in enumerate(G.nodes())}
node2id = {node: idx for idx, node in enumerate(G.nodes())}

adj = nx.adjacency_matrix(G)
pkl.dump(adj, open("GCN_data/graph.list", "wb"))
# pkl.dump(node_feature, open("GCN_data/feature.list", "wb" ))
pkl.dump(node2label, open("GCN_data/node2label.dict", "wb"))
pkl.dump(id2node, open("GCN_data/id2node.dict", "wb"))
pkl.dump(node2id, open("GCN_data/node2id.dict", "wb"))
pkl.dump(test_id, open("GCN_data/test_id.dict", "wb"))