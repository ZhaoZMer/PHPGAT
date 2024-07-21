import os
import pickle as pkl
import numpy as np
from Bio import SeqIO

os.makedirs('node_feature/', exist_ok=True)

# def node2id(file_in_fn):
#     file_list = os.listdir(file_in_fn)
#     file2idx = {}
#     for idx, file in enumerate(file_list):
#         file2idx[file.rsplit('.', 1)[0]] = idx
#     return file2idx
#
#
# virus2id = node2id('train_phage/')
# pkl.dump(virus2id, open('node_feature/phages.dict', 'wb'))
#
# prokaryote2id = node2id('prokaryote/')
# pkl.dump(prokaryote2id, open('node_feature/prokaryote.dict', 'wb'))
#
# test_virus2id = node2id('single_contig/')
# pkl.dump(test_virus2id, open('node_feature/test_phages.dict', 'wb'))

def return_4mer(file_in_fn):
    # alphbet
    k_list = ["A", "C", "G", "T"]
    nucl_list = ["A", "C", "G", "T"]
    for i in range(3):
        tmp = []
        for item in nucl_list:
            for nucl in k_list:
                tmp.append(nucl+item)
        k_list = tmp
    mer2dict = {mer: idx for idx, mer in enumerate(k_list)}
    file_list = os.listdir(file_in_fn)
    num_file = len(file_list)
    file2idx = {}
    feature = np.zeros((num_file, 256))  #创建一个全零矩阵
    for idx, file in enumerate(file_list):  #遍历文件列表，同时给每个文件一个唯一的索引
        file2idx[file.rsplit('.', 1)[0]] = idx  #使用file2idx字典存储每个文件的名字和它对应的索引
        for record in SeqIO.parse(file_in_fn + file, 'fasta'):  #读取Fasta文件
            seq = str(record.seq)  #转换为字符串
            seq = seq.upper()  #转换为大写，以确保碱基的一致性
            for pos in range(len(seq)-3):   #遍历每个序列，以4-mer为窗口
                try:
                    feature[idx][mer2dict[seq[pos:pos+4]]] += 1  #若找到了该索引，则在该位置的计数加1；若pos=0，指从第1个元素开始取到第4个（[0,4]）
                except:
                    #print(seq[pos:pos+4])
                    pass    #若碱基序列不在字典中，则跳过
    # nomarlization  归一化
    norm_feature = np.zeros((num_file, 256))
    for i in range(len(feature)):
        norm_feature[i] = (feature[i] - np.min(feature[i]))/(np.max(feature[i]) - np.min(feature[i]))  #归一化公式：（X-Xmin）/（Xmax-Xmin） np.min(feature[i])指第i行中的最小值
    return norm_feature, file2idx

phages, phages2id = return_4mer('train_phage/')
pkl.dump(phages2id, open('node_feature/phages.dict', 'wb'))  #保存位置索引
pkl.dump(phages, open('node_feature/phages.F', 'wb'))   #保存特征矩阵

prokaryote, prokaryote2id = return_4mer('prokaryote/')
pkl.dump(prokaryote2id, open('node_feature/prokaryote.dict', 'wb'))
pkl.dump(prokaryote, open('node_feature/prokaryote.F', 'wb'))

test_phages, test_phages2id = return_4mer('single_contig/')
pkl.dump(test_phages2id, open('node_feature/test_phages.dict', 'wb'))
pkl.dump(test_phages, open('node_feature/test_phages.F', 'wb'))