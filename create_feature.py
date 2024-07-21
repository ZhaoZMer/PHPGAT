import os
import pickle as pkl
import numpy as np
from Bio import SeqIO

os.makedirs('node_feature/', exist_ok=True)

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
    feature = np.zeros((num_file, 256))  
    for idx, file in enumerate(file_list):  
        file2idx[file.rsplit('.', 1)[0]] = idx  
        for record in SeqIO.parse(file_in_fn + file, 'fasta'):  
            seq = str(record.seq)  
            seq = seq.upper()
            for pos in range(len(seq)-3):  
                try:
                    feature[idx][mer2dict[seq[pos:pos+4]]] += 1  
                except:
                    #print(seq[pos:pos+4])
                    pass    
    # nomarlization  归一化
    norm_feature = np.zeros((num_file, 256))
    for i in range(len(feature)):
        norm_feature[i] = (feature[i] - np.min(feature[i]))/(np.max(feature[i]) - np.min(feature[i])) 
    return norm_feature, file2idx

phages, phages2id = return_4mer('train_phage/')
pkl.dump(phages2id, open('node_feature/phages.dict', 'wb'))  
pkl.dump(phages, open('node_feature/phages.F', 'wb'))   

prokaryote, prokaryote2id = return_4mer('prokaryote/')
pkl.dump(prokaryote2id, open('node_feature/prokaryote.dict', 'wb'))
pkl.dump(prokaryote, open('node_feature/prokaryote.F', 'wb'))

test_phages, test_phages2id = return_4mer('single_contig/')
pkl.dump(test_phages2id, open('node_feature/test_phages.dict', 'wb'))
pkl.dump(test_phages, open('node_feature/test_phages.F', 'wb'))
