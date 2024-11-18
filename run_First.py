import numpy as np
import pandas as pd
import os
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import subprocess
import argparse
import re
import utils
from utils import check_folder

#####################################################################
##########################  Input Params  ###########################
#####################################################################

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--contigs', type=str, default='test.fasta', help=' the file name of the input file')
parser.add_argument('--len', type=int, default=8000, help='the minimum input length limit, that is, predict only for sequence >= len bp (default 0)')
parser.add_argument('--gpus', type=int, default=0, help='GPU device ID (only used if --use-cpu is not set)')
parser.add_argument('--model', type=str, default='pretrain', choices=['pretrain', 'retrain'], help='predicting host with pretrained parameters or retrained paramters (default pretrain)')
parser.add_argument('--lr', default=0.001, type=float, help='Learning rate of training')
parser.add_argument('--weight-decay', default=0, type=float, help='Weight Decay of optimizer')
parser.add_argument('--topk',  type=int, default=1, help='topK')
parser.add_argument('--t',  type=float, default=0.98, help='The confident threshold for predicting phages, the higier the threshold the higher the precision. (default 0.98)')
parser.add_argument('--epochs', type=int, default=4000, help='Epochs to train')
parser.add_argument('--batch-size', type=int, default=512, help='Batch size of Training')
parser.add_argument('--use-cpu', action='store_true', help='Force the use of CPU (default: False)')

inputs = parser.parse_args()

check_folder("input")
check_folder("pred")
check_folder("Split_files")
check_folder("tmp_pred")
check_folder("train_phage/")

#####################################################################
#########################   processing    ###########################
#####################################################################

for record in SeqIO.parse('dataset/phages_train.fasta', 'fasta'):
    _ = SeqIO.write(record, 'train_phage/'+record.id, 'fasta')

#####################################################################
#############################  Start  ###############################
#####################################################################

cnt = 0
file_id = 0
records = []
for record in SeqIO.parse(inputs.contigs, 'fasta'):
    seq = str(record.seq).upper()
    record.seq = Seq(seq)
    if len(record.seq) > inputs.len:
        records.append(record)
        cnt += 1
    if cnt != 0 and cnt % 2000 == 0:
        SeqIO.write(records, f"Split_files/contig_{file_id}.fasta", "fasta")
        records = []
        file_id += 1
        cnt = 0
SeqIO.write(records, f"Split_files/contig_{file_id}.fasta","fasta")
file_id += 1

for i in range(file_id):

    cmd = f"mv Split_files/contig_{i}.fasta input/"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"Moving file Error")
        exit()

    cmd = "python edge_phage_phage.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"phage_phage Error")
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        exit()

    cmd = f"python edge_host_host.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"host_host Error")
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        exit()

    cmd = f"python edge_phage_host.py --gpus {inputs.gpus}"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"phage_host Error")
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        exit()

    cmd = f"python create_feature.py"
    try:
        check_folder("node_feature")
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"create feature Error")
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        exit()

    cmd = f"python multimodal_graph.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print(f"multimodal Graph Error for file contig_{i}")
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        exit()

    cmd = f"python train.py --model {inputs.model} --lr {inputs.lr} --weight-decay {inputs.weight_decay} --gpus {inputs.gpus}  --topk {inputs.topk} --t {inputs.t} --epochs {inputs.epochs} --batch-size {inputs.batch_size} "
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Error for file contig_{i}")
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        exit()

