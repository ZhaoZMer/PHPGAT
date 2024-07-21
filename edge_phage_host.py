import os
import subprocess
import glob
import pickle as pkl
import pandas as pd
from Bio import SeqIO
from utils import check_folder


prokaryote_fn = "prokaryote/"
blast_tab_out = "blast_tab/"
Knowledge_graph = "Cyber_data/"
crt_out_fn = "out/crt_out/"

################################################################################
############################  Check the folder #################################
################################################################################
check_folder(blast_tab_out)
check_folder(crt_out_fn)
check_folder("crispr_db/")


_ = subprocess.check_call("cat dataset/phages_train.fasta single_contig/* > out/query.fa", shell=True)

################################################################################
###############################  Run CRISRP   ##################################
################################################################################

print("\n\n" + "{:-^80}".format("Run CRISPR"))

try:
    genome_list = os.listdir(prokaryote_fn)
    for genome in genome_list:
        make_crt_cmd = 'java -cp CRT1.2-CLI/.:bin crt ' + prokaryote_fn + genome + ' ' + crt_out_fn + genome.split(".")[0] + '.out'
        _ = subprocess.check_call(make_crt_cmd, shell=True)
except:
    print("CRISPR Error")


for filename in glob.glob(os.path.join(crt_out_fn, '*.out')):
    with open(filename, 'r') as f:
        if 'No CRISPR' in f.read():
            os.remove(filename)

for filename in glob.glob(os.path.join(crt_out_fn, '*.out')):
    with open(filename, 'r') as f:
        lines = f.readlines()
        repeat_base = [line.split()[2] for line in lines if line[0].isdigit() and len(line.split()) >= 3]

    repeat_file_name = filename.replace('.out', '.fa')

    with open(repeat_file_name, 'w') as f:
        for i, repeat in enumerate(repeat_base):
            name = repeat_file_name.split('/')[2]
            print(name)
            name = name.split('.')[0]
            f.write('>' + name + '.' + 'spacer_locus' + str(i + 1) + '_' + str(i + 1) + '\n')
            f.write(repeat + '\n')

_ = subprocess.check_call("rm out/crt_out/*.out", shell=True)



_ = subprocess.check_call("cat out/crt_out/*.fa > out/crt_out/allCRISPR.fasta", shell=True)

make_blast_cmd = 'makeblastdb -in out/crt_out/allCRISPR.fasta -dbtype nucl -parse_seqids -out crispr_db/allCRISPRs'
subprocess.check_call(make_blast_cmd, shell=True)

crispr_output_file = "out/crispr_out.tab"
crispr_call_cmd = 'blastn -query out/test.fa -db crispr_db/allCRISPRs -outfmt 6 -out out/crispr_out.tab -num_threads 16 -evalue 1e-10 -gapopen 10 -penalty -1 -gapextend 2 -word_size 7 -dust no -task megablast -perc_identity 90'

_ = subprocess.check_call(crispr_call_cmd, shell=True)

crispr_pred = {}
with open(crispr_output_file) as file_out:
    for line in file_out.readlines():
        parse = line.replace("\n", "").split("\t")
        phages = parse[0]
        prokaryote = parse[1].split('|')[1]
        prokaryote = prokaryote.split('.')[0]
        if phages not in crispr_pred:
            crispr_pred[phages] = []
            crispr_pred[phages].append(prokaryote)
        elif phages in crispr_pred:
            crispr_pred[phages].append(prokaryote)

for key in crispr_pred:
    crispr_pred[key] = list(set(crispr_pred[key]))

with open("out/crispr_pred.ntw", 'w') as file_out:
    for phage, prokaryotes in crispr_pred.items():
        for prokaryote in prokaryotes:
            file_out.write(phage + "," + prokaryote + "\n")

pkl.dump(crispr_pred, open('out/crispr_pred.dict', 'wb'))

################################################################################
################################  Run RBP   ####################################
################################################################################

print("\n\n" + "{:-^80}".format("Run RBP"))
try:
    cmd = 'cp out/all_testProteins.fasta data/sequences.fasta'
    _ = subprocess.check_call(cmd, shell=True)
    cmd = 'python PhageRBPdetect_v3_inference.py'
    _ = subprocess.check_call(cmd, shell=True)
except:
    print("PhageRBPdetect Error")

RBPpredict_df = pd.read_csv("out/data/predictions.csv")
RBP_csv_df = RBPpredict_df[RBPpredict_df.iloc[:, 1] == 1]
RBP_values = RBP_csv_df.iloc[:, 0].values.tolist()

with open('out/RBPID.txt', 'w', encoding='utf-8') as f:
    for value in RBP_values:
        f.write("%s\n" % value)

with open('out/RBPID.txt', 'r') as f:
    accessions_to_extract = set(line.strip() for line in f)

RBP_dict = {}

with open('out/data/sequences.fasta', 'r') as f:
    for line in f:
        if line.startswith('>'):
            accession = line.strip()[1:]
            accession = accession.split(" ")[0]
            sequence = ''
        else:
            sequence += line.strip()
        if accession:
            RBP_dict[accession] = sequence

with open('out/RBP_sequences.fasta', 'w') as f:
    for accession in accessions_to_extract:
        if accession in RBP_dict:
            f.write(f'>{accession}\n')
            f.write(RBP_dict[accession] + '\n')

make_diamond_cmd = 'diamond makedb --threads 8 --in out/all_host_proteins.fa -d out/host_database.dmnd'
_ = subprocess.check_call(make_diamond_cmd,shell=True)
diamond_cmd = 'diamond blastp --threads 8 --sensitive -d out/host_database.dmnd -q out/RBP_sequences.fasta -o out/RBPdata-diamond.tab --outfmt 6'
_ = subprocess.check_call(diamond_cmd, shell=True)
RBP_diamond_out_fp = "out/RBPdata-diamond.tab"

RBP_pred = {}
with open(RBP_diamond_out_fp) as file_out:
    for line in file_out.readlines():
        parse = line.replace("\n", "").split("\t")
        phages = parse[0].rsplit("_", 1)[0]
        prokaryote = parse[1].rsplit("_", 1)[0]
        if phages not in RBP_pred:
            RBP_pred[phages] = []
            RBP_pred[phages].append(prokaryote)
        elif phages in RBP_pred:
            RBP_pred[phages].append(prokaryote)

for key in RBP_pred:
    RBP_pred[key] = list(set(RBP_pred[key]))

with open("out/RBP_pred.ntw", 'w') as file_out:
    for phage, prokaryotes in RBP_pred.items():
        for prokaryote in prokaryotes:
            file_out.write(phage + "," + prokaryote + "\n")
pkl.dump(RBP_pred, open('out/RBP_pred.dict', 'wb'))

################################################################################
##############################  Run 16S_rRNA   #################################
################################################################################
print("\n\n" + "{:-^80}".format("Run 16S_rRNA"))

blast_cmd = 'blastn -query out/test.fa -db rRNAblastout/all16S_rRNA -outfmt 6 -out out/rRNA_pred.tab -num_threads 16 -evalue 1e-10 -gapopen 10 -penalty -1 -gapextend 2 -word_size 7 -dust no -task megablast -perc_identity 90'
_ = subprocess.check_call(blast_cmd, shell=True)

rRNA_pred = {}
with open('out/rRNA_pred.tab') as file_out:
    for line in file_out.readlines():
        parse = line.replace("\n", "").split("\t")
        phages = parse[0].rsplit("_", 1)[0]
        prokaryote = parse[1].rsplit("_", 1)[0]
        if phages not in rRNA_pred:
            rRNA_pred[phages] = []
            rRNA_pred[phages].append(prokaryote)
        elif phages in rRNA_pred:
            rRNA_pred[phages].append(prokaryote)

for key in rRNA_pred:
    rRNA_pred[key] = list(set(rRNA_pred[key]))

with open("out/rRNA_pred.ntw", 'w') as file_out:
    for phage, prokaryotes in rRNA_pred.items():
        for prokaryote in prokaryotes:
            file_out.write(phage + "," + prokaryote + "\n")
pkl.dump(rRNA_pred, open('out/rRNA_pred.dict', 'wb'))

################################################################################
###################################  Run BLASTN  ###############################
################################################################################

genome_list = os.listdir(prokaryote_fn)
for genome in genome_list:
    blast_cmd = 'blastn -query out/query.fa -db blast_db/'+genome.split(".")[0]+' -outfmt 6 -out ' + blast_tab_out + genome.split(".")[0]+'.tab -num_threads 16 -evalue 1e-10 -gapopen 10 -penalty -1 -gapextend 2 -word_size 7 -dust no -task megablast -perc_identity 90'
    print("Running blastn...")
    _ = subprocess.check_call(blast_cmd, shell=True)

tab_file_list = os.listdir(blast_tab_out)
prokaryote2phages = {}
for file in tab_file_list:
    prokaryote_id = file.split('.')[0]
    with open(blast_tab_out+file) as file_in:
        for line in file_in.readlines():
            tmp = line.split('\t')
            phages_id = tmp[0]
            if prokaryote_id not in prokaryote2phages:
                prokaryote2phages[prokaryote_id] = []
                prokaryote2phages[prokaryote_id].append(phages_id)
            elif prokaryote_id in prokaryote2phages:
                prokaryote2phages[prokaryote_id].append(phages_id)

for key in prokaryote2phages:
    prokaryote2phages[key] = list(set(prokaryote2phages[key]))

with open("out/phage_host.ntw", 'w') as file_out:
    for prokaryote , phages in prokaryote2phages.items():
        for phage in phages:
            _ = file_out.write(prokaryote + "," + phage + "\n")
