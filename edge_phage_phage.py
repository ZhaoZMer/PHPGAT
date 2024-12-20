import os
import subprocess
import pandas as pd
import pickle as pkl
import networkx as nx
from Bio import SeqIO
from utils import check_folder, make_protein_clusters_mcl, build_clusters, build_pc_matrices, create_network, to_clusterer

out_fn = "out/"
contig_in = "input/"
contig_out = "single_contig/"
file_in_fn = "single_contig/"
file_out_fnn = "all_proteins/"
Knowledge_graph = "Cyber_data/"

check_folder(out_fn)
check_folder(file_in_fn)
check_folder(file_out_fnn)
check_folder("out/phage_db")

################################################################################
############################  Rename the files  ################################
################################################################################

file_list = sorted(os.listdir(contig_in)) 
seq = []
old_file_id = 0
contig_id = 0
with open("name_list.csv",'w') as list_out:
    list_out.write("contig_name,idx\n") 
    for file_n in file_list: 
        for record in SeqIO.parse(contig_in+file_n, "fasta"):
            name = "newphage_"+str(old_file_id) + "_" + str(contig_id)
            _ = list_out.write(record.id + "," + name + "\n")
            record.id = "newphage_"+str(old_file_id) + "_" + str(contig_id)
            _ = SeqIO.write(record, contig_out+name+".fasta", "fasta")
            contig_id += 1
        old_file_id += 1

################################################################################
###################### Translate contigs into 6 ORFs ###########################
################################################################################

file_list = os.listdir(file_in_fn)
for file in file_list:
    prodigal_cmd = 'prodigal -i ' + file_in_fn + file + ' -a ' + file_out_fnn + file +' -f gff -p meta'
    print("Running prodigal...")
    _ = subprocess.check_call(prodigal_cmd, shell=True)

all_testprotein_f = out_fn + "all_testProteins.fasta"
_ = subprocess.check_call("cat {0} > {1}".format(file_out_fnn+"*",  all_testprotein_f), shell=True)

all_protein_f = out_fn + "all_proteins.fa"
_ = subprocess.check_call("cat {0} {1} > {2}".format(all_testprotein_f, "dataset/protein.fasta", all_protein_f), shell=True)


################################################################################
############################## Run diamond BLASTp  #############################
################################################################################

print("\n\n" + "{:-^80}".format("Diamond BLASTp"))
print("Creating Diamond database and running Diamond...")

try:
    make_diamond_cmd = 'diamond makedb --threads 8 --in out/all_proteins.fa -d out/database.dmnd'
    print("Creating Diamond database...")
    _ = subprocess.check_call(make_diamond_cmd, shell=True)

    diamond_cmd = 'diamond blastp --threads 8 --sensitive -d out/database.dmnd -q out/all_proteins.fa -o out/database.self-diamond.tab'
    print("Running Diamond...")
    _ = subprocess.check_call(diamond_cmd, shell=True)
    with open("out/database.self-diamond.tab", 'r') as infile:
        filtered_lines = [line for line in infile if float(line.split()[2]) >= 90]
    with open("out/database.diamond.tab", 'w') as outputfile:
        outputfile.writelines(filtered_lines)
    diamond_out_fp = "out/database.diamond.tab"
    database_abc_fp = "out/database.self-diamond.tab.abc"
    _ = subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fp, database_abc_fp),shell=True)
except:
    print("create database failed")
    exit(1)

blastp = pd.read_csv(database_abc_fp, sep=' ', names=["contig", "ref", "e-value"])
protein_id = sorted(list(set(blastp["contig"].values) | set(blastp["ref"].values)))
contig_protein = [item for item in protein_id if "newphage" == item.split("_")[0]]
contig_id = [item.rsplit("_", 1)[0] for item in contig_protein]
description = ["hypothetical protein" for item in contig_protein]
gene2genome = pd.DataFrame({"protein_id": contig_protein, "contig_id": contig_id,"keywords": description})
gene2genome.to_csv(out_fn + "contig_gene_to_genome.csv", index=None)

_ = subprocess.check_call("cat dataset/database_gene_to_genome.csv {0}contig_gene_to_genome.csv > {1}gene_to_genome.csv".format(out_fn, out_fn), shell=True)

################################################################################
################################## Run MCL #####################################
################################################################################

print("\n\n" + "{:-^80}".format("Protein clustering"))
print("Loading proteins...")
gene2genome_fp = out_fn+"gene_to_genome.csv"
gene2genome_df = pd.read_csv(gene2genome_fp, sep=',', header=0)

pc_overlap, pc_penalty, pc_haircut, pc_inflation = 0.8, 2.0, 0.1, 2.0
pcs_fp = make_protein_clusters_mcl(database_abc_fp, out_fn, pc_inflation)
print("Building the cluster and profiles (this may take some time...)")

protein_df, clusters_df, profiles_df, contigs_df = build_clusters(pcs_fp, gene2genome_df)
print("Saving files")
dfs = [gene2genome_df, contigs_df, clusters_df]
names = ['proteins', 'contigs', 'pcs']
output_dir = out_fn

for name, df in zip(names, dfs):
    fn = "Cyber_{}.csv".format(name)
    fp = os.path.join(output_dir, fn)
    index_id = name.strip('s') + '_id'
    if not os.path.exists(fp):
        df.set_index(index_id).to_csv(fp)
    else:
        print("File {} exists and will be used. Use -f to overwrite.".format(fn))

profiles_fn = "Cyber_profiles.csv"
profiles_fp = os.path.join(out_fn, profiles_fn)
if not os.path.exists(profiles_fp):
    profiles_df.to_csv(profiles_fp, index=False)
else:
    print("File {} exists and will be used. Use -f to overwrite.".format(profiles_fn))

contigs_df = pd.read_csv("out/Cyber_contigs.csv")
clusters_df = pd.read_csv("out/Cyber_pcs.csv")
profiles_df = pd.read_csv("out/Cyber_profiles.csv")

contigs_csv_df = contigs_df.copy()
print("Read {} entries from {}".format(len(contigs_csv_df), os.path.join(output_dir, '{}_contigs.csv'.format(name))))
contigs_csv_df.index.name = "pos"
contigs_csv_df.reset_index(inplace=True)

pcs_csv_df = clusters_df.copy()
profiles = profiles_df.copy()

before_filter = len(profiles)
cont_by_pc = profiles.groupby("pc_id").count().contig_id.reset_index()

cont_by_pc.columns = ["pc_id", "nb_proteins"]
pcs_csv_df = pd.merge(pcs_csv_df, cont_by_pc, left_on="pc_id", right_on="pc_id", how="left")
pcs_csv_df.fillna({"nb_proteins": 0}, inplace=True)

pcs_csv_df = pcs_csv_df[pcs_csv_df['nb_proteins'] > 1]
at_least_a_cont = cont_by_pc[cont_by_pc['nb_proteins'] > 1]
profiles = profiles[profiles['pc_id'].isin(at_least_a_cont.pc_id)]
print("Read {} entries (dropped {} singletons) from {}".format(len(profiles), (before_filter - len(profiles)), profiles_fp))
pcs_csv_df = pcs_csv_df.reset_index(drop=True)
pcs_csv_df.index.name = "pos"
pcs_csv_df = pcs_csv_df.reset_index()

matrix, singletons = build_pc_matrices(profiles, contigs_csv_df, pcs_csv_df)
profiles_csv = {"matrix": matrix, "singletons": singletons}
merged_df = contigs_csv_df
merged_fp = os.path.join(output_dir, 'merged_df.csv')
merged_df.to_csv(merged_fp)

ntw = create_network(matrix, singletons, thres=1, max_sig=300)
fi = to_clusterer(ntw, out_fn+"intermediate.ntw", merged_df.copy())

################################################################################
################################# Run BLASTN ###################################
################################################################################

_ = subprocess.check_call("cat single_contig/* > out/test.fa", shell=True)
_ = subprocess.check_call("cat out/test.fa dataset/phages_trainfasta > out/allphage.fasta", shell=True)

make_blast_cmd = 'makeblastdb -in out/allphage.fasta -dbtype nucl -parse_seqids -out out/phage_db/allPHAGE'
subprocess.check_call(make_blast_cmd, shell=True)

output_file = "out/phage_out.tab"
phage_call_cmd = 'blastn -query out/test.fa -db out/phage_db/allPHAGE -outfmt 6 -out out/phage_out.tab -num_threads 16 -evalue 1e-10 -gapopen 10 -penalty -1 -gapextend 2 -word_size 7 -dust no -task megablast -perc_identity 90'
print("Running BLASTN...")
_ = subprocess.check_call(phage_call_cmd, shell=True)


phage_pred = {}
with open(output_file) as file_out:
    for line in file_out.readlines():
        parse = line.replace("\n", "").split("\t")
        phage = parse[0]
        ref_phage = parse[1]
        if phage not in phage_pred:
            phage_pred[phage] = []
            phage_pred[phage].append(ref_phage)
        if phage in phage_pred:
            phage_pred[phage].append(ref_phage)


for key in phage_pred:
    phage_pred[key] = list(set(phage_pred[key]))

with open("out/phage_pred.ntw", 'w') as file_out:
    for phage, ref_phages in phage_pred.items():
        for ref_phage in ref_phages:
            file_out.write(phage + "," + ref_phage + "\n")

################################################################################
############################### Dump the graph #################################
################################################################################

with open("out/intermediate.ntw", 'r') as file_in:
    with open("out/phage_phage.ntw", 'w') as file_out:
        for line in file_in:
            tmp = line.strip().split(" ")
            node1 = tmp[0]
            node2 = tmp[1]
            file_out.write(node1 + "," + node2 + "\n")
