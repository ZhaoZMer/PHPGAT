import os
import subprocess
from utils import check_folder

check_folder("out/rRNA/")
check_folder("out/rRNA_out")
check_folder("rRNAblastout/")
check_folder("blast_db")
check_folder("host_host_blast_out/")
host_fn = "prokaryote/"
################################################################################
##############################  Run 16S_rRNA   #################################
################################################################################

print("\n\n" + "{:-^80}".format("Run 16S_rRNA"))
genome_list = os.listdir(host_fn)
for genome in genome_list:
    cmd = 'barrnap prokaryote/' + genome + ' -o out/rRNA/' + genome.split('.')[0] + '.fasta'
    _ = subprocess.check_call(cmd, shell=True)

cmd = 'rm prokaryote/*.fai'
_ = subprocess.check_call(cmd,shell=Ture)

rRNA_fn = "out/rRNA/"
rRNA_fn_out = "out/rRNA_out/"
genome_16S_rRNA_list = os.listdir(rRNA_fn)
for genome in genome_16S_rRNA_list:
    rRNA_dict = {}
    current_accession = None
    current_sequence = ''

    with open(rRNA_fn + genome, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_accession:
                    rRNA_dict[current_accession] = current_sequence
                    current_sequence = ''

                accession_full = line.strip()
                if "16S_rRNA" in accession_full:
                    accession_parts = accession_full[8:]
                    current_accession = accession_parts.split(":")[2]
                else:
                    current_accession = None
            else:
                if current_accession:
                    current_sequence += line.strip()

    if current_accession and current_sequence:
        rRNA_dict[current_accession] = current_sequence

    i = 0
    with open(rRNA_fn_out + genome.split('.')[0] + '.fasta', 'w') as f:
        for accession, sequence in rRNA_dict.items():
            if sequence:
                f.write(f'>{genome.split(".")[0] + "." + str(i)}\n')
                f.write(sequence + '\n')
                i += 1

_ = subprocess.check_call('cat out/rRNA_out/* > out/rRNA.fasta', shell=True)

make_blast_cmd = 'makeblastdb -in out/rRNA.fasta -dbtype nucl -parse_seqids -out rRNAblastout/all16S_rRNA'
_ = subprocess.check_call(make_blast_cmd, shell=True)
blast_cmd = 'blastn -query out/rRNA.fasta -db rRNAblastout/all16S_rRNA -outfmt 6 -out out/rRNAHHI.tab -num_threads 16'
_ = subprocess.check_call(blast_cmd, shell=True)

rRNA_pred = {}
with open('out/rRNAHHI.tab') as file_out:
    for line in file_out.readlines():
        parse = line.replace("\n", "").split("\t")
        prokaryote1 = parse[0].rsplit(".", 1)[0]
        prokaryote2 = parse[1].rsplit(".", 1)[0]
        ident = float(parse[2])
        pre = float(parse[-1])
        e = float(parse[-2])
        if prokaryote1 not in rRNA_pred and ident > 97 and e < 0.001 and pre > 97:
            rRNA_pred[prokaryote1] = []
            rRNA_pred[prokaryote1].append(prokaryote2)
        elif prokaryote1 in rRNA_pred and ident > 97 and e < 0.001 and pre > 97:
            rRNA_pred[prokaryote1].append(prokaryote2)

for key in rRNA_pred:
    rRNA_pred[key] = list(set(rRNA_pred[key]))

with open("out/rRNAHHI.ntw", 'w') as file_out:
    for prokaryote1 in rRNA_pred:
        for prokaryote2 in rRNA_pred[prokaryote1]:
            if prokaryote1 != prokaryote2:
                file_out.write(prokaryote1 + "," + prokaryote2 + "\n")

################################################################################
#################################  Run BLASTN  #################################
################################################################################
host_host_blast_tab_out = "host_host_blast_out/"

genome_list = os.listdir(host_fn)

for genome in genome_list:
    make_blast_cmd = 'makeblastdb -in ' + host_fn + genome + ' -dbtype nucl -parse_seqids -out blast_db/' + genome.split(".")[0]
    print("Creating blast database...")
    _ = subprocess.check_call(make_blast_cmd, shell=True)

for genome in genome_list:
    for genome2 in genome_list:
        blast_cmd = 'blastn -query prokaryote/' + genome2 + ' -db blast_db/' + genome.split(".")[0] + ' -outfmt 6 -out ' + host_host_blast_tab_out + genome.split(".")[0] + '.tab -num_threads 16 -evalue 1e-10 -gapopen 10 -penalty -1 -gapextend 2 -dust no -perc_identity 90'
        print("Running blastn...")
        _ = subprocess.check_call(blast_cmd, shell=True)

################################################################################
################################  host-host   #################################
################################################################################
tab_file_list = os.listdir(host_host_blast_tab_out)
host2host = {}
for file in tab_file_list:
    host_id = file.split('.')[0] 
    host2_id_list = []
    with open(host_host_blast_tab_out + file) as file_in:
        for line in file_in.readlines():
            tmp = line.split('\t') 
            host2_id = tmp[0]  
            if host_id not in host2host:
                host2host[host_id] = []
                host2host[host_id].append(host2_id)  
            elif host_id in host2host:
                host2host[host_id] = [host2_id] 

for key in host2host:
    host2host[key] = list(set(host2host[key]))

with open("out/blastNHHI.ntw", 'w') as file_out:
    for host in host2host: 
        for host2 in host2host[host]: 
            _ = file_out.write(host + "," + host2 + "\n")

_ = subprocess.check_call("cat out/16S_rRNAHHI.ntw out/blastNHHI.ntw > out/host_host.ntw")