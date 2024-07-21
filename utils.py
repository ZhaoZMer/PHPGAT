import os
import sys
import subprocess
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sparse
from torch.utils.data import Dataset


def check_folder(file_name):
    if not os.path.exists(file_name):  
        _ = os.makedirs(file_name)   
    else:                             
        print("folder {0} exist... cleaning dictionary".format(file_name))
        if os.listdir(file_name): 
            try:   
                _ = subprocess.check_call(f"rm -rf {file_name}", shell=True) 
                _ = os.makedirs(file_name)  
                print("Dictionary cleaned")
            except:    
                print("Cannot clean your folder... permission denied")
                exit(1)  

def load_mcl_clusters(fi):
    with open(fi) as f:
        c = [line.rstrip("\n").split("\t") for line in f]
    c = [x for x in c if len(c) > 1]
    nb_clusters = len(c)
    formatter = "PC_{{:>0{}}}".format(int(round(np.log10(nb_clusters))+1))
    name = [formatter.format(str(i)) for i in range(nb_clusters)]
    size = [len(i) for i in c]
    clusters_df = pd.DataFrame({"size": size, "pc_id": name}).set_index("pc_id")
    return clusters_df, name, c

def make_protein_clusters_mcl(abc_fp, out_p, inflation=2):
    print("Running MCL...")
    abc_fn = "merged"
    mci_fn = '{}.mci'.format(abc_fn)
    mci_fp = os.path.join(out_p, mci_fn)
    mcxload_fn = '{}_mcxload.tab'.format(abc_fn)
    mcxload_fp = os.path.join(out_p, mcxload_fn)
    subprocess.check_call("mcxload -abc {0} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {1} -write-tab {2}".format(abc_fp, mci_fp, mcxload_fp), shell=True)
    mcl_clstr_fn = "{0}_mcl{1}.clusters".format(abc_fn, int(inflation*10))
    mcl_clstr_fp = os.path.join(out_p, mcl_clstr_fn)
    subprocess.check_call("mcl {0} -I {1} -use-tab {2} -o {3}".format(mci_fp, inflation, mcxload_fp, mcl_clstr_fp), shell=True)
    return mcl_clstr_fp

def build_clusters(fp, gene2genome):
    clusters_df, name, c = load_mcl_clusters(fp)
    print("Using MCL to generate PCs.")
    gene2genome.set_index("protein_id", inplace=True)
    for prots, clust in zip(c, name):
        try:
            gene2genome.loc[prots, "cluster"] = clust
        except KeyError:
            prots_in = [p for p in prots if p in gene2genome.index]
            not_in = frozenset(prots) - frozenset(prots_in)
            print("{} protein(s) without contig: {}".format(len(not_in), not_in))
            gene2genome.loc[prots_in, "cluster"] = clust
    for clust, prots in gene2genome.groupby("cluster"):
        clusters_df.loc[clust, "annotated"] = prots.keywords.count()
        if prots.keywords.count():
            keys = ";".join(prots.keywords.dropna().values).split(";")
            key_count = {}
            for k in keys:
                k = k.strip()
                try:
                    key_count[k] += 1
                except KeyError:
                    key_count[k] = 1
            clusters_df.loc[clust, "keys"] = "; ".join(["{} ({})".format(x, y) for x, y in key_count.items()])
    gene2genome.reset_index(inplace=True)
    clusters_df.reset_index(inplace=True)
    profiles_df = gene2genome.loc[:, ["contig_id", "cluster"]].drop_duplicates()
    profiles_df.columns = ["contig_id", "pc_id"]
    contigs_df = pd.DataFrame(gene2genome.fillna(0).groupby("contig_id").count().protein_id)
    contigs_df.index.name = "contig_id"
    contigs_df.columns = ["proteins"]
    contigs_df.reset_index(inplace=True)
    return gene2genome, clusters_df, profiles_df, contigs_df

def build_pc_matrices(profiles, contigs, pcs):
    pc_by_cont = profiles.groupby("contig_id").count().pc_id
    pc_by_cont = pd.merge(contigs.sort_values("pos").loc[:, ["pos", "contig_id", "proteins"]], pc_by_cont.to_frame(), how="left", left_on="contig_id", right_on="contig_id").fillna(0)
    singletons = (pc_by_cont.proteins - pc_by_cont.pc_id).values
    singletons = sparse.lil_matrix(singletons).transpose()
    profiles.index.name = "pos"
    profiles.reset_index(inplace=True)
    profiles = pd.merge(profiles, pcs.loc[:, ["pc_id", "pos"]], left_on="pc_id", right_on="pc_id", how="inner", suffixes=["", "_pc"])
    profiles = pd.merge(profiles, contigs.loc[:, ["contig_id", "pos"]], left_on="contig_id", right_on="contig_id", how="inner", suffixes=["", "_contig"])
    profiles = profiles.loc[:, ["pos_contig", "pos_pc"]]
    matrix = sparse.coo_matrix(([1]*len(profiles), (zip(*profiles.values))), shape=(len(contigs), len(pcs)), dtype="bool")
    return matrix.tocsr(), singletons.tocsr()

def create_network(matrix, singletons, thres=1, max_sig=1000):
    contigs, pcs = matrix.shape
    pcs += singletons.sum()
    T = 0.5 * contigs * (contigs - 1)
    logT = np.log10(T)
    number_of_pc = matrix.sum(1) + singletons
    number_of_pc = number_of_pc.A1
    commons_pc = matrix.dot(sparse.csr_matrix(matrix.transpose(), dtype=int))
    S = sparse.lil_matrix((contigs, contigs))
    total_c = float(commons_pc.getnnz())
    i = 0
    for A, B in zip(*commons_pc.nonzero()):
        if A != B:
            a, b = sorted([number_of_pc[A], number_of_pc[B]])
            pval = stats.hypergeom.sf(commons_pc[A, B] - 1, pcs, a, b)
            sig = min(max_sig, np.nan_to_num(-np.log10(pval) - logT))
            if sig > thres:
                S[min(A, B), max(A, B)] = sig
            i += 1
            if i % 1000 == 0:
                sys.stdout.write(".")
            if i % 10000 == 0:
                sys.stdout.write("{:6.2%} {}/{}\n".format(i / total_c, i, total_c))
    S += S.T
    S = S.tocsr()
    if len(S.data) != 0:
        print("Hypergeometric contig-similarity network:\n {0:10} contigs,\n {1:10} edges (min:{2:.2} max:{3:.2}, threshold was {4})".format(contigs, S.getnnz(), S.data.min(), S.data.max(), thres))
    else:
        raise ValueError("No edge in the similarity network !")
    return S

def to_clusterer(matrix, fi, contigs=None,names=None):
    names = contigs if names is None else names
    names = names.set_index("pos").contig_id
    with open(fi, "wt") as f:
        matrix = sparse.dok_matrix(matrix)
        for r, c in zip(*matrix.nonzero()):
            f.write(" ".join([str(x) for x in (names[r], names[c], matrix[r, c])]))
            f.write("\n")
    print("Saving network in file {0} ({1} lines).".format(fi, matrix.getnnz()))
    return fi

class PosEdge(Dataset):
    def __init__(self, edge_index_pos):
        self.edge_index_pos = edge_index_pos

    def __len__(self):
        return self.edge_index_pos.shape[1]

    def __getitem__(self, idx):
        return self.edge_index_pos[:, [idx]]

class AverageMeter(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count