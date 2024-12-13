import torch
import torch.nn as nn
from torch.nn import functional as F
from torch_geometric.nn import GATv2Conv


class DotDecoder(torch.nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x1, x2):
        out = (x1*x2).sum(dim=-1)
        return out.reshape(-1)

#########################################################################
######################## GATV2 #########################################
#########################################################################
class GATV2_Encoder(torch.nn.Module):
    def __init__(self, num_node, num_feat, device):
        super().__init__()
        self.device = device
        self.gat1 = GATv2Conv(256, 64, 4, concat=False)
        self.ELU = nn.ELU()
        self.emb = torch.nn.Embedding(num_node, num_feat)
        torch.nn.init.xavier_normal_(self.emb.weight)


    def forward(self, x, adj_t, perturbed=False):
        all_embeddings = []
        x_out = self.emb.weight
        x_out = self.gat1(x_out, adj_t)
        x_out = self.ELU(x_out)
        if (perturbed):
            random_noise = torch.rand_like(x_out, device=self.device)
            x_out = x_out + torch.sign(x_out) * F.normalize(random_noise, dim=-1) * self.eps
        all_embeddings.append(x_out)

        all_embeddings = torch.stack(all_embeddings, dim=1)
        all_embeddings = torch.mean(all_embeddings, dim=1)

        return all_embeddings
