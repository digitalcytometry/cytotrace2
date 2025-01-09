import torch
import torch.nn as nn
import torch.autograd as autograd
import math
import torch.nn.functional as F

class STEFunction(autograd.Function):
    @staticmethod
    def forward(ctx, input):
        return (input > 0).float()

    @staticmethod
    def backward(ctx, grad_output):
        return F.hardtanh(grad_output)


class StraightThroughEstimator(nn.Module):
    def __init__(self):
        super(StraightThroughEstimator, self).__init__()

    def forward(self, x):
        x = STEFunction.apply(x)
        return x


class Binary_module(nn.Module):
    def __init__(self, input_size: int = 14271, hidden_size: int = 24, dropout: float = 0.5,
                *args, **kwargs):
        super(Binary_module, self).__init__()

        self.input_size=input_size
        self.hidden_size=hidden_size
        self.drop=dropout
        self.dropout = nn.Dropout(self.drop)
        self.maxrank = nn.Parameter(torch.FloatTensor(1,1).uniform_(0,1))
        self.maxrank_param = 1000
        self.weight = nn.Parameter(torch.FloatTensor(self.input_size, self.hidden_size).normal_(mean=-0.1,std=0.055))
        self.ste = StraightThroughEstimator()
        num_enrich = 2
        num_labels = 1
        self.batchnorm = nn.BatchNorm1d(self.hidden_size*num_enrich, affine=False)
        self.out = nn.Linear(self.hidden_size*num_enrich, num_labels)

    def forward(self, x_rank, x_log2, B):
        # UCell calculation using x_rank
        W = self.ste(self.weight)
        n = W.sum(0).unsqueeze(0)
        maxrank = n.max()+10+torch.clamp(self.maxrank,min=0)*self.maxrank_param

        x_rank = torch.minimum(x_rank,maxrank)
        R = torch.matmul(x_rank,W)
        R_UCell = 1 - ((R - (n*(n+1))/2)/(n*maxrank))

        # AMS calculation using x_log2
        gs_backgrounds = (B @ W) 
        bg_score = ((x_log2 @ gs_backgrounds) / gs_backgrounds.sum(0))
        raw_score = (x_log2 @ W) / W.sum(0)
        R_AMS = raw_score - bg_score

        # Concatenated gene set scores
        R_all = torch.cat((R_UCell, R_AMS),1)

        # Batch normalization to transfer scores into shared space
        R_out_norm = self.batchnorm(R_all)

        # Apply gene set level dropout
        R_out_norm = self.dropout(R_out_norm)

        # Generate prediction
        pred = self.out(R_out_norm)

        return pred

class BinaryEncoder(nn.Module):

    def __init__(self, num_layers=6, **block_args):
        super().__init__()
        self.ste = StraightThroughEstimator()
        self.layers = nn.ModuleList([Binary_module(**block_args) for i in range(num_layers)])

    def forward(self, x_rank, x_log2, B):
        pred_list = []
        for l in self.layers:
            pred = l(x_rank, x_log2, B)
            pred_list.append(pred.unsqueeze(1))
        pred_out = torch.cat(pred_list,1)

        return pred_out


