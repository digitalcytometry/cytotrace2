import concurrent.futures
import anndata as ad
import scanpy as sc
import sys
import os
import math
import numpy as np 
import pandas as pd
import torch
import scipy
import random
import warnings
from sklearn.preprocessing import MinMaxScaler, scale
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
import pkg_resources
from cytotrace2_py.common import models


use_dt = True
try:
    import datatable as dt
except ImportError as e:
    use_dt = False
    pass

# functions for one-to-one best match mapping of orthology
def top_hit_human(x):
    r = x.sort_values('%id. target Mouse gene identical to query gene')
    return r.iloc[-1]
def top_hit_mouse(x):
    r = x.sort_values('%id. query gene identical to target Mouse gene')
    return r.iloc[-1]


# data loader
def generic_data_loader(rank_expression, log2_expression, batch_size):
    
    rank_tensor = torch.Tensor(rank_expression)
    log2_tensor = torch.Tensor(log2_expression)
    train_tensor = torch.utils.data.TensorDataset(rank_tensor, log2_tensor)
    data_loader = torch.utils.data.DataLoader(train_tensor, batch_size=batch_size, shuffle=False, drop_last=False)
    
    return data_loader


# one model predictor
def validation(data_loader, B_in, model, device):
    prob_pred_list = []
    order_pred_list = []

    model.eval()
    softmax = torch.nn.Softmax(dim=1)
    order_vector = torch.Tensor(np.arange(6).reshape(6,1)/5).to(device)
    for batch_idx, tensor in enumerate(data_loader):
        X_rank, X_log2 = tensor
        X_rank = X_rank.to(device)
        X_log2 = X_log2.to(device)

        model_output = model(X_rank, X_log2, B_in)
        prob_pred = model_output
        prob_pred = prob_pred.squeeze(2)
        prob_pred = softmax(prob_pred)
        prob_order = torch.matmul(prob_pred,order_vector)
        prob_pred_list.append(prob_pred.detach().cpu().numpy())
        order_pred_list.append(prob_order.squeeze(1).detach().cpu().numpy())
        
    return np.concatenate(prob_pred_list,0), np.concatenate(order_pred_list,0)


# dispersion function
def disp_fn(x):
    if len(np.unique(x)) == 1:
        return 0
    else:
        return np.var(x)/np.mean(x)


# choosing top variable genes
def top_var_genes(log2_data):

    dispersion_index = [disp_fn(log2_data[:, i]) for i in range(log2_data.shape[1])]
    top_col_inds = np.argsort(dispersion_index)[-1000:]
                        
    return top_col_inds


def build_mapping_dict():
    fn_mart_export = pkg_resources.resource_filename("cytotrace2_py", "resources/mart_export.txt")
    human_mapping = pd.read_csv(fn_mart_export,sep='\t').dropna().reset_index()
    fn_features = pkg_resources.resource_filename("cytotrace2_py", "resources/features_model_training.csv")
    features = pd.read_csv(fn_features)['0']
    mapping_unique = human_mapping.groupby('Gene name').apply(top_hit_human)
    mapping_unique = mapping_unique.groupby('Mouse gene name').apply(top_hit_mouse)
    mt_dict = dict(zip(mapping_unique['Gene name'].values,mapping_unique['Mouse gene name'].values))
    fn_alias_list = pkg_resources.resource_filename("cytotrace2_py", "resources/human_alias_list.txt")
    human_mapping_alias_and_previous_symbols = pd.read_csv(fn_alias_list,sep='\t')
    mt_dict_alias_and_previous_symbols = dict(zip(human_mapping_alias_and_previous_symbols['Alias or Previous Gene name'].values,
                                                  human_mapping_alias_and_previous_symbols['Mouse gene name'].values))
    
    for gene in features[~features.isin(mt_dict.values())].values:
        if gene.upper() in mt_dict.keys():
            print(gene)
        mt_dict[gene.upper()] = gene
        
    # **New Addition**: Load mouse alias mapping
    fn_mouse_alias = pkg_resources.resource_filename("cytotrace2_py", "resources/mouse_alias_list.txt")
    mouse_alias_mapping = pd.read_csv(fn_mouse_alias, sep='\t')
    mt_mouse_alias_dict = dict(zip(
        mouse_alias_mapping['alias'].values,
        mouse_alias_mapping['mmgene'].values
    ))
    
    return mt_dict, mt_dict_alias_and_previous_symbols, mt_mouse_alias_dict, features


def load(input_path):
    print("cytotrace2: Loading dataset")
    expression = pd.read_csv(input_path,sep='\t',index_col=0).T # read data
    # expression = expression.loc[:,~expression.columns.duplicated()].copy() #drop duplicate gene names, to avoid random information loss make sure the genes are unique
    if expression.columns.duplicated().any():
        raise ValueError("   Please make sure the gene names are unique.")
    if expression.index.duplicated().any():
        raise ValueError("   Please make sure the cell names are unique.")
    return expression

def read_file(file_path):
    
    if use_dt:
        try:
            file_delim = "," if file_path.lower().endswith(".csv") else "\t"

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=pd.errors.ParserWarning)
                file_data = dt.fread(file_path, header=True)
                colnames = pd.read_csv(file_path, sep=file_delim, nrows=1, index_col=0).columns
                rownames = file_data[:, 0].to_pandas().values.flatten()
                file_data = file_data[:, 1:].to_pandas()
                file_data.index = rownames
                file_data.columns = colnames
                file_data = file_data.astype(float)
                if file_data.columns.duplicated().any():
                    raise ValueError("   Please make sure the gene names are unique.")
                if file_data.index.duplicated().any():
                    raise ValueError("   Please make sure the cell names are unique.")

        except Exception as e:
            print("Error encountered while reading input: {}".format(file_path))
            print("Please make sure that you provided the correct path to the input files.",
                    "The following input file formats are supported:",
                    ".csv with comma ',' as delimiter,",
                    "and .txt or .tsv with tab '\\t' as delimiter.")
            raise

        return file_data.transpose()
    
    else:

        file_data = load(file_path)
        return file_data

    


def preprocess(expression, species, cores_to_use=1):
    gene_names = expression.columns 
    mt_dict, mt_dict_alias_and_previous_symbols, mt_mouse_alias_dict, features = build_mapping_dict()

    # If human, map genes to orthologs
    if species == "human":
        mapped_genes = gene_names.map(mt_dict)
        mapped_genes = mapped_genes.astype(object)
    
        unmapped_genes = {value: index for index, value in enumerate(gene_names) if value in gene_names[mapped_genes.isna()]}
        mapped_genes.values[mapped_genes.isna()] = gene_names[mapped_genes.isna()].map(mt_dict_alias_and_previous_symbols)
        expression.columns = mapped_genes.values
        num_genes_mapped = len([i for i in mapped_genes if i in features.values])
        print("    Mapped "+str(num_genes_mapped)+" input gene names to mouse orthologs")    
        duplicate_genes = expression.columns[expression.columns.duplicated()].values
        duplicate_genes = [i for i in duplicate_genes if i is not np.nan]
        idx = [unmapped_genes[i.upper()] for i in duplicate_genes if i.upper() in unmapped_genes.keys()]
        expression = expression.iloc[:, [j for j, c in enumerate(expression.columns) if j not in idx]]
        
    else:   
        # Handle mouse gene aliases
        mapped_genes = gene_names.to_numpy()
        unmapped_genes = set(mapped_genes) - set(features)
        unmapped_genes_in_alias = unmapped_genes.intersection(set(mt_mouse_alias_dict.keys()))
        valid_unmapped_genes = [gene for gene in unmapped_genes_in_alias if mt_mouse_alias_dict[gene] not in mapped_genes]
        is_valid_unmapped_gene = np.isin(mapped_genes, list(valid_unmapped_genes))
        indices_to_replace = np.where(is_valid_unmapped_gene)[0]

        for idx in indices_to_replace:
            alias = mapped_genes[idx]
            mapped_genes[idx] = mt_mouse_alias_dict.get(alias, alias)
            
        expression.columns = mapped_genes
            
    expression = expression[expression.columns[~expression.columns.isna()]]
    
    # check the number of input genes mapped to model features
    intersection = set(expression.columns).intersection(features)
    print("    "+str(len(intersection))+" input genes are present in the model features.")
    if len(intersection) < 9000:
        warnings.warn("    Please verify the input species is correct.\n    In case of a correct species input, be advised that model performance might be compromised due to gene space differences.")
    expression = pd.DataFrame(index=features).join(expression.T).T
    expression = expression.fillna(0)

    cell_names = expression.index
    gene_names = expression.columns
    adata_X = expression.to_numpy()
    log2_data = np.log2(1000000*adata_X.transpose()/adata_X.sum(1)+1).transpose()
    rank_data = scipy.stats.rankdata(expression.values*-1,axis=1,method='average')

    return cell_names, gene_names, rank_data, log2_data
    
def predict(rank_data, log2_data, B_in, cell_names, model_dir, batch_size = 1000):
    
    device = "cpu"

    all_preds_test = []
    all_order_test = []
    all_models_path = pd.Series(np.array([os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(model_dir)) for f in fn]))
    all_models_path = all_models_path[all_models_path.str.endswith('.pt')]

    data_loader = generic_data_loader(rank_data, log2_data, batch_size)
            
    #print('    Started prediction')
    for model_path in all_models_path:
        pytorch_model = torch.load(model_path, map_location=torch.device('cpu'))
        
        model = models.BinaryEncoder()
        model = model.to(device)
        model.load_state_dict(pytorch_model)
    
        prob_test, order_test = validation(data_loader, B_in, model, device)
    
        all_preds_test.append(prob_test.reshape(-1,6,1))
        all_order_test.append(order_test.reshape(-1,1))

    predicted_order = np.mean(np.concatenate(all_order_test,1),1)
    predicted_potency = np.argmax(np.concatenate(all_preds_test,2).mean(2),1)

    labels = ['Differentiated','Unipotent','Oligopotent','Multipotent','Pluripotent','Totipotent']
    labels_dict = dict(zip([0,1,2,3,4,5],labels))
    predicted_df = pd.DataFrame({'preKNN_CytoTRACE2_Score':predicted_order,'preKNN_CytoTRACE2_Potency':predicted_potency})
    predicted_df['preKNN_CytoTRACE2_Potency'] = predicted_df['preKNN_CytoTRACE2_Potency'].map(labels_dict)

    ## GSBN scores
    predicted_df[labels] = np.concatenate(all_preds_test,2).mean(2)
    predicted_df.index = cell_names

    return predicted_df


# CytoTRACE-based smoothing functions
def get_markov_matrix(log2_data, top_col_inds):
    num_samples, num_genes = log2_data.shape
    
    sub_mat = log2_data[:, top_col_inds]
    with np.errstate(divide="ignore", invalid="ignore"): 
        D = np.corrcoef(sub_mat)  # Pairwise pearson-r corrs

    D[np.arange(num_samples), np.arange(num_samples)] = 0
    D[np.where(D != D)] = 0
    cutoff = max(np.mean(D), 0)
    D[np.where(D < cutoff)] = 0

    A = D / (D.sum(1, keepdims=True) + 1e-5)
    return A

def smooth_subset(chunk_log2_data,chunk_predicted_df,top_col_inds,maxiter):
    markov_mat = get_markov_matrix(chunk_log2_data, top_col_inds)
    score = chunk_predicted_df["preKNN_CytoTRACE2_Score"]
    init_score = score.copy()
    prev_score = score.copy()
    traj = []

    for _ in range(int(maxiter)):
        cur_score = 0.9 * markov_mat.dot(prev_score) + 0.1 * init_score
        traj.append(np.mean(np.abs(cur_score - prev_score)) / (np.mean(init_score) + 1e-6))
        if np.mean(np.abs(cur_score - prev_score)) / (np.mean(init_score) + 1e-6) < 1e-6:
            break
        prev_score = cur_score

    return cur_score


def smoothing_by_diffusion(predicted_df, log2_data, top_col_inds, smooth_batch_size=1000, smooth_cores_to_use=1, seed = 14,
                           maxiter=1e4,rescale=True, rescale_deg=1, rescale_ratio=None):
    # Set seed for reproducibility
    np.random.seed(seed)

    if smooth_batch_size > len(log2_data):
        print("The passed subsample size is greater than the number of cells in the subsample. \n    Now setting subsample size to "+str(len(log2_data))+". ")
        chunk_number = 1
    else:
        chunk_number = math.ceil(len(log2_data) / smooth_batch_size)

    original_names = predicted_df.index
    subsamples_indices = np.arange(len(log2_data))
    np.random.shuffle(subsamples_indices)
    subsamples = np.array_split(subsamples_indices, chunk_number)

    smoothed_scores = []
    smooth_results = []
    # Process each chunk separately
    with concurrent.futures.ProcessPoolExecutor(max_workers=smooth_cores_to_use) as executor:
        for subsample in subsamples:
            chunk_log2_data = log2_data[subsample, :]
            chunk_predicted_df = predicted_df.iloc[subsample, :]
            smooth_results.append(executor.submit(smooth_subset,chunk_log2_data,chunk_predicted_df,top_col_inds,maxiter))
        for f in concurrent.futures.as_completed(smooth_results):
            cur_score = f.result()
            smoothed_scores.append(cur_score)

    # Concatenate the smoothed scores for all chunks
    smoothed_scores_concatenated = pd.concat(smoothed_scores)
   
    return smoothed_scores_concatenated[original_names]


def binning(predicted_df, scores): # scores is smoothed scores
    labels = ['Differentiated',
      'Unipotent',
      'Oligopotent',
      'Multipotent',
      'Pluripotent',
      'Totipotent']
    
    #print('    Started binning')
    pred_potencies = predicted_df["preKNN_CytoTRACE2_Potency"]
    unique_potency = np.unique(pred_potencies)
    score = 'preKNN_CytoTRACE2_Score'
    df_pred_potency = pd.DataFrame({'preKNN_CytoTRACE2_Potency':pred_potencies,'preKNN_CytoTRACE2_Score':scores})
    limits = np.arange(7)/6
    for potency_i, potency in enumerate(labels):
        lower = limits[potency_i]
        upper = limits[potency_i+1]
        if potency in unique_potency:
            data_order =  df_pred_potency[df_pred_potency['preKNN_CytoTRACE2_Potency']==potency]['preKNN_CytoTRACE2_Score'].sort_values()
            index = data_order.index
            n = len(index)
            scaler = MinMaxScaler(feature_range=(lower+1e-8, upper-1e-8))
            order = scaler.fit_transform(np.arange(n).reshape(-1,1))[:,0]
            df_pred_potency.loc[index,score] = order 

    predicted_df["preKNN_CytoTRACE2_Score"] = df_pred_potency[score][predicted_df.index]

    return predicted_df


# Adaptive nearest neighbors smoothing

def map_score_to_potency(score):
    labels = ['Differentiated','Unipotent','Oligopotent','Multipotent','Pluripotent','Totipotent']
    ranges = np.linspace(0, 1, 7)  
    if score <= ranges[1]:
        return labels[0]
    elif score <= ranges[2]:
        return labels[1]
    elif score <= ranges[3]:
        return labels[2]
    elif score <= ranges[4]:
        return labels[3]   
    elif score <= ranges[5]:
        return labels[4]
    elif score <= ranges[6]:
        return labels[5]
    else:
        return np.nan
        

def shortest_consensus(neighbor_scores):
    idx_use = 2
    last_part = False
    for i in range(2,math.floor(len(neighbor_scores)/2+1)):
        if map_score_to_potency(np.mean(neighbor_scores[:i]))==map_score_to_potency(np.mean(neighbor_scores[i:2*i])) and not last_part:
            idx_use = i
            last_part = True
    return 2*idx_use


def neighborhood_smoothing_single_chunk(df_in, df_pca, chunk_cell_names):
    cell_names = df_pca.index        
    new_scores = []
    for cell in chunk_cell_names:
        cell_dist = pairwise_distances(df_pca.loc[cell,:].values.reshape(1, -1),df_pca)[0]
        cell_dist = cell_dist / np.max(cell_dist)
        neighbor_cells = np.argsort(cell_dist)[:30]
        neighbor_dists = cell_dist[neighbor_cells]
        neighbor_scores = df_in.loc[cell_names[neighbor_cells], 'preKNN_CytoTRACE2_Score'].values
        num_neighbors_keep = shortest_consensus(neighbor_scores)
        if num_neighbors_keep > 1:
            new_neighbor_cells = neighbor_cells[:num_neighbors_keep]
            new_neighbor_dists = neighbor_dists[:num_neighbors_keep]
            new_neighbor_scores = neighbor_scores[:num_neighbors_keep]
            neighbor_score_weights = ((1 - new_neighbor_dists) ** 2) / (((1 - new_neighbor_dists) ** 2).sum())
            proposed_new_score = (new_neighbor_scores * (1 - new_neighbor_dists) ** 2).sum() / (
                        (1 - new_neighbor_dists) ** 2).sum()
            new_scores.append(proposed_new_score)
        else:
             new_scores.append(-1)
    return pd.DataFrame(new_scores,columns=['Score'],index=chunk_cell_names)


def neighborhood_smoothing(df_in, log2_data, cores_to_use = 10):
    labels = ['Differentiated','Unipotent','Oligopotent','Multipotent','Pluripotent','Totipotent']
    ranges = np.linspace(0, 1, 7)  
    if len(df_in) < 100:
        print('Fewer than 100 cells, neighborhood smoothing disabled')
        df_in['CytoTRACE2_Score'] = df_in['preKNN_CytoTRACE2_Score']
        df_in['CytoTRACE2_Potency'] = df_in['preKNN_CytoTRACE2_Potency']
        
    data_scale = scale(log2_data,axis=1)
    df_out = df_in.copy()
    df_out['CytoTRACE2_Score'] = 0.0
    df_out['CytoTRACE2_Potency'] = ''

    df_out['Smoothed_Score'] = df_in['preKNN_CytoTRACE2_Score'].copy()
    
    num_pcs = min(30, data_scale.shape[0] - 1)
    pca = PCA(n_components=num_pcs, svd_solver='arpack')
    cell_names = df_in.index
    df_pca = pd.DataFrame(pca.fit_transform(data_scale),index=cell_names,columns=['PC'+str(i+1) for i in range(num_pcs)])
    num_cells = df_in.shape[0]
    num_chunks = cores_to_use
    chunk_size = math.ceil(num_cells/num_chunks)
    results = []
    all_scores_out = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores_to_use) as executor:

        for chunk in range(num_chunks):
            chunk_cell_names = cell_names[(chunk_size*chunk):np.min([chunk_size*(chunk+1),num_cells])]
            results.append(executor.submit(neighborhood_smoothing_single_chunk, df_in, df_pca, chunk_cell_names))
        for f in concurrent.futures.as_completed(results):
            chunk_scores = f.result()
            all_scores_out.append(chunk_scores)

    df_temp = pd.concat(all_scores_out,ignore_index=False)
    df_out.loc[df_temp.index,'CytoTRACE2_Score'] = df_temp['Score'].values
    df_out.loc[df_out['CytoTRACE2_Score']<-0.1,'CytoTRACE2_Score'] = df_out.loc[df_out['CytoTRACE2_Score']<-0.1,'CytoTRACE2_Score']
    df_out['CytoTRACE2_Potency'] = ''
    for i in range(len(labels)):
        range_min = ranges[i]
        range_max = ranges[i + 1]
        df_out.loc[(range_min < df_out['CytoTRACE2_Score']) * (
                    df_out['CytoTRACE2_Score'] <= range_max), 'CytoTRACE2_Potency'] = labels[i]
    return df_out


# String formatting for plots
def format_string(s,max_line_length=25):
    all_chunks = s.split(' ')
    total_num_chunks = len(all_chunks)
    new_string = ''
    current_line_length = 0
    i = 0
    while i < total_num_chunks:  
        chunk = all_chunks[i]
        if len(chunk)>max_line_length:
            line = chunk[:max_line_length]+'...\n'
            new_string += line
            i += 1
        else:
            current_line_length = len(chunk)
            line = chunk
            while (current_line_length < max_line_length) and ((i+1) < total_num_chunks):
                candidate_line_length = len(line+' '+all_chunks[i+1])
                if candidate_line_length < max_line_length:
                    line += ' '+all_chunks[i+1]
                    current_line_length = len(line)
                    i += 1
                else:
                    current_line_length = max_line_length
            new_string += line+'\n'
            i += 1
            
    if new_string.endswith('\n'):
        new_string = new_string[:-1]
    
    return new_string


# Data processing for scanpy UMAPs
def process_with_scanpy(df_exp,df_anno=None):
    adata = sc.AnnData(df_exp)
    if df_anno is not None:
        pheno_col = df_anno.columns[0]
        adata.obs['phenotype'] = df_anno.loc[df_exp.index,pheno_col].copy()
        adata.obs['phenotype_txtwrap'] = [format_string(p) for p in adata.obs['phenotype']]
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    sc.pp.scale(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    del adata.X, adata.obsp, adata.var
    return adata
    


