import concurrent.futures
import math
import numpy as np
import os
import pandas as pd
import scanpy as sc
import subprocess
import torch
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from sklearn.model_selection import train_test_split

from cytotrace2_py.common.gen_utils import *
from cytotrace2_py.common.argument_parser import *
from cytotrace2_py.common.plot import *

def process_subset(idx, chunked_expression, B_in, smooth_batch_size, smooth_cores_to_use, species, use_model_dir, output_dir, seed, disable_verbose):

    # map and rank
    cell_names, gene_names, rank_data, log2_data = preprocess(chunked_expression, species)

    # Check gene counts as QC measure
    gene_counts = (log2_data>0).sum(1)
    low_gc_frac = (gene_counts < 500).mean()
    if low_gc_frac >= 0.2:
        warn_str = "{:.1f}".format(100*low_gc_frac)+"% of input cells express fewer than 500 genes. For best results, a minimum gene count of 500-1000 is recommended. \n    Please see FAQ for guidelines at https://github.com/digitalcytometry/cytotrace2#frequently-asked-questions"
        warnings.warn(warn_str,stacklevel=2)
    # top variable genes
    top_col_inds = top_var_genes(log2_data)
    top_col_names = gene_names[top_col_inds]
    
    # predict by unrandomized chunked batches
    if not disable_verbose:
        print('cytotrace2: Performing initial model prediction')
    predicted_df = predict(rank_data, log2_data, B_in, cell_names, use_model_dir, chunked_expression.shape[0])
    predicted_df['Raw_Score'] = predicted_df['preKNN_CytoTRACE2_Score'].copy()
    if not disable_verbose:
        print('cytotrace2: Performing smoothing by diffusion')
    smooth_score = smoothing_by_diffusion(predicted_df, log2_data, top_col_inds, smooth_batch_size, smooth_cores_to_use, seed) 
    binned_score_pred_df = binning(predicted_df, smooth_score)

    if chunked_expression.shape[0] < 100:
        print('cytotrace2: Fewer than 100 cells in dataset. Skipping KNN smoothing step.')
        binned_score_pred_df['CytoTRACE2_Score'] = binned_score_pred_df['preKNN_CytoTRACE2_Score'].copy()
        return binned_score_pred_df
    else:
        if not disable_verbose:
            print('cytotrace2: Performing smoothing by adaptive KNN')
        smooth_by_knn_df = neighborhood_smoothing(binned_score_pred_df, log2_data, smooth_cores_to_use)
        return smooth_by_knn_df

def calculate_cores_to_use(chunk_number,smooth_chunk_number,max_cores,disable_parallelization):

    pred_cores_to_use = 1
    smooth_cores_to_use = 1
    if smooth_chunk_number == 1:
        print("cytotrace2: The number of cells in your dataset is less than the specified smoothing batch size.\n")
        print("    Model prediction will not be parallelized.")

    if not disable_parallelization:
        # Calculate number of available processors
        num_proc = os.cpu_count()
        #print("cytotrace2: "+str(num_proc)+" cores detected")
        if num_proc == 1:
            print("cytotrace2: Only one core detected. CytoTRACE 2 will not be run in parallel.")
        elif max_cores == None:
            pred_cores_to_use = max(1, num_proc // 2)
            smooth_cores_to_use = min(smooth_chunk_number,max(1, num_proc // 2))
            print('cytotrace2: Running '+str(chunk_number)+' prediction batch(es) sequentially using '+str(smooth_cores_to_use)+' cores per batch.')
        else:
            max_cores = min(max_cores,max(1, num_proc // 2))
            pred_cores_to_use = min(chunk_number,max_cores)
            smooth_cores_to_use = min(smooth_chunk_number,max_cores)
            print('cytotrace2: Running '+str(chunk_number)+' prediction batch(es) sequentially using '+str(smooth_cores_to_use)+' cores per batch.')

    return pred_cores_to_use, smooth_cores_to_use


def cytotrace2(input_path, 
               annotation_path = "",
               species = "mouse",
               batch_size = 20000,
               smooth_batch_size = 1000,
               disable_parallelization = False,
               max_cores = None,
               seed = 14,
               output_dir = 'cytotrace2_results',
               disable_plotting = False,
               disable_verbose = False):

    # Make output directory 
    out = os.system('mkdir -p '+output_dir)

    # Load data
    print('cytotrace2: Input parameters')
    print('    Input file: '+input_path)
    print('    Species: '+species)
    print('    Parallelization enabled: '+str(not disable_parallelization))
    print('    Batch size: '+str(batch_size))
    print('    Smoothing batch size: '+str(smooth_batch_size))
    print('    Seed: '+str(seed))
    print('    Output directory: '+output_dir)
    print('    Plotting enabled: '+str(not disable_plotting))
    print('    Verbose mode enabled: '+str(not disable_verbose))

    print('    User-provided limit for number of cores to use: '+str(max_cores))
    if max_cores is None:
        cpus_detected = os.cpu_count()
        print('       ...'+str(cpus_detected)+' cores detected. CytoTRACE 2 will run using up to '+str(max(cpus_detected // 2, 1))+'/'+str(cpus_detected)+' cores.')

    expression = read_file(input_path)
    print('cytotrace2: Dataset characteristics')
    print('    Number of input genes: ',str(expression.shape[1]))
    print('    Number of input cells: ',str(expression.shape[0]))

    # Raise warning if the data seems to be log2-transformed already
    if expression.max().max() <= 20:
        warnings.warn("Input expression data seem to be log2-transformed. Please provide data as raw counts or CPM/TPM.",stacklevel=1)
    
    if not disable_plotting:
        if not disable_verbose:
            print('cytotrace2: Computing UMAP embeddings from full expression')
        if len(annotation_path)>0:
            df_anno = pd.read_csv(annotation_path,index_col=0,sep='\t')
            adata_plot = process_with_scanpy(expression, df_anno)
        else:
            adata_plot = process_with_scanpy(expression)

    # Check if the input species is accurate
    # Calculate the proportion of row names that are all uppercase (assumed to be human) or not all uppercase (assumed to be mouse)
    is_human = sum([name.isupper() for name in expression.columns]) / expression.shape[1] > 0.9
    is_mouse = sum([not name.isupper() for name in expression.columns]) / expression.shape[1] > 0.9

    if is_human and species == 'mouse':
        warnings.warn("Species is most likely human. Please revise the 'species' input to the function.",stacklevel=1)

    if is_mouse and species == 'human':
        warnings.warn("Species is most likely mouse. Please revise the 'species' input to the function.",stacklevel=1)

    np.random.seed(seed)
    if batch_size > len(expression):
        print("cytotrace2: The passed batch_size is greater than the number of cells in the subsample. \n    Now setting batch_size to "+str(len(expression))+".")
        batch_size <- len(expression)
    elif len(expression) > 50000 and batch_size > 50000:
        print("cytotrace2: Please consider reducing the batch_size to 50000 for runtime and memory efficiency.")
    
    print('cytotrace2: Preprocessing')
    
    # Calculate chunk number
    chunk_number = math.ceil(len(expression) / batch_size)
    smooth_chunk_number = math.ceil(batch_size/smooth_batch_size)
    if len(expression) < 1000:
        chunk_number = 1
        smooth_chunk_number = 1

    # Determine multiprocessing parameters
    pred_cores_to_use, smooth_cores_to_use = calculate_cores_to_use(chunk_number, smooth_chunk_number, max_cores, disable_parallelization)
    torch.set_num_threads(pred_cores_to_use)

    use_model_dir = pkg_resources.resource_filename("cytotrace2_py","resources/models/")
    background_path = pkg_resources.resource_filename("cytotrace2_py","resources/background.pt")
    B = torch.load(background_path)
    B = B.to_dense().T
    original_names = expression.index
    subsamples_indices = np.arange(len(expression)) 
    if chunk_number > 1:
        np.random.shuffle(subsamples_indices)
    subsamples = np.array_split(subsamples_indices, chunk_number)
    
    predictions = []
   
    # Process each chunk separately
    results = []
    for idx in range(chunk_number):
        chunked_expression = expression.iloc[subsamples[idx], :]
        print('cytotrace2: Initiated processing batch '+str(idx+1)+'/'+str(chunk_number)+' with '+str(chunked_expression.shape[0])+' cells')
        smooth_by_knn_df = process_subset(idx, chunked_expression, B, smooth_batch_size, smooth_cores_to_use, species, use_model_dir, output_dir, seed, disable_verbose)
        predictions.append(smooth_by_knn_df)
    
    predicted_df_final = pd.concat(predictions, ignore_index=False)
    predicted_df_final = predicted_df_final.loc[original_names]
    ranges = np.linspace(0, 1, 7)  
    labels = [
        'Differentiated',
        'Unipotent',
        'Oligopotent',
        'Multipotent',
        'Pluripotent',
        'Totipotent']
    
    predicted_df_final['CytoTRACE2_Potency'] = pd.cut(predicted_df_final['CytoTRACE2_Score'], bins=ranges, labels=labels, include_lowest=True)

    all_scores = predicted_df_final['CytoTRACE2_Score'].values
    predicted_df_final['CytoTRACE2_Relative'] = (all_scores-min(all_scores))/(max(all_scores)-min(all_scores))
    predicted_df_final = predicted_df_final[["CytoTRACE2_Score", "CytoTRACE2_Potency" , "CytoTRACE2_Relative", "preKNN_CytoTRACE2_Score", "preKNN_CytoTRACE2_Potency"]]
    predicted_df_final.to_csv(output_dir+'/cytotrace2_results.txt',sep='\t')

    if not disable_plotting:
        print('cytotrace2: Plotting outputs')
        adata_plot.obs = pd.concat([adata_plot.obs,predicted_df_final.loc[adata_plot.obs_names,:]],axis=1)
        #adata_plot.write_h5ad(dir_out+'saved_metadata_for_plots.h5ad')
        plot_all_outputs(adata_plot,output_dir)
    else:
        print('cytotrace2: Plotting disabled')

    print('cytotrace2: Finished.')
    return predicted_df_final

def run_cytotrace2():
    arguments = argument_parser()
    cytotrace2(**arguments)

if __name__ == "__main__":
    arguments = argument_parser()
    cytotrace2(**arguments)
