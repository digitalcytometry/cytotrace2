import anndata as ad
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

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


def plot_all_outputs(adata,dir_out):
    if 'phenotype' in adata.obs.columns:
        ##################################################
        # Plot UMAP colored by phenotype
        ##################################################
        sc.set_figure_params(figsize=(4, 4))
        sc.set_figure_params(scanpy=True, fontsize=14,figsize=(4,4))
        plt.rc('legend',fontsize=14)
        fig = sc.pl.umap(adata,color='phenotype_txtwrap',title='Phenotype',return_fig=True)
        ax = fig.axes[0]
        ax.spines[['right', 'top']].set_visible(False)
        ax.set_aspect(1 / ax.get_data_ratio())
        plt.savefig(dir_out+'/phenotype_UMAP.png',bbox_inches='tight')

        ##################################################
        # Plot boxplot of potency scores by phenotype
        ##################################################
        pheno_median_score = adata.obs.groupby("phenotype",observed=True)['CytoTRACE2_Score'].median().sort_values()
        adata.obs['CytoTRACE2_Score_Phenotype_Median'] = list(pheno_median_score.loc[adata.obs['phenotype'].values].values)
        cytotrace2_pal = [plt.get_cmap('Spectral', 51)(round(max(0, min(5.5 - 6 * cytotrace2, 5)) * 10)) for cytotrace2 in pheno_median_score.values]
        
        violinfont = {'fontname':'Helvetica','fontsize':14}
        xlabel_font = {'fontname':'Helvetica','fontsize':10}
        
        pheno_order = list(adata.obs.groupby("phenotype",observed=True)['CytoTRACE2_Score'].median().sort_values(ascending=False).index)
        plt.figure(figsize=(6+len(pheno_order)/5, 6))
        pheno_order_format = [format_string(s) for s in pheno_order]
        ax = sns.boxplot(data=adata.obs, x='phenotype', y='CytoTRACE2_Score',hue='CytoTRACE2_Score_Phenotype_Median',showfliers=False,
                         order=pheno_order,linewidth=1,palette=cytotrace2_pal,legend=False)
        ax.grid(False)
        sns.stripplot(data=adata.obs, x='phenotype', y='CytoTRACE2_Score', alpha=0.25, dodge=True,
                      jitter=True, size=1, order=pheno_order,legend=False, color='k')
        plt.xticks(pheno_order,pheno_order_format,rotation=90,**xlabel_font)
        plt.yticks(**violinfont)
        plt.ylabel('Potency score',**violinfont)
        plt.xlabel('Phenotype')
        plt.ylim(0,1.0)
        plt.xlim(-0.75,len(pheno_order)-0.25)
        # Plot lines demarcating potency categories and add labels
        for i in range(6):
            plt.hlines((i + 1)/6, -0.75, len(pheno_order)-0.25, linestyles='-', colors='grey', linewidth=0.15)
        
        ax_right = ax.twinx()
        ax_right.grid(False)
        ax_right.set_ylim([0,1])
        
        ax_right.set_yticks(np.linspace(1/12,11/12,6),['Differentiated','Unipotent','Oligopotent','Multipotent','Pluripotent','Totipotent'],**violinfont)
        ax_right.tick_params(axis='y', left=False, right=False, labelright=True, pad=0)
        
        plt.title("Developmental potential by phenotype",**violinfont,pad=10)
        sns.despine()
        plt.tight_layout()
        plt.savefig(dir_out+'/cytotrace2_potency_score_by_phenotype.png')
        
    ##################################################
    # Plot UMAP colored by predicted potency category
    ##################################################
    sc.set_figure_params(scanpy=True, fontsize=14)
    plt.rc('legend',fontsize=14)
    labels = ['Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent', 'Totipotent']
    colors = ["#5E4FA2", "#66C2A5", "#E6F598", "#FEE08B", "#F46D43", "#9E0142"]
    potency_palette = dict(zip(labels,colors))
    adata.obs['CytoTRACE2_Potency'] = adata.obs['CytoTRACE2_Potency'].astype(str)
    fig = sc.pl.umap(adata,color='CytoTRACE2_Potency',palette=potency_palette,title='CytoTRACE 2',return_fig=True)
    ax = fig.axes[0]
    ax.legend_.set_title("Potency category")
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_aspect(1 / ax.get_data_ratio())
    plt.savefig(dir_out+'/cytotrace2_potency_category_UMAP.png',bbox_inches='tight')

    ##################################################
    # Plot UMAP colored by predicted potency score
    ##################################################
    # Rescale values for UMAP color scheme 
    adata.obs['CytoTRACE2_Score_Plot'] = 6*adata.obs['CytoTRACE2_Score']-5.5
    
    colors = ["#5E4FA2", "#66C2A5", "#E6F598", "#FEE08B", "#F46D43", "#9E0142"]
    pad_labels = ['','Differentiated', 'Unipotent','Oligopotent','Multipotent','Pluripotent','Totipotent']
    
    sc.set_figure_params(scanpy=True, fontsize=14)
    plt.rc('legend',fontsize=14)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("potency_cols",list(zip(np.linspace(0,1,6),colors)))
    fig = sc.pl.umap(adata,color='CytoTRACE2_Score_Plot',cmap=cmap,title='CytoTRACE 2',return_fig=True,vmin=-5,vmax=0,colorbar_loc=None)
    ax = plt.gca()
    clb = plt.colorbar(ax.get_children()[0],ax=ax,fraction=0.035);
    clb.ax.set_title('Potency score',pad=10,loc='left')
    clb.ax.set_yticks(list(np.linspace(-5,0,7)))
    clb.ax.set_yticklabels(pad_labels,va='top')
    ax.spines[['right', 'top']].set_visible(False)
    
    # Offset potency category labels to center of regions
    offset = matplotlib.transforms.ScaledTranslation(0,-1/12, fig.dpi_scale_trans)
    for label in clb.ax.yaxis.get_majorticklabels():
        label.set_transform(label.get_transform() + offset)
    ax.set_aspect(1 / ax.get_data_ratio())
    
    plt.savefig(dir_out+'/cytotrace2_potency_score_UMAP.png',bbox_inches='tight')

    ##################################################
    # Plot UMAP colored by predicted potency score
    ##################################################

    rel_colors = ["#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF"]
    rel_labels = ["0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"]
    
    sc.set_figure_params(scanpy=True, fontsize=14)
    plt.rc('legend',fontsize=14)
    rel_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("potency_cols",list(zip(np.linspace(0,1,6),rel_colors)))
    fig = sc.pl.umap(adata,color='CytoTRACE2_Relative',cmap=rel_cmap,title='CytoTRACE 2',return_fig=True,vmin=0,vmax=1,colorbar_loc=None)
    ax = plt.gca()
    clb = plt.colorbar(ax.get_children()[0],ax=ax,fraction=0.035);
    clb.ax.set_title('Relative order',pad=15,loc='left')
    clb.ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
    clb.ax.set_yticklabels(rel_labels,va='center')
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_aspect(1 / ax.get_data_ratio())
    
    plt.savefig(dir_out+'/cytotrace2_relative_order_UMAP.png',bbox_inches='tight')

