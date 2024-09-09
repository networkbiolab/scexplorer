import mygene
import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import rpy2.robjects as robjects
import re
import scanpy as sc
import scipy.sparse as sp

from fastapi import HTTPException
from dotenv import load_dotenv
import plotly.io as pio
from plotly.offline import plot
from plotly.subplots import make_subplots
from scipy.sparse import issparse
from scanpy import preprocessing
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail
import matplotlib.colors as mcolors


colors = ["lightgray", "red"]  # Define the transition colors from light gray to red
n_bins = [3]  # Number of bins for each segment
cmap_name = 'custom_lightgray_to_red'

# Create the colormap
cm = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors, N=100)

def create_plotly_colorscale(colormap, N=100):
    # Generate a list of colors from the colormap
    scale = np.linspace(0, 1, N)
    colors = colormap(scale)

    # Convert to Plotly color scale format
    plotly_scale = [(s, mcolors.rgb2hex(c)) for s, c in zip(scale, colors)]
    return plotly_scale

plotly_cm = create_plotly_colorscale(cm)

load_dotenv()

sendgrid_api_key = os.getenv("SENDGRID_API_KEY")


def get_dataset_description(adata):
    num_observations = adata.n_obs
    num_variables = adata.n_vars
    mean_n_genes_by_counts = adata.obs['n_genes_by_counts'].mean()
    mean_total_counts = adata.obs['total_counts'].mean()
    mean_pct_counts_mito = adata.obs['pct_counts_mito'].mean()
    std_n_genes_by_counts = adata.obs['n_genes_by_counts'].std()
    std_total_counts = adata.obs['total_counts'].std()
    std_pct_counts_mito = adata.obs['pct_counts_mito'].std()   
    cell_observations = len(adata.obs.columns.to_list())
    gene_observations = len(adata.var.columns.to_list())
    dataset_description = {
    'num_observations': num_observations,
    'num_variables': num_variables,
    'mean_n_genes_by_counts': f"{mean_n_genes_by_counts:0.2f}",
    'mean_total_counts': f"{mean_total_counts:0.2f}",
    'mean_pct_counts_mito': f"{mean_pct_counts_mito:0.2f}",
    'std_n_genes_by_counts': f"{std_n_genes_by_counts:0.2f}",
    'std_total_counts': f"{std_total_counts:0.2f}",
    'std_pct_counts_mito': f"{std_pct_counts_mito:0.2f}",
    'cell_observations' : cell_observations,
    'gene_observations': gene_observations
    }
    return dataset_description

def get_quality_control_plots(file_path,save_path,species,flavor,root_path,unique_id):
    scrna = sc.read_h5ad(file_path)
    scrna.obs.to_csv(f"{root_path}/{unique_id}/{unique_id}_obs.csv")
    scrna.var.to_csv(f"{root_path}/{unique_id}/{unique_id}_var.csv")
    scanpy_adata = mygene_converter(scrna, species,flavor)
    # print(f"qc_check_names_my_gene = {scanpy_adata.var.index}")
    if species == "human":
        if flavor == "symbol":
            scanpy_adata.var["mito"] = scanpy_adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(scanpy_adata, qc_vars=["mito"], inplace=True)
            plot_paths = plot_violin(scanpy_adata,save_path)
            dataset_description = get_dataset_description(scanpy_adata)
            scanpy_adata.write_h5ad(file_path)

            return plot_paths,dataset_description
        else:
            # adata.var["mito"] = adata.var_names.str.startswith('MT-')
            # sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
            # plot_paths = plot_violin(adata,save_path)
            # dataset_description = get_dataset_description(adata)
            # return plot_paths,dataset_description
            pass
            #revisar el ensembl id de genes mitocondriales
    else:
        if flavor == "symbol":
            scanpy_adata.var["mito"] = scanpy_adata.var_names.str.startswith('mt-')
            sc.pp.calculate_qc_metrics(scanpy_adata, qc_vars=["mito"], inplace=True)
            plot_paths = plot_violin(scanpy_adata,save_path)
            dataset_description = get_dataset_description(scanpy_adata)
            scanpy_adata.write_h5ad(file_path)

            return plot_paths,dataset_description
        else:
            # adata.var["mito"] = adata.var_names.str.startswith('MT-')
            # sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)
            # plot_paths = plot_violin(adata,save_path)
            # dataset_description = get_dataset_description(adata)
            # return plot_paths,dataset_description
            pass
            #revisar el ensembl id de genes mitocondriales   

def load_and_preprocess_data(file_path: str,save_path:str, doublet_detection: bool,root_path: str, unique_id:str, min_genes: int = 200, min_cells: int = 3, mito_threshold: int = 5):
    """
    Load and preprocess the data.

    Parameters:
    - file_path: Path to the .h5ad file.
    - min_genes: Minimum number of genes expressed required for a cell to pass filtering.
    - min_cells: Minimum number of cells for a gene to be kept.
    - mito_threshold: Threshold for mitochondrial gene content. Cells with a percentage above this threshold will be filtered out.

    Returns:
    - info_dic: Dictionary with preprocessing information.
    - adata: Processed AnnData object.
    """
    adata = sc.read_h5ad(file_path)
    n_cells_pre = adata.n_obs
    n_genes_pre = adata.n_vars
    
    # Filter cells and genes based on min_genes and min_cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Identify Mitochondrial Genes
    if len(adata.var_names[adata.var_names.str.startswith('mt-')]) > 0:
        adata.var["mito"] = adata.var_names.str.startswith('mt-')
    elif len(adata.var_names[adata.var_names.str.startswith('MT-')]) > 0:
        adata.var["mito"] = adata.var_names.str.startswith('MT-')    
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)    # Filter cells based on mitochondrial content using the provided threshold
    adata = adata[adata.obs["pct_counts_mito"] < mito_threshold, :]
    
    if doublet_detection:
        sc.external.pp.scrublet(adata)
        adata = adata[adata.obs["predicted_doublet"] == False]
        adata.layers["counts"] = adata.X.copy()
        adata.raw = adata
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.write_h5ad(file_path)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        adata.obs.to_csv(f"{root_path}/{unique_id}/{unique_id}_obs.csv")
        adata.var.to_csv(f"{root_path}/{unique_id}/{unique_id}_var.csv")
        all_paths = plot_violin(adata,save_path)
        plot_total_counts_n_genes_by_count_path = plot_scatter(adata,save_path,"total_counts_n_genes_by_count","total_counts",["n_genes_by_counts"], x_label="Total Counts",y_label="N° Genes by Counts",color="orange")
        plot_total_counts_mitochondrial_counts_path = plot_scatter(adata,save_path,"total_counts_mitochondrial_counts","total_counts",["pct_counts_mito"], x_label="Total Counts",y_label="% Mitochondrial Counts",color="green")
        plot_hvg_path = plot_scatter_hvg(adata,save_path,"hvg","means",["dispersions_norm"],x_label="Normalized Dispersion",y_label="Mean Expression")
        all_paths.extend([plot_total_counts_n_genes_by_count_path,plot_total_counts_mitochondrial_counts_path,plot_hvg_path])
        n_cells_post_mito_filter = adata.n_obs 
        info_dic = {
            "cells_number_pre_processing": n_cells_pre,
            "genes_number_pre_processing": n_genes_pre,
            "cells_number_post_processing": n_cells_post_mito_filter,
            "genes_number_post_processing": adata.n_vars,
            "processing_plot_paths": all_paths,
            "message": "Data processed successfully."
        }
        return info_dic
    else:
        adata.layers["counts"] = adata.X.copy()
        adata.raw = adata
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.write_h5ad(file_path)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        adata.obs.to_csv(f"{root_path}/{unique_id}/{unique_id}_obs.csv")
        adata.var.to_csv(f"{root_path}/{unique_id}/{unique_id}_var.csv")
        all_paths = plot_violin(adata,save_path)
        plot_total_counts_n_genes_by_count_path = plot_scatter(adata,save_path,"total_counts_n_genes_by_count","total_counts",["n_genes_by_counts"], x_label="Total Counts",y_label="N° Genes by Counts",color="orange")
        plot_total_counts_mitochondrial_counts_path = plot_scatter(adata,save_path,"total_counts_mitochondrial_counts","total_counts",["pct_counts_mito"], x_label="Total Counts",y_label="% Mitochondrial Counts",color="green")
        plot_hvg_path = plot_scatter_hvg(adata,save_path,"hvg","means",["dispersions_norm"],x_label="Normalized Dispersion",y_label="Mean Expression")
        all_paths.extend([plot_total_counts_n_genes_by_count_path,plot_total_counts_mitochondrial_counts_path,plot_hvg_path])
        n_cells_post_mito_filter = adata.n_obs 
        info_dic = {
            "cells_number_pre_processing": n_cells_pre,
            "genes_number_pre_processing": n_genes_pre,
            "cells_number_post_processing": n_cells_post_mito_filter,
            "genes_number_post_processing": adata.n_vars,
            "processing_plot_paths": all_paths,
            "message": "Data processed successfully."
        }
        return info_dic
    
def run_pca_on_data(file_path,save_path,top_genes,hvg_flavor,root_path,unique_id):
    adata = sc.read_h5ad(file_path)
    adata.var["mito"] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True)    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=int(top_genes),flavor=hvg_flavor)
    adata_hvg= adata.copy()
    adata_hvg= adata_hvg[:, adata_hvg.var.highly_variable]
    sc.tl.pca(adata_hvg, svd_solver='arpack')
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.uns["pca"] = adata_hvg.uns["pca"]
    adata.write(file_path)
    adata_hvg.write_h5ad(f"{save_path}/embedding.h5ad")
    adata.obs.to_csv(f"{root_path}/{unique_id}/{unique_id}_obs.csv")
    adata.var.to_csv(f"{root_path}/{unique_id}/{unique_id}_var.csv")
    plot_pca_percentage_mito_path = plot_pca_and_obs(pca_coords=adata.obsm["X_pca"][:, :2],fig_name="pca_percentage_mito",color_column="pct_counts_mito",title=f"% of Mitochondrial Genes",data_df=adata.obs,save_path=save_path,output_format="html")
    plot_pca_total_counts_path  = plot_pca_and_obs(pca_coords=adata.obsm["X_pca"][:, :2],fig_name="pca_total_counts",color_column="total_counts",title="Total Counts",data_df=adata.obs,save_path=save_path,output_format="html")
    plot_pca_n_genes_per_count_path = plot_pca_and_obs(pca_coords=adata.obsm["X_pca"][:, :2],fig_name="pca_n_genes_by_counts",color_column="n_genes_by_counts",title="N° Genes per Counts",data_df=adata.obs,save_path=save_path,output_format="html")
    elbow_plot = plot_elbow(adata=adata,save_path=save_path,fig_name="pca_elbow_plot")
    
    return [plot_pca_percentage_mito_path,plot_pca_total_counts_path,plot_pca_n_genes_per_count_path,elbow_plot]

def run_clustering_on_data(file_path,save_path,n_neighbors,n_pcs,resolution,root_path,unique_id):
    adata_all_genes = sc.read_h5ad(file_path)
    adata = sc.read_h5ad(f"{save_path}/embedding.h5ad")
    # adata.raw = adata
    # adata = adata[:, adata.var.highly_variable]
    # # sc.pp.scale(adata, max_value=10)#ADD AS PARAMETER;FOR THE MOMENT DEFAULT 
    # sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, int(n_neighbors), int(n_pcs))
    sc.tl.umap(adata)
    sc.tl.leiden(adata,resolution=resolution)
    adata_all_genes.obs['leiden'] = adata.obs['leiden']
    adata_all_genes.obsm['X_umap'] = adata.obsm['X_umap'] 
    adata_all_genes.uns["umap"] = adata.uns["umap"]
    adata_all_genes.obs.to_csv(f"{root_path}/{unique_id}/{unique_id}_obs.csv")
    adata_all_genes.var.to_csv(f"{root_path}/{unique_id}/{unique_id}_var.csv")
    adata_all_genes.write_h5ad(file_path)
    sc.tl.dendrogram(adata,groupby="leiden")
    sc.tl.rank_genes_groups(adata,groupby="leiden", method='wilcoxon')
    # clustering_info = len(adata.obs['leiden'].unique()) #usar para dar info
    plot_umap_leiden_path = plot_umap_and_obs(umap_coords=adata.obsm["X_umap"],fig_name="umap_leiden",color_column="leiden",title="Leiden Clusters",data_df=adata.obs,save_path=save_path)
    plot_umap_percentage_mito_path = plot_umap_and_obs(umap_coords=adata.obsm["X_umap"],fig_name="umap_percentage_mito",color_column="pct_counts_mito",title=f"% of Mitochondrial Genes",data_df=adata.obs,save_path=save_path)
    plot_umap_total_counts_path = plot_umap_and_obs(umap_coords=adata.obsm["X_umap"],fig_name="umap_total_counts",color_column="total_counts",title="Total Counts",data_df=adata.obs,save_path=save_path)
    plot_umap_n_genes_per_count_path = plot_umap_and_obs(umap_coords=adata.obsm["X_umap"],fig_name="umap_n_genes_by_counts",color_column="n_genes_by_counts",title="N° Genes per Counts",data_df=adata.obs,save_path=save_path)
    plot_embedding_paths = [plot_umap_leiden_path,plot_umap_percentage_mito_path,plot_umap_total_counts_path,plot_umap_n_genes_per_count_path]
    return plot_embedding_paths

def genes_visualization(file_path,save_path,dim_red,gene_list):
    adata = sc.read_h5ad(file_path)
    genes_vis_plots_path = []
    for gene in gene_list:
        gene_expression = adata[:, gene].X.toarray()  
        adata.obs[gene] = gene_expression.ravel() 
        if dim_red =="X_umap":
            template_plot = plot_umap_and_obs(umap_coords=adata.obsm[dim_red],fig_name=f"{dim_red}_{gene}",color_column=gene,title=f"{gene}",data_df=adata.obs,save_path=save_path)
            genes_vis_plots_path.append(template_plot)
        else:
            template_plot = plot_pca_and_obs(pca_coords=adata.obsm[dim_red],fig_name=f"{dim_red}_{gene}",color_column=gene,title=f"{gene}",data_df=adata.obs,save_path=save_path,output_format="svg")
            genes_vis_plots_path.append(template_plot)       
    return genes_vis_plots_path

def perform_batch_correction(file_path,save_path, batch_correction_method, unique_id, preprocess):
    adata_concat = sc.read_h5ad(file_path)
    if preprocess:
        sc.pp.normalize_per_cell(adata_concat, counts_per_cell_after=1e4)
        sc.pp.log1p(adata_concat)
        sc.pp.highly_variable_genes(adata_concat, n_top_genes=2000)
        # sc.pp.regress_out(adata_concat, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata_concat, max_value=10)
        sc.tl.pca(adata_concat, svd_solver='arpack')

    if batch_correction_method == 'combat':
        sc.pp.combat(adata_concat, key='batch')
        sc.pp.neighbors(adata_concat,use_rep="X_pca")
        sc.tl.umap(adata_concat)
        sc.tl.leiden(adata_concat,resolution=0.5)
        sc.pl.umap(adata_concat,color=["leiden"])
        plot_umap_batch_path = plot_umap_and_obs(umap_coords=adata_concat.obsm["X_umap"],fig_name="combat_integration_umap_batch",color_column="batch",title="Combat Integration",data_df=adata_concat.obs,save_path=save_path)
        plot_umap_leiden_path = plot_umap_and_obs(umap_coords=adata_concat.obsm["X_umap"],fig_name="combat_integration_umap_leiden",color_column="leiden",title="Leiden Clusters",data_df=adata_concat.obs,save_path=save_path)
        integration_plots = [plot_umap_batch_path,plot_umap_leiden_path]
        
    elif batch_correction_method == 'scanorama':
        sc.external.pp.scanorama_integrate(adata_concat, key='batch')
        sc.pp.neighbors(adata_concat,use_rep="X_scanorama")
        sc.tl.umap(adata_concat)
        sc.tl.leiden(adata_concat,resolution=0.5)
        sc.pl.umap(adata_concat,color=["leiden"])
        plot_umap_batch_path = plot_umap_and_obs(umap_coords=adata_concat.obsm["X_umap"],fig_name="scanorama_integration_umap_batch",color_column="batch",title="Scanorama Integration",data_df=adata_concat.obs,save_path=save_path)
        plot_umap_leiden_path = plot_umap_and_obs(umap_coords=adata_concat.obsm["X_umap"],fig_name="scanorama_integration_umap_leiden",color_column="leiden",title="Combat Integration",data_df=adata_concat.obs,save_path=save_path)
        integration_plots = [plot_umap_batch_path,plot_umap_leiden_path]
        
    elif batch_correction_method == 'bbknn':
        sc.external.pp.bbknn(adata_concat, batch_key='batch')
        sc.pp.neighbors(adata_concat,use_rep="X_pca")
        sc.tl.umap(adata_concat)
        sc.tl.leiden(adata_concat,resolution=0.5)
        sc.pl.umap(adata_concat,color=["leiden"])
        plot_umap_batch_path = plot_umap_and_obs(umap_coords=adata_concat.obsm["X_umap"],fig_name="bbknn_integration_umap_batch",color_column="batch",title="BBKNN Integration",data_df=adata_concat.obs,save_path=save_path)
        plot_umap_leiden_path = plot_umap_and_obs(umap_coords=adata_concat.obsm["X_umap"],fig_name="bbknn_integration_umap_leiden",color_column="leiden",title="Combat Integration",data_df=adata_concat.obs,save_path=save_path)
        integration_plots = [plot_umap_batch_path,plot_umap_leiden_path]

    elif batch_correction_method == 'harmony':
        sc.external.pp.harmony_integrate(adata_concat, key='batch')
        sc.pp.neighbors(adata_concat,use_rep="X_pca_harmony")
        sc.tl.umap(adata_concat)
        sc.tl.leiden(adata_concat,resolution=0.5)
        plot_umap_batch_path = plot_umap_and_obs(umap_coords=adata_concat.obsm["X_umap"],fig_name="harmony_integration_umap_batch",color_column="batch",data_df=adata_concat.obs,save_path=save_path)
        plot_umap_leiden_path = plot_umap_and_obs(umap_coords=adata_concat.obsm["X_umap"],fig_name="harmony_integration_umap_leiden",color_column="leiden",data_df=adata_concat.obs,save_path=save_path)
        integration_plots = [plot_umap_batch_path,plot_umap_leiden_path]
    else:
        raise HTTPException(status_code=400, detail="Invalid batch correction method")
    
    adata_concat.write_h5ad(os.path.join(save_path, f"{unique_id}.h5ad"))
    
    return integration_plots
    
def dea(file_path,save_path,n_genes,flavor,gene_list,method,root_path,unique_id):
    adata = sc.read_h5ad(file_path)
    sc.tl.dendrogram(adata,groupby="leiden")
    sc.tl.rank_genes_groups(adata,groupby="leiden", method=method)
    adata.obs.to_csv(f"{root_path}/{unique_id}/{unique_id}_obs.csv")
    adata.var.to_csv(f"{root_path}/{unique_id}/{unique_id}_var.csv")
    save_rank_genes_groups(uuid=unique_id,adata=adata)
    adata = adata.raw.to_adata()
    if isinstance(gene_list, list) and len(gene_list) > 0:        
        if flavor == "mean_expression":##FIX LATER
            plot_mean_expression_list = plotly_dotplot_unified(adata,save_path,"mean_expression_list",n_genes,gene_list=gene_list,group_by="leiden",flavor="mean_expression")
            plot_fold_change_list = plotly_dotplot_unified(adata,save_path,"fold_change_list",n_genes,gene_list=gene_list,group_by="leiden",flavor="fold_change")
            plot_mean_expression = plotly_dotplot_unified(adata,save_path,"mean_expression",n_genes,group_by="leiden",flavor="mean_expression")
            plot_fold_change_expression = plotly_dotplot_unified(adata,save_path,"fold_change",n_genes,group_by="leiden",flavor="fold_change")
            dea_plots = [plot_mean_expression_list,plot_fold_change_list,plot_mean_expression,plot_fold_change_expression]
            return dea_plots
        else:
            plot_fold_change_expression_list = plotly_dotplot_unified(adata,save_path,"fold_change_list",n_genes,gene_list=gene_list,group_by="leiden",flavor="fold_change")
            plot_fold_change_list = plotly_dotplot_unified(adata,save_path,"fold_change_list",n_genes,gene_list=gene_list,group_by="leiden",flavor="fold_change")
            plot_mean_expression = plotly_dotplot_unified(adata,save_path,"mean_expression",n_genes,group_by="leiden",flavor="mean_expression")
            plot_fold_change_expression = plotly_dotplot_unified(adata,save_path,"fold_change",n_genes,group_by="leiden",flavor="fold_change")
            dea_plots = [plot_fold_change_expression_list,plot_fold_change_list,plot_mean_expression_list,plot_mean_expression,plot_fold_change_expression]
            return dea_plots 
    else:
        plot_mean_expression = plotly_dotplot_unified(adata,save_path,"mean_expression",n_genes,group_by="leiden",flavor="mean_expression")
        plot_fold_change_expression = plotly_dotplot_unified(adata,save_path,"fold_change",n_genes,group_by="leiden",flavor="fold_change")
        dea_plots = [plot_mean_expression,plot_fold_change_expression]
        return dea_plots
    
def var_names_adata(file_path):
    adata = sc.read_h5ad(file_path)
    gene_list = adata.var.index.tolist()
    return gene_list

def fetch_analysis_results(adata):
    n_clusters = len(adata.obs['leiden'].unique())
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    return {
        "clusters": n_clusters,
        "cells_number": n_cells,
        "genes_number":n_genes
    }

def plot_elbow(adata,save_path,fig_name):
    variance_ratio = adata.uns["pca"]["variance_ratio"].tolist()
    log_variance_ratio = np.log(variance_ratio[:25])
    percentage_change = np.diff(log_variance_ratio) / log_variance_ratio[:-1] * 100
    pca_labels = ['PCA_{}'.format(i+1) for i in range(len(log_variance_ratio))]
    colors = ['DarkGreen'] * len(log_variance_ratio)
    found = False
    for i, change in enumerate(percentage_change, start=1):
        if abs(change) < 1 and not found:
            colors[i] = 'red'
            found = True
            break  

    trace = go.Scatter(
        x=list(range(1, len(log_variance_ratio) + 1)),
        y=log_variance_ratio,
        mode='markers',
        marker=dict(
            color=colors,  
            size=15
        ),
        hoverinfo='text',  
        hovertext=pca_labels  
    )

    text_annotations = [
        go.layout.Annotation(
            x=xi,
            y=yi,
            text=text,
            showarrow=False,
            xanchor='center',
            yanchor='bottom',
            textangle=-90,
            font=dict(
                size=9 
            ),
            yshift=10  
        )
        for xi, yi, text in zip(range(1, len(log_variance_ratio) + 1), log_variance_ratio, pca_labels)
    ]

    layout = go.Layout(
        xaxis=dict(
            title='Ranking',
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=False
        ),
        yaxis=dict(
            title='Log of Variance Ratio',
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=False
        ),
        annotations=text_annotations,
        paper_bgcolor='white', 
        plot_bgcolor='white',
        width=600,  # Set the width of the plot
        height=400     
    )

    fig = go.Figure(data=[trace], layout=layout)
    file_path = plot(fig, filename=f"{save_path}/{fig_name}.html", auto_open=False)
    return file_path

def plot_density_contour(adata, 
                         template="plotly_white", 
                         x_hist_color='lightblue', 
                         y_hist_color='lightblue', 
                         contour_line_color='deepskyblue'):
    """
    Generate a density contour plot using Plotly.
    
    Parameters:
    - adata: the data to plot
    - template: string for the plotly template
    - x_hist_color: color for the x-axis histogram
    - y_hist_color: color for the y-axis histogram
    - contour_line_color: color for the contour lines
    
    Returns:
    - Displays the plot
    """
    
    fig = px.density_contour(adata.obs, 
                             x="log1p_total_counts", 
                             y="log1p_n_genes_by_counts", 
                             marginal_x="histogram", 
                             marginal_y="histogram", 
                             template=template)

    # Change histogram colors
    fig.data[-2]['marker']['color'] = x_hist_color  # For x-axis histogram
    fig.data[-1]['marker']['color'] = y_hist_color  # For y-axis histogram

    # Update contour colors
    for data in fig.data:
        if data['type'] == 'contour':
            data['line']['color'] = contour_line_color
            data['line']['width'] = 0.5

    # Adjust the overall figure width
    fig.update_layout(
        xaxis_title='log(1 + Total counts)',
        yaxis_title='log(1 + N° Genes by counts)',
        width=600  # Adjust the value as needed
    )

    fig.show()

def plot_violin(adata, save_path, template="plotly_white", features=['n_genes_by_counts', 'total_counts', 'pct_counts_mito'], y_labels=["N° genes by counts", "Total counts", f"% of Mitochondrial Genes"]):
    """
    Plot violin plots for the specified features using Plotly and save them as HTML files.
    
    Parameters:
    - adata: Data object, usually from anndata library.
    - save_path: Path where to save the HTML plots.
    - template: Plotly template for styling the plot.
    - features: List of features (columns) in adata.obs to plot.
    - y_labels: List of y-axis labels corresponding to each feature.
    
    Note:
    - Ensure that the length of 'features' and 'y_labels' are the same.
    """
    paths_plots = []
    # Pre-defined color sequence for the violin plots
    colors = [
        "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", 
        "#bcbd22", "#17becf", "#1f77b4", "#ff7f0e", "#2ca02c"
    ]
    output_format = "png" if len(adata.obs) > 15000 else "html"

    # Iterate through each feature and its corresponding y-axis label
    for idx, (feature, y_label) in enumerate(zip(features, y_labels)):
        # Create violin plot for the current feature
        fig = px.violin(adata.obs, y=feature, box=True, points="all", template=template, color_discrete_sequence=[colors[idx]])
        
        # Adjust the size of the dots inside the violin plot
        for trace in fig.data:
            trace.marker.size = 2
        
        # Update layout and settings for the plot
        fig.update_layout(
            width=250,
            height=400,
            yaxis_title=y_label,
            showlegend=False,
            margin=dict(l=80, r=20, b=40, t=20)  # Setting margins for the plot
        )
        # Save the plot as an HTML file
        file_path = f"{save_path}/{feature}_violin.{output_format}"
        if output_format == "html":
            plot(fig, filename=file_path, auto_open=False)
        else:
            fig.write_image(file_path, scale=10)  # High definition PNG
        paths_plots.append(file_path)

    return paths_plots

def plot_scatter(adata,save_path,fig_name, x, y, template="plotly_white",x_label=None,y_label=None,color="black"):
    """
    Plot scatter plots based on provided x and y columns from adata.obs.
    
    Parameters:
    - adata: Data object, usually from anndata library.
    - x: Column name in adata.obs for x-axis.
    - y: List of column names in adata.obs for y-axis.
    - template: Plotly template for styling the plot.
    """
    output_format = "png" if len(adata.obs) > 15000 else "html"

    for idx, y_col in enumerate(y):
        fig = go.Figure(data=go.Scatter(x=adata.obs[x], 
                                        y=adata.obs[y_col],
                                        mode='markers',
                                        marker=dict(
                                            color=color,
                                            opacity=0.5  
                                        )))
        
        fig.update_layout(
            xaxis_title=x if x_label is None else x_label,
            yaxis_title=y_col if y_label is None else y_label,
            width=280,
            template=template,
            showlegend=False,
            margin=dict(l=60, r=60, b=40, t=40)  
        )
        file_path = f"{save_path}/{fig_name}.{output_format}"
        if output_format == "html":
            plot(fig, filename=file_path, auto_open=False)
        else:
            fig.write_image(file_path, scale=10)  # High definition PNG
    return file_path        

def highest_expr_genes_plotly(
    adata,
    n_top=30,
    gene_symbols=None,
    log=False,
    template="white",
    **kwargs
):
    """
    Display the highest expressed genes in a dataset using a box plot.

    Parameters:
    ----------
    - adata: An AnnData object containing expression data.
    - n_top: Number of top genes to display.
    - gene_symbols: Column name in adata.var that stores gene symbols.
    - log: Boolean indicating whether y-axis should be log-scaled.
    - template: Style template for the plot. Either "white" or "dark".
    - **kwargs: Additional keyword arguments passed to the box plot.

    Returns:
    --------
    A Plotly figure object.
    """

    # Define a color palette for the box plots
    colors = [
        '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', 
        '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', 
        '#ccebc5', '#ffed6f', '#1f78b4', '#33a02c', '#e31a1c'
    ]

    # Normalize the expression data to compute the percentage of each gene per cell
    norm_dict = preprocessing.normalize_total(adata, target_sum=100, inplace=False)

    # Determine the genes with the highest mean expression
    if issparse(norm_dict['X']):
        mean_percent = norm_dict['X'].mean(axis=0).A1
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict['X'][:, top_idx].A
    else:
        mean_percent = norm_dict['X'].mean(axis=0)
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict['X'][:, top_idx]
    
    # Extract the names of the top genes, using gene symbols if provided
    columns = (
        adata.var_names[top_idx]
        if gene_symbols is None
        else adata.var[gene_symbols][top_idx]
    )
    counts_top_genes = pd.DataFrame(
        counts_top_genes, index=adata.obs_names, columns=columns
    )

    # Initialize a Plotly figure
    fig = go.Figure()

    # Add each gene as a separate trace (box plot) to the figure
    for idx, col in enumerate(counts_top_genes.columns):
        fig.add_trace(go.Box(x=counts_top_genes[col], name=col, boxpoints='outliers', marker_color=colors[idx % len(colors)], **kwargs))

    # Set the template based on the input
    plot_template = "plotly_white" if template == "white" else "plotly_dark"

    # Set the layout and styling details for the figure
    fig.update_layout(
        xaxis_title="% of total counts",
        xaxis=dict(categoryorder='total descending'),
        height=600,
        width=550,
        paper_bgcolor='white' if template == "white" else 'black',
        plot_bgcolor='white' if template == "white" else 'black',
        template=plot_template,
        showlegend=False,
        margin=dict(l=60, r=60, b=40, t=40)  # Padding: left, right, bottom, top
    )

    # If log scale is requested for the y-axis, apply it
    if log:
        fig.update_layout(yaxis_type="log")

    return fig


        
def prepare_data(adata, cell_types, n_genes, gene_list, group_by, flavor, min_logfc=1):
    """
    Prepare gene expression data for plotting. 

    This function processes the AnnData object to extract the mean expression 
    and fraction of cells expressing specific genes, either provided through 
    a gene list or by identifying top expressed genes.

    Parameters:
    ----------
    adata : AnnData
        An object containing gene expression data and potentially 
        differential expression results.
        
    cell_types : list
        List of cell types to consider.
        
    n_genes : int
        Number of top genes to consider if gene_list is not provided.
        
    gene_list : list, optional
        List of specific genes to consider. If provided, n_genes is ignored.
        
    group_by : str
        Column name in adata.obs that contains cell type information.
        
    flavor : str
        Specifies the method for obtaining gene values. Should be one of 
        'mean_expression' or 'fold_change'. If 'mean_expression', the function 
        will retrieve the mean expression of each gene. If 'fold_change', 
        the function will fetch the fold change value from differential 
        expression results stored in adata.uns.
        
    min_logfc : float, optional
        Minimum log fold change to consider when selecting genes based 
        on fold change. Only used if flavor is 'fold_change'.

    Returns:
    -------
    tuple
        A tuple containing three arrays:
        - Mean expression data across cell types for selected genes.
        - Fraction of cells expressing each gene across cell types.
        - List of all genes considered.
    """
    data = []
    fraction_expressing_data = []
    all_genes = []
    
    if gene_list is not None:  # If gene_list is provided, collect data once for each gene in gene_list
        selected_genes = gene_list
        all_genes.extend(selected_genes)
        
        for gene in selected_genes:
            values = []
            fraction_values = []
            
            for ct in cell_types:
                subset = adata[adata.obs[group_by] == ct]
                
                if flavor == "mean_expression":
                    value = subset[:, gene].X.mean()
                elif flavor == "fold_change":
                    value = adata.uns["rank_genes_groups"]["logfoldchanges"][ct][
                        np.where(adata.uns["rank_genes_groups"]["names"][ct] == gene)
                    ][0]
                else:
                    raise ValueError("Invalid flavor. Choose 'mean_expression' or 'fold_change'.")
                
                fraction_value = np.sum(subset[:, gene].X > 0) / subset.shape[0]
                
                values.append(value)
                fraction_values.append(fraction_value)
            
            data.append(values)
            fraction_expressing_data.append(fraction_values)
    
    else:  # Original logic when gene_list is not provided
        for cell_type in cell_types:
            top_genes = adata.uns["rank_genes_groups"]["names"][cell_type][:n_genes * 2].tolist()  # Taking 2x to have a buffer
            selected_genes = []
            
            for gene in top_genes:
                if len(selected_genes) >= n_genes:
                    break  # Stop if we already have enough genes
                
                if flavor == "fold_change":
                    logfc_value = adata.uns["rank_genes_groups"]["logfoldchanges"][cell_type][
                        np.where(adata.uns["rank_genes_groups"]["names"][cell_type] == gene)
                    ][0]
                    
                    if logfc_value >= min_logfc:
                        selected_genes.append(gene)
                else:
                    selected_genes.append(gene)
            
            # Check if we have enough genes; if not, fill with top genes irrespective of fold change
            remaining_genes_needed = n_genes - len(selected_genes)
            for gene in top_genes:
                if gene not in selected_genes:
                    selected_genes.append(gene)
                    remaining_genes_needed -= 1
                
                if remaining_genes_needed <= 0:
                    break
            
            all_genes.extend(selected_genes)
            
            for gene in selected_genes:
                values = []
                fraction_values = []
                
                for ct in cell_types:
                    subset = adata[adata.obs[group_by] == ct]
                    
                    if flavor == "mean_expression":
                        value = subset[:, gene].X.mean()
                    elif flavor == "fold_change":
                        value = adata.uns["rank_genes_groups"]["logfoldchanges"][ct][
                            np.where(adata.uns["rank_genes_groups"]["names"][ct] == gene)
                        ][0]
                    else:
                        raise ValueError("Invalid flavor. Choose 'mean_expression' or 'fold_change'.")
                    
                    fraction_value = np.sum(subset[:, gene].X > 0) / subset.shape[0]
                    
                    values.append(value)
                    fraction_values.append(fraction_value)
                
                data.append(values)
                fraction_expressing_data.append(fraction_values)
    
    return np.array(data), np.array(fraction_expressing_data), all_genes

def add_annotations(fig, cell_types, font_family, n_genes, gene_list=None):
    """
    Annotate the Plotly figure with cell type labels and horizontal bars.

    This function adds text annotations to the figure, representing the cell types. 
    If a gene list is not provided, it also adds horizontal bars with cell type 
    labels above the dot plot, providing a visual summary of the number of top genes 
    for each cell type.

    Parameters:
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which the annotations will be added.

    cell_types : list
        List of cell types to annotate.

    font_family : str
        Font family to use for the text annotations.

    n_genes : int
        Number of top genes considered. Used to determine the length of the horizontal bars.

    gene_list : list, optional
        List of specific genes considered. If provided, horizontal bars are not added.

    Returns:
    -------
    None
    """
    # Initialize variables for bar plotting
    bar_y_start = len(cell_types) + 0.1  # The y-position where the bars will start
    bar_x_start = 0  # The x-position where the first bar will start
    for idx, cell_type in enumerate(cell_types):
        fig.add_annotation(
            x=-1.25,
            y=idx,
            xref="x",
            yref="y",
            text=cell_type,
            showarrow=False,
            font=dict(
                size=12,
                color="black",
                family=font_family
            ),
            row=1,
            col=1
        )

    if gene_list is None:

        # Add the horizontal bars and vertical lines on top of the scatter plot
        for j, cell_type in enumerate(cell_types):
            bar_x_end = bar_x_start + n_genes - 1  # The x-position where the bar will end
            mid_x = (bar_x_start + bar_x_end) / 2  # Midpoint for text positioning

            # Add the horizontal bar
            fig.add_trace(
                go.Scatter(
                    x=[bar_x_start, bar_x_end],
                    y=[bar_y_start, bar_y_start],
                    mode='lines',
                    line=dict(color="black", width=1),
                    hoverinfo='none'
                ),
                row=1, col=1
            )

            # Add the cell type text at the midpoint of the bar
            fig.add_trace(
                go.Scatter(
                    x=[mid_x],
                    y=[bar_y_start],
                    mode='text',
                    text=[cell_type],
                    textposition="top center",
                    hoverinfo='none'
                ),
                row=1, col=1
            )

            # Add vertical lines (square brackets)
            for x_pos in [bar_x_start, bar_x_end]:
                fig.add_trace(
                    go.Scatter(
                        x=[x_pos, x_pos],
                        y=[bar_y_start - 0.2, bar_y_start + 0.02],
                        mode='lines',
                        line=dict(color="black", width=1),
                        hoverinfo='none'
                    ),
                    row=1, col=1
                )

            bar_x_start = bar_x_end + 1

    fig.add_annotation(
        x=-2.5,
        y=len(cell_types) / 2,
        text="Cell Types",
        textangle=-90,
        showarrow=False,
        font=dict(
            size=14,
            color="black",
            family=font_family
        ),
        row=1,
        col=1
    )
    fig.update_xaxes(fixedrange=True, domain=[0, 0.99], row=1, col=1)
    fig.update_yaxes(fixedrange=True, range=[-1, bar_y_start + 1], row=1, col=1)
    fig.update_xaxes(fixedrange=True, domain=[0.95, 0.985], row=1, col=2)
    fig.update_yaxes(fixedrange=True, range=[-5, 105], row=1, col=2)

def set_plot_layout(fig, all_genes, cell_types, font_family, base_height, base_width, scaling_factor):
    """
    Configure the layout and style of a Plotly dot plot figure.

    This function sets the appearance and layout of the figure, such as axis titles, 
    tick values, figure dimensions, and font styling. 

    Parameters:
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure whose layout will be configured.

    all_genes : list
        List of all genes considered, used to set the x-axis tick values.

    cell_types : list
        List of cell types considered, used for layout scaling calculations.

    font_family : str
        Font family to be used for all text elements in the figure.

    base_height : int
        Base height of the figure.

    base_width : int
        Base width of the figure before scaling.

    scaling_factor : float
        Factor by which the figure width should be scaled, based on the number of genes.

    Returns:
    -------
    dict
        A dictionary of configuration options for interactive features of the figure.
    """
    adjusted_width = int(base_width * scaling_factor)
    config = {
        'toImageButtonOptions': {
            'format': 'svg',
            'filename': 'custom_image',
            'height': base_height,
            'width': adjusted_width,
            'scale': 1
        }
    }
    fig.update_layout(
        xaxis1=dict(tickvals=list(range(len(all_genes))), ticktext=all_genes, title="Genes", showgrid=False, zeroline=False, tickangle=270),
        yaxis=dict(showgrid=False,zeroline=False,showline=False,ticks='',showticklabels=False),
        yaxis2=dict(showticklabels=False, showgrid=False, zeroline=False),
        xaxis2=dict(showticklabels=False, showgrid=False, zeroline=False),
        showlegend=False,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=1, r=30, t=20, b=50),
        font=dict(family=font_family, color="black"),
        width=adjusted_width,
        height=base_height
    )
    return config

def add_scatter_plot(fig, mean_expression_data, fraction_expressing_data, all_genes, cell_types, colorscale, max_size, flavor):
    """
    Add a scatter plot representation of gene expression data to a Plotly figure.

    This function generates and adds to the figure a scatter plot showing gene 
    expression across cell types. The size of each point represents the fraction 
    of cells expressing the gene, while the color indicates the mean expression 
    level or fold change, depending on the flavor.

    Parameters:
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which the scatter plot will be added.

    mean_expression_data : numpy.ndarray
        2D array containing mean expression data of genes across cell types.

    fraction_expressing_data : numpy.ndarray
        2D array containing the fraction of cells expressing each gene across cell types.

    all_genes : list
        List of all genes considered, used for x-axis positioning.

    cell_types : list
        List of cell types, used for y-axis positioning.

    colorscale : str
        Colorscale for the scatter plot points, representing expression level or fold change.

    max_size : int
        Maximum size for the scatter plot points.

    flavor : str
        Specifies the method for obtaining gene values: 'mean_expression' or 'fold_change'.
        Determines the data visualized by the scatter plot's color.

    Returns:
    -------
    None
    """
    normalized_sizes = (fraction_expressing_data / np.max(fraction_expressing_data)) * max_size
    x_data, y_data, sizes, colors, hover_texts = [], [], [], [], []
    
    for i, gene in enumerate(all_genes):
        for j, cell_type in enumerate(cell_types):
            x_data.append(i)
            y_data.append(j)
            sizes.append(normalized_sizes[i, j])
            colors.append(mean_expression_data[i, j])
            hover_texts.append(f"Gene: {gene}<br>Fraction: {fraction_expressing_data[i, j]*100:.2f}%<br>Value: {mean_expression_data[i, j]:.2f}")
    
    final_colorscale = colorscale if flavor == "mean_expression" else "balance"
    color_range = [0, np.max(mean_expression_data)] if flavor == "mean_expression" else [-4, 4]
    
    fig.add_trace(
        go.Scatter(
            x=x_data,
            y=y_data,
            mode='markers',
            marker=dict(
                size=sizes,
                color=colors,
                colorscale=final_colorscale,
                showscale=False,
                colorbar=dict(title="", orientation='h'),
                cmin=color_range[0],
                cmax=color_range[1],
                line=dict(width=1, color="Black"),
                sizemode="diameter",
                opacity=1
            ),
            hoverinfo="text",
            hovertext=hover_texts
        ),
        row=1, col=1
    )

def add_dendrogram_to_plot(fig, adata, cell_types, color_mapping=None, shift_y=0):
    """
    Add a dendrogram to an existing Plotly figure.
    
    This function extracts dendrogram information from the AnnData object 
    and plots it onto the provided figure. The dendrogram typically represents
    the hierarchical clustering structure of the data.

    Parameters:
    ----------
    fig: plotly.graph_objects.Figure
        The Plotly figure to which the dendrogram will be added.
        
    adata: AnnData
        An AnnData object, which is expected to contain dendrogram data 
        under the 'dendrogram_leiden' key in its uns attribute.
        
    color_mapping: dict, optional
        A dictionary mapping the original colors in the dendrogram to desired 
        colors for visualization. If a color is not provided in this mapping,
        it defaults to black.

    Returns:
    -------
    None
    """
    # Extract dendrogram information from the AnnData object
    icoord = adata.uns["dendrogram_leiden"]['dendrogram_info']['icoord']
    dcoord = adata.uns["dendrogram_leiden"]['dendrogram_info']['dcoord']
    color_list = adata.uns["dendrogram_leiden"]['dendrogram_info']['color_list']

    y_min = min(min(xs) for xs in icoord)
    y_max = max(max(xs) for xs in icoord)
    scale_factor = 1.1

    # Plot dendrogram lines onto the figure
    for xs, ys, color in zip(icoord, dcoord, color_list):
        scaled_xs = [(x - y_min) * scale_factor for x in xs]  # Scale the y-coordinates
        shifted_xs = [x + shift_y for x in scaled_xs]  # Shift the y-coordinates
        line_color = color_mapping.get(color, 'black')
        fig.add_trace(
            go.Scatter(
                x=ys,
                y=shifted_xs,
                mode='lines',
                line=dict(color=line_color, width=1),
                hoverinfo='none'
            ),
            row=1, col=2
        )
        
def plotly_dotplot_unified(adata,save_path,fig_name, n_genes=10, gene_list=None, group_by="leiden", 
                           max_size=14, show_dendrogram=True, dendrogram_color='black', 
                           font_family='Arial', colorscale="Reds", flavor="mean_expression", min_logfc=1):
    """
    Create a unified dot plot visualization of gene expression using Plotly.

    This function produces a dot plot showcasing gene expression across different cell types.
    Each dot's size represents the fraction of cells expressing the gene, while its color 
    indicates the mean expression level or fold change. Additional annotations and an optional 
    dendrogram can also be added to provide context about the hierarchical structure of the data.

    Parameters:
    ----------
    adata : AnnData
        An object containing gene expression data, differential expression results, and 
        potentially dendrogram data.

    n_genes : int, optional
        Number of top genes to consider if gene_list is not provided. Default is 10.

    gene_list : list, optional
        List of specific genes to consider. If provided, n_genes is ignored.

    group_by : str, optional
        Column name in adata.obs that contains cell type or cluster information. Default is 'leiden'.

    max_size : int, optional
        Maximum size for the scatter plot points. Default is 15.

    show_dendrogram : bool, optional
        Whether to overlay a dendrogram on the plot, indicating hierarchical clustering. Default is True.

    dendrogram_color : str, optional
        Color for the dendrogram lines. Default is 'black'.

    font_family : str, optional
        Font family to use for all text elements in the figure. Default is 'Arial'.

    colorscale : str, optional
        Colorscale for the scatter plot points, representing expression level or fold change. Default is "Reds".

    flavor : str, optional
        Specifies the method for obtaining gene values: 'mean_expression' or 'fold_change'. Default is 'mean_expression'.

    min_logfc : float, optional
        Minimum log fold change to consider when selecting genes based on fold change. Only used if flavor is 'fold_change'. Default is 1.

    Returns:
    -------
    plotly.graph_objects.Figure
        The fully constructed Plotly figure.
    """
    
    cell_types = list(adata.uns["dendrogram_leiden"]['categories_ordered'])
    
    if gene_list is not None:
        n_genes = len(gene_list)  
    
    mean_expression_data, fraction_expressing_data, all_genes = prepare_data(
        adata, cell_types, n_genes, gene_list, group_by, flavor, min_logfc)   # Pass min_logfc here
    
    fig = make_subplots(rows=1, cols=2, shared_yaxes=False, column_widths=[0.98, 0.02])
    
    add_scatter_plot(fig, mean_expression_data, fraction_expressing_data, all_genes, cell_types, colorscale, max_size, flavor)
    
    base_height = 350
    
    if gene_list is not None:
        if len(gene_list) < 20:
            base_width = 75
            scaling_factor = 5
        else:
            base_width = 55
            scaling_factor = 15
    else:
        scaling_factor = max(1, n_genes / 2)
        base_width = 500
    config = set_plot_layout(fig, all_genes, cell_types, font_family, base_height, base_width, scaling_factor)
    
    add_annotations(fig, cell_types, font_family, n_genes, gene_list=gene_list)
    
    if show_dendrogram:
        add_dendrogram_to_plot(fig, adata, cell_types, color_mapping={dendrogram_color: 'black'}, shift_y=4)  # Adjust shift_y as needed
    
    output_format = "png" if len(adata.obs) > 15000 else "html"

    file_path = f"{save_path}/{fig_name}.{output_format}"
    if output_format == "html":
        plot(fig, filename=file_path, auto_open=False, config=config)
    else:
        fig.write_image(file_path, scale=10) 
    
    return file_path


def plot_umap_and_obs(umap_coords, save_path, fig_name, data_df,title, color_column=None, template="plotly_white", output_format="html"):
    """
    Plot UMAP scatter plot using Plotly based on colors from a specified column.

    Parameters:
    - umap_coords: numpy array containing UMAP coordinates.
    - data_df: DataFrame containing metadata and columns for coloring.
    - color_column: Column name in data_df to color the scatter plot.
    - template: Plotly template for styling the plot.
    """
    
    # Convert UMAP coordinates to DataFrame
    umap_df = pd.DataFrame(umap_coords, columns=['UMAP1', 'UMAP2'])
    
    # Concatenate UMAP DataFrame with data_df side by side
    merged_df = pd.concat([umap_df, data_df.reset_index(drop=True)], axis=1)
    
    # Sort the DataFrame based on the condition
    if color_column:
        merged_df = merged_df.sort_values(by=color_column)
    
    num_unique_classes = merged_df[color_column].nunique() if color_column else 0

    # Create the scatter plot
    fig = px.scatter(
        merged_df,
        x='UMAP1',
        y='UMAP2',
        color=color_column,
        color_continuous_scale=plotly_cm,
        template=template
    )

    # Get the label for the color column and set as title
    if title is None:
        title = data_df[color_column].name if color_column else "UMAP Scatter Plot"

    show_legend = True if num_unique_classes < 30 else False

    # Update the layout and settings for the plot
    fig.update_layout(
        title=title,
        title_x=0.4,
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            linecolor='gray',
            linewidth=1.5,
            mirror=True,
            showline=True,
            showspikes=False,  # Hide axis spikelines
            title_text=""  # Empty title to hide it
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            linecolor='gray',
            linewidth=1.5,
            mirror=True,
            showline=True,
            showspikes=False,  # Hide axis spikelines
            title_text=""  # Empty title to hide it
        ),
        showlegend=show_legend,
        margin=dict(l=40, r=40, b=60, t=60),  # Adjusted bottom margin
        plot_bgcolor='white',
        width=400,
        height=380,
        annotations=[
            # X-axis label
            dict(
                xref="paper",
                yref="paper",
                x=0.5,
                y=-0.1,  # Adjusted y position
                showarrow=False,
                text="UMAP1",
                font=dict(size=14)
            ),
            # Y-axis label
            dict(
                xref="paper",
                yref="paper",
                x=-0.08,
                y=0.5,
                showarrow=False,
                text="UMAP2",
                textangle=-90,
                font=dict(size=14)
            )
        ]
    )
    
    # Remove colorbar title
    fig.update_layout(coloraxis_colorbar_title_text="")
    
    # Display the plot
    fig.update_traces(marker=dict(size=4))
    output_format = "png" if len(data_df) > 15000 else "html"

    file_path = f"{save_path}/{fig_name}.{output_format}"
    if output_format == "html":
        pio.write_html(fig, file=file_path, auto_open=False)
    elif output_format == "svg":
        fig.write_image(file_path)
    else:
        fig.write_image(file_path, scale=10)  

    return file_path

def plot_pca_and_obs(pca_coords, save_path, fig_name, data_df,title, color_column=None, template="plotly_white", output_format="html"):
    """
    Plot PCA scatter plot using Plotly based on colors from a specified column.

    Parameters:
    - pca_coords: numpy array containing PCA coordinates.
    - data_df: DataFrame containing metadata and columns for coloring.
    - color_column: Column name in data_df to color the scatter plot.
    - template: Plotly template for styling the plot.
    """
    
    # Convert PCA coordinates to DataFrame
    pca_df = pd.DataFrame(pca_coords[:, :2], columns=['PC1', 'PC2'])
    
    # Concatenate PCA DataFrame with data_df side by side
    merged_df = pd.concat([pca_df, data_df.reset_index(drop=True)], axis=1)
    
    # Sort the DataFrame based on the condition
    if color_column:
        merged_df = merged_df.sort_values(by=color_column)
    
    num_unique_classes = merged_df[color_column].nunique() if color_column else 0

    # Create the scatter plot
    fig = px.scatter(
        merged_df,
        x='PC1',
        y='PC2',
        color=color_column,
        color_continuous_scale=plotly_cm,
        template=template
    )

    # # Get the label for the color column and set as title
    if title is None:
        title = data_df[color_column].name if color_column else "PCA Scatter Plot"

    show_legend = True if num_unique_classes < 30 else False

    # Update the layout and settings for the plot
    fig.update_layout(
        title=title,
        title_x=0.4,
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            linecolor='gray',
            linewidth=1.5,
            mirror=True,
            showline=True,
            showspikes=False,  # Hide axis spikelines
            title_text=""  # Empty title to hide it
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            linecolor='gray',
            linewidth=1.5,
            mirror=True,
            showline=True,
            showspikes=False,  # Hide axis spikelines
            title_text=""  # Empty title to hide it
        ),
        showlegend=show_legend,
        margin=dict(l=40, r=40, b=60, t=60),  # Adjusted bottom margin
        plot_bgcolor='white',
        width=400,
        height=380,
        annotations=[
            # X-axis label
            dict(
                xref="paper",
                yref="paper",
                x=0.5,
                y=-0.1,  # Adjusted y position
                showarrow=False,
                text="PC1",
                font=dict(size=14)
            ),
            # Y-axis label
            dict(
                xref="paper",
                yref="paper",
                x=-0.08,
                y=0.5,
                showarrow=False,
                text="PC2",
                textangle=-90,
                font=dict(size=14)
            )
        ]
    )
    
    # Remove colorbar title
    fig.update_layout(coloraxis_colorbar_title_text="")
    
    # Display the plot
    fig.update_traces(marker=dict(size=4))
    output_format = "png" if len(data_df) > 15000 else "html"

    file_path = f"{save_path}/{fig_name}.{output_format}"
    if output_format == "html":
        pio.write_html(fig, file=file_path, auto_open=False)
    elif output_format == "svg":
        fig.write_image(file_path)
    else:
        fig.write_image(file_path, scale=10)  
        
    return file_path

def plot_scatter_hvg(adata,save_path,fig_name, x, y,  template="plotly_white", x_label=None, y_label=None):
    """
    Plot scatter plots based on provided x and y columns from adata.

    Parameters:
    - adata: Data object, usually from anndata library.
    - x: Column name in adata for x-axis.
    - y: List of column names in adata for y-axis.
    - source: Specifies where to extract data from. Either "obs" or "var".
    - template: Plotly template for styling the plot.
    - x_label: Custom label for x-axis. If not provided, defaults to column name.
    - y_label: Custom label for y-axis. If not provided, defaults to column name.
    """
    output_format = "png" if len(adata.obs) > 15000 else "html"

    # Determine the source based on the provided parameter
    data_source = adata.var
    
    for idx, y_col in enumerate(y):
        # Extract highly variable and other genes
        hv_indices = adata.var["highly_variable"]
        
        # Create scatter plot for Highly Variable Genes
        fig = go.Figure(data=go.Scatter(x=data_source[x][hv_indices], 
                                        y=data_source[y_col][hv_indices],
                                        mode='markers',
                                        marker=dict(opacity=0.5),
                                        name="Highly Variable Genes"))
        
        # Add scatter plot for Other Genes
        fig.add_trace(go.Scatter(x=data_source[x][~hv_indices], 
                                 y=data_source[y_col][~hv_indices],
                                 mode='markers',
                                 marker=dict(opacity=0.5),
                                 name="Other Genes"))
        
        # Adjust layout for legend position
        fig.update_layout(
            width=280,
            template=template,
            xaxis_title=x if x_label is None else x_label,
            yaxis_title=y_col if y_label is None else y_label,
            margin=dict(l=60, r=60, b=40, t=40),
            legend=dict(
                x=1,      # 1 is the far right, so this pushes the legend to the right edge
                y=1,      # 1 is the top, so this pushes the legend to the top edge
                xanchor="right",   # Anchor the legend's right edge at x=1
                yanchor="top"      # Anchor the legend's top at y=1
            )
        )
        file_path = f"{save_path}/{fig_name}.{output_format}"
        if output_format == "html":
            plot(fig, filename=file_path, auto_open=False)
        else:
            fig.write_image(file_path, scale=10)        
    return file_path
        
def send_email(email: str, subject: str, content: str):
    message = Mail(
        from_email='scrnaexplorer@gmail.com',
        to_emails=email,
        subject=subject,
        html_content=content
    )
    try:
        sg = SendGridAPIClient(sendgrid_api_key)
        response = sg.send(message)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
    
def upload_notification(mail, name_analysis, id_for):
    template_path = "mail/templates/upload_notification.html"
    content = read_template(template_path)
    content = content.replace("name_analysis", name_analysis) \
                        .replace("id_for_email", id_for)
    send_email(mail, "Sc Explorer Upload Confirmation 🧬", content)
    
def read_template(template_path):
    with open(template_path, 'r', encoding='utf-8') as file:
        content = file.read()
    return content

def sanitize_column_names(dataframe):
    """Sanitize column names to have only alphanumeric characters and underscores."""
    dataframe.columns = [re.sub(r'\W+', '', col) for col in dataframe.columns]
    return dataframe

# def run_r_script(uuid):

#     r_script = f"""
#     library(Seurat)
#     library(SeuratDisk)

#     # Specify the path to your h5ad file
#     h5ad_path <- 'uploads/{uuid}/seurat.h5ad'
#     h5seurat_path <- 'uploads/{uuid}/seurat.h5seurat' # This is the corrected part

#     # Specify the path where you want to save the RDS file
#     rds_path <- 'uploads/{uuid}/{uuid}.rds'
#     Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE) # Specify the output format and overwrite if exists

#     # Load the converted h5Seurat file as a Seurat object
#     seuratObject <- LoadH5Seurat(h5seurat_path,assays = "RNA") 

#     # Save the Seurat object as an RDS file
#     saveRDS(seuratObject, file = rds_path)
#     """
    
#     input = sc.read_h5ad(f"uploads/{uuid}/{uuid}.h5ad")
#     adata = input.copy()
#     if "counts" in adata.layers:
#         matrix = adata.layers["counts"]
#         adata.obs.columns = [sub.replace('(', '') for sub in adata.obs.columns]
#         adata.obs.columns = [sub.replace(')', '') for sub in adata.obs.columns]
#         adata.obs.columns = [sub.replace('/', '') for sub in adata.obs.columns]
#         adata.obs.columns = [sub.replace('=', '.') for sub in adata.obs.columns]
#         adata.obs.columns = [sub.replace(' ', '_') for sub in adata.obs.columns]
#         adata.obs.columns = [sub.replace('-', '_') for sub in adata.obs.columns]
#         adata.obs_names = [sub.replace('-', '_') for sub in adata.obs_names]
#         new_adata = sc.AnnData(X=matrix, obs=sanitize_column_names(adata.obs))
#         new_adata.write_h5ad(f"uploads/{uuid}/seurat.h5ad")

#         robjects.r(r_script)
#         os.remove(f"uploads/{uuid}/seurat.h5ad")
#         os.remove(f"uploads/{uuid}/seurat.h5seurat")
#         path = f"uploads/{uuid}/{uuid}.rds"
#         return path
#     else:
#         input = sc.read_h5ad(f"uploads/{uuid}/{uuid}.h5ad")
#         input.write_h5ad(f"uploads/{uuid}/seurat.h5ad")
#         robjects.r(r_script)
#         os.remove(f"uploads/{uuid}/seurat.h5ad")
#         os.remove(f"uploads/{uuid}/seurat.h5seurat")
#         path = f"uploads/{uuid}/{uuid}.rds"
#         return path
def run_r_script(uuid):

    r_script = f"""
    library(Seurat)
    library(anndata)
    library(Matrix)

    h5ad_path <- 'uploads/{uuid}/seurat.h5ad'
    rds_path <- 'uploads/{uuid}/{uuid}.rds'
    data <- read_h5ad(h5ad_path)
    colnames(data$X) <- gsub('_', '-', colnames(data$X))
    if (inherits(data$X, 'sparseMatrix')) data$X <- as.matrix(data$X)
    data <- CreateSeuratObject(counts = t(data$X))
    saveRDS(data, file = rds_path)
    """

    clean_data = sc.read_h5ad(f"uploads/{uuid}/{uuid}.h5ad")
    adata = clean_data.copy()
 
    adata.write_h5ad(f"uploads/{uuid}/seurat.h5ad")

    robjects.r(r_script)
    os.remove(f"uploads/{uuid}/seurat.h5ad")
    path = f"uploads/{uuid}/{uuid}.rds"
    return path

    
def clustree_seurat(uuid):
    input = sc.read_h5ad(f"uploads/{uuid}/{uuid}.h5ad")
    adata = input.copy()
    if not sp.issparse(adata.X):
        print("matrix to csr for clustree")
        adata.X = sp.csr_matrix(adata.X)
    if "counts" in adata.layers:
        matrix = adata.layers["counts"]
        adata.obs.columns = [sub.replace('(', '') for sub in adata.obs.columns]
        adata.obs.columns = [sub.replace(')', '') for sub in adata.obs.columns]
        adata.obs.columns = [sub.replace('/', '') for sub in adata.obs.columns]
        adata.obs.columns = [sub.replace('=', '.') for sub in adata.obs.columns]
        adata.obs.columns = [sub.replace(' ', '_') for sub in adata.obs.columns]
        adata.obs.columns = [sub.replace('-', '_') for sub in adata.obs.columns]
        adata.obs_names = [sub.replace('-', '_') for sub in adata.obs_names]
        new_adata = sc.AnnData(X=matrix, obs=sanitize_column_names(adata.obs))
        new_adata.write_h5ad(f"uploads/{uuid}/clustree.h5ad")
    
    r_script = f"""
    
    library(Seurat)
    library(anndata)
    library(Matrix)
    library(clustree)
    set.seed(1234)

    h5ad_path <- 'uploads/{uuid}/clustree.h5ad'
    data <- read_h5ad(h5ad_path)
    colnames(data$X) <- gsub('_', '-', colnames(data$X))
    if (inherits(data$X, 'sparseMatrix')) data$X <- as.matrix(data$X)
    seuratObject <- CreateSeuratObject(counts = t(data$X))

    seuratObject <- NormalizeData(seuratObject)
    seuratObject <- FindVariableFeatures(seuratObject)
    seuratObject <- ScaleData(seuratObject)
    seuratObject <- RunPCA(seuratObject)
    seuratObject <- FindNeighbors(seuratObject)
    seuratObject <- FindClusters(seuratObject, res = seq(0, 2, by = 0.1))
    clustreePlot <- clustree(seuratObject, prefix = "RNA_snn_res.")

    ggsave("uploads/{uuid}/embedding/clustreePlot.svg", plot = clustreePlot, width = 10, height = 10)
    """
    robjects.r(r_script)
    os.remove(f"uploads/{uuid}/clustree.h5ad")

    path = f"uploads/{uuid}/embedding/clustreePlot.svg"
    return path


def mygene_converter(adata, species, flavor="symbol"):
    mg = mygene.MyGeneInfo()
    gene_ids = adata.var.index.tolist()
    first_id = gene_ids[0]

    if (flavor == "symbol" and not first_id.startswith(('ENSG', 'ENSM'))) or \
       (flavor == "ensembl" and first_id.startswith(('ENSG', 'ENSM'))):
        return adata

    if flavor == "symbol":
        scopes = 'ensembl.gene'
        fields = 'symbol'
    elif flavor == "ensembl":
        scopes = 'refseq,symbol'
        fields = 'ensembl.gene'
    else:
        return "wrong flavor"
    
    out = mg.querymany(gene_ids, scopes=scopes, fields=fields, species=species, returnall=True,as_dataframe=True)
    out["out"]["symbol"] = out["out"]["symbol"].fillna(pd.Series(out["out"].index, index=out["out"].index)) 
    df_unique = out["out"][~out["out"].index.duplicated(keep='first')]

    # try:
    #     out = mg.querymany(gene_ids, scopes=scopes, fields=fields, species=species, returnall=True)
    # except Exception as e:
    #     print(f"Error querying MyGeneInfo: {e}")
    #     return "Error Mygene"
    
    # output_dict = {}
    # for item in out['out']:
    #     if 'notfound' in item:
    #         output_dict[item['query']] = item['query']
    #     elif 'ensembl' in item:
    #         ensembl = item['ensembl']
    #         if isinstance(ensembl, list):
    #             output_dict[item['query']] = ensembl[0]['gene']
    #         else:
    #             output_dict[item['query']] = ensembl['gene']

    # final_output = [output_dict.get(tid, tid) for tid in gene_ids]
    adata.var.index = df_unique["symbol"].tolist()
    return adata
    


def save_rank_genes_groups(uuid, adata):
    def recarray_to_dict_with_dtype(rec_array):
        data_dict = {}
        for field in rec_array.dtype.names:
            data_dict[field] = rec_array[field].tolist()
        return data_dict

    def flatten_dict_with_dtype(data_dict):
        flat_list = []
        dtype_list = []
        for field, values in data_dict.items():
            flat_list.extend(values)
            dtype_list.extend([field] * len(values))
        return flat_list, dtype_list

    # Ensure the output directory exists
    output_dir = f'uploads/{uuid}/dea/results'
    os.makedirs(output_dir, exist_ok=True)

    adata_uns_rank_genes_groups = {
        'names': adata.uns["rank_genes_groups"]["names"],
        'scores': adata.uns["rank_genes_groups"]["scores"],
        'pvals': adata.uns["rank_genes_groups"]["pvals"],
        'pvals_adj': adata.uns["rank_genes_groups"]["pvals_adj"],
        'logfoldchanges': adata.uns["rank_genes_groups"]["logfoldchanges"]
    }

    names_dict = recarray_to_dict_with_dtype(adata_uns_rank_genes_groups['names'])
    scores_dict = recarray_to_dict_with_dtype(adata_uns_rank_genes_groups['scores'])
    pvals_dict = recarray_to_dict_with_dtype(adata_uns_rank_genes_groups['pvals'])
    pvals_adj_dict = recarray_to_dict_with_dtype(adata_uns_rank_genes_groups['pvals_adj'])
    logfoldchanges_dict = recarray_to_dict_with_dtype(adata_uns_rank_genes_groups['logfoldchanges'])

    # Flatten the lists and add dtype information
    names_list, names_dtype = flatten_dict_with_dtype(names_dict)
    scores_list, scores_dtype = flatten_dict_with_dtype(scores_dict)
    pvals_list, pvals_dtype = flatten_dict_with_dtype(pvals_dict)
    pvals_adj_list, pvals_adj_dtype = flatten_dict_with_dtype(pvals_adj_dict)
    logfoldchanges_list, logfoldchanges_dtype = flatten_dict_with_dtype(logfoldchanges_dict)

    # Creating a DataFrame with dtype column
    df = pd.DataFrame({
        'cluster': names_dtype,
        'names': names_list,
        'scores': scores_list,
        'pvals': pvals_list,
        'pvals_adj': pvals_adj_list,
        'logfoldchanges': logfoldchanges_list,
    })

    combined_dict = {
        'names': names_dict,
        'scores': scores_dict,
        'pvals': pvals_dict,
        'pvals_adj': pvals_adj_dict,
        'logfoldchanges': logfoldchanges_dict
    }

    dataframes = {}
    for dtype_index in names_dict.keys():
        data = {
            'names': combined_dict['names'][dtype_index],
            'scores': combined_dict['scores'][dtype_index],
            'pvals': combined_dict['pvals'][dtype_index],
            'pvals_adj': combined_dict['pvals_adj'][dtype_index],
            'logfoldchanges': combined_dict['logfoldchanges'][dtype_index]
        }
        dataframes[dtype_index] = pd.DataFrame(data)

    for dtype_index, df in dataframes.items():
        df.to_csv(f'{output_dir}/cluster_{dtype_index}.csv', index=False)

