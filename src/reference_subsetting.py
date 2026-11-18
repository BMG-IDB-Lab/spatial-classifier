import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import ast

def load_neighbors(file_path):
    neighbors = {}
    with open(file_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            key, values = line.split(":", 1)
            neighbors[key.strip()] = ast.literal_eval(values.strip())
    return neighbors



def create_meta_and_adata(region_of_interest, write_reference=False):
    INPUT_REGION = region_of_interest
  
    dir_list = ['Zhuang-ABCA-1', 'Zhuang-ABCA-2','Zhuang-ABCA-3', 'Zhuang-ABCA-4']
    
    cell = {}
    for d in dir_list:
        cell[d] = pd.read_csv(f'data/{d}/metadata/20241115/cell_metadata.csv')
        
        cell[d].set_index('cell_label', inplace=True)
        
        sdf = cell[d].groupby('brain_section_label')
        
        print(d,":","Number of cells = ", len(cell[d]), ", ", "Number of sections =", len(sdf))
    
    
    cluster_colors= pd.read_csv('data/WMB_taxonomy/cluster_to_cluster_annotation_membership.csv')
    cluster_colors.set_index('cluster_alias', inplace=True)
    cluster_details = pd.read_csv('data/WMB_taxonomy/views/cluster_to_cluster_annotation_membership_pivoted.csv')
    cluster_details.set_index('cluster_alias', inplace=True)
    
    
    cell_extended = {}
    
    for d in dir_list:
        cell_extended[d] = cell[d].join(cluster_details, on='cluster_alias')
        cell_extended[d] = cell_extended[d].join(cluster_colors, on='cluster_alias')
    
    
    ccf_coordinates = {}
    
    for d in dir_list :
    
        ccf_coordinates[d] = pd.read_csv(f'data/{d}/ccf_coordinates/20230830//ccf_coordinates.csv')
        ccf_coordinates[d].set_index('cell_label', inplace=True)
        ccf_coordinates[d].rename(columns={'x': 'x_ccf',
                                           'y': 'y_ccf',
                                           'z': 'z_ccf'},
                                  inplace=True)
        
        cell_extended[d] = cell_extended[d].join(ccf_coordinates[d], how='inner')
    
    parcellation_annotation = pd.read_csv('/raid/zstarinnov/merf_base_data/abc_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv')
    parcellation_annotation.set_index('parcellation_index', inplace=True)
    
    parcellation_color = pd.read_csv('/raid/zstarinnov/merf_base_data/abc_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_color.csv')
    parcellation_color.set_index('parcellation_index', inplace=True)
    
    parcellation_color.columns = ['parcellation_%s'% x for x in  parcellation_color.columns]
    
    for d in dir_list :
        cell_extended[d] = cell_extended[d].join(parcellation_annotation, on='parcellation_index')
        cell_extended[d] = cell_extended[d].join(parcellation_color, on='parcellation_index')
    
    for i in dir_list:
        del cell_extended[f'{i}']['brain_section_label']
        cell_extended[f'{i}'] = cell_extended[f'{i}'].loc[~cell_extended[f'{i}'].index.duplicated(keep='first'), :]
    
    adatas = []
    for i in dir_list:
        adata = sc.read_h5ad(f'data/{i}/expression_matrices/20230830/{i}-raw.h5ad')
        genes = pd.read_csv(f'data/{i}/metadata/20241115/gene.csv')
        genes.set_index('gene_identifier', inplace =True)
        del genes['gene_symbol']
        adata.var = adata.var.join(genes, how='inner')
        print(len(adata))
        
        cell_extended[f'{i}'].index = cell_extended[f'{i}'].index.astype('str')
        adata = adata[adata.obs.index.isin(cell_extended[f'{i}'].index)].copy()
        adata.obs = adata.obs.join(cell_extended[f'{i}'], how='inner')
        print(len(adata))
        adatas.append(adata)
    
    
    adata = sc.concat(adatas)
    
    adata_neurons = adata[~adata.obs['class'].isin(['30 Astro-Epen',
     '31 OPC-Oligo',
     '32 OEC',
     '33 Vascular',
     '34 Immune'])].copy()
    
    
    neighbors = load_neighbors("allen_original/division_neighbors.txt")
    neigbouring_struct = neighbors[INPUT_REGION]
    
    adata_neurons = adata_neurons[((adata_neurons.obs.division==INPUT_REGION)|(adata_neurons.obs.structure.isin(neigbouring_struct)))].copy()
    
    if write_reference == True:
        adata_neurons.write_h5ad(f'data/reference_dataset_{INPUT_REGION}.h5ad')
    
    
    return adata_neurons
