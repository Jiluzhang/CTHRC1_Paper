## workdir: /fs/home/jiluzhang/xiaoqing20231204/Science/human_sc

# scp -P 6122 scRNA-1_EmptyDrops_CR_matrix.tar jiluzhang@10.11.41.108:/fs/home/jiluzhang/xiaoqing20231204/Science/human_sc
# scp -P 6122 scRNA-2_EmptyDrops_CR_matrix.tar jiluzhang@10.11.41.108:/fs/home/jiluzhang/xiaoqing20231204/Science/human_sc
# scp -P 6122 scRNA1210_EmptyDrops_CR_matrix.tar jiluzhang@10.11.41.108:/fs/home/jiluzhang/xiaoqing20231204/Science/human_sc
# scp -P 6122 scRNA1224-2_EmptyDrops_CR_matrix.tar jiluzhang@10.11.41.108:/fs/home/jiluzhang/xiaoqing20231204/Science/human_sc

## scRNA1210: wt-1
## scRNA1224-2: wt-2
## scRNA-1: ex-1
## scRNA-2: ex-2

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import harmonypy as hm
from matplotlib import rcParams

sc.settings.set_figure_params(dpi=600)
plt.rcParams['pdf.fonttype'] = 42

wt_1 = sc.read_10x_mtx('./raw_data/wt_1')   # 15902 × 38606
wt_2 = sc.read_10x_mtx('./raw_data/wt_2')   # 13299 × 38606
ex_1 = sc.read_10x_mtx('./raw_data/ex_1')   # 14399 × 38606
ex_2 = sc.read_10x_mtx('./raw_data/ex_2')   # 17808 × 38606

wt_1.obs_names = ['wt-1_'+i for i in wt_1.obs_names]
wt_2.obs_names = ['wt-2_'+i for i in wt_2.obs_names]
ex_1.obs_names = ['ex-1_'+i for i in ex_1.obs_names]
ex_2.obs_names = ['ex-2_'+i for i in ex_2.obs_names]

wt_1_gene_info = wt_1.var.copy()
wt_1_gene_info['gene_name'] = wt_1_gene_info.index.values
wt_1_gene_info.reset_index(drop=True, inplace=True)

wt_2_gene_info = wt_2.var.copy()
wt_2_gene_info['gene_name'] = wt_2_gene_info.index.values
wt_2_gene_info.reset_index(drop=True, inplace=True)

ex_1_gene_info = ex_1.var.copy()
ex_1_gene_info['gene_name'] = ex_1_gene_info.index.values
ex_1_gene_info.reset_index(drop=True, inplace=True)

ex_2_gene_info = ex_2.var.copy()
ex_2_gene_info['gene_name'] = ex_2_gene_info.index.values
ex_2_gene_info.reset_index(drop=True, inplace=True)

gene_info = pd.concat([wt_1_gene_info, wt_2_gene_info, ex_1_gene_info, ex_2_gene_info])
gene_info.drop_duplicates(inplace=True)  # 38606*2
gene_info.index = gene_info['gene_name'].values

adata = ad.concat([wt_1, wt_2, ex_1, ex_2], join='outer')   # 61408 × 38606
adata.var['gene_ids'] = gene_info.loc[adata.var.index]['gene_ids'].values

sc.pp.filter_cells(adata, min_genes=500)     # 52723 × 38606
sc.pp.filter_genes(adata, min_cells=10)      # 52723 × 25929

adata.obs['sample_groups'] = [i.split('_')[0] for i in adata.obs.index]
adata.obs['wt_or_ex'] = [i.split('-')[0] for i in adata.obs['sample_groups']]

## qc for mt genes
adata.var['MT'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['MT'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_MT<10, :]  # 45501 × 25929

## normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

## hvg identification
sc.pp.highly_variable_genes(adata, batch_key='sample_groups')
adata.raw = adata
adata = adata[:, adata.var.highly_variable]  # 45501 × 3022

## plot qc metrics
# n_genes_by_counts: the number of genes expressed in the count matrix
sc.pl.violin(adata, keys='n_genes_by_counts', groupby='sample_groups', rotation=90, stripplot=False,
             order=['wt-1', 'wt-2', 'ex-1', 'ex-2'], save='_n_genes.pdf')

# total_counts: the total counts per cell
sc.pl.violin(adata, keys='total_counts', groupby='sample_groups', rotation=90, stripplot=False,
             order=['wt-1', 'wt-2', 'ex-1', 'ex-2'], save='_nCounts.pdf')

# pct_counts_mt: the percentage of counts in mitochondrial genes
sc.pl.violin(adata, keys='pct_counts_MT', groupby='sample_groups', rotation=90, stripplot=False,
             order=['wt-1', 'wt-2', 'ex-1', 'ex-2'], save='_pct_mt.pdf')

# actb expression level
sc.pl.violin(adata, keys='ACTB', groupby='sample_groups', rotation=90, stripplot=False,
             order=['wt-1', 'wt-2', 'ex-1', 'ex-2'], save='_actb_exp.pdf')

## PCA
sc.tl.pca(adata)

## Clustering
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'sample_groups')  # Converged after 5 iterations
adata.obsm['X_pca_harmony'] = ho.Z_corr.T
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.pl.umap(adata, color='sample_groups', legend_fontsize='5', legend_loc='right margin', size=2, 
           title='', frameon=True, save='_batch_harmony_sample_groups.pdf')
# sc.pl.umap(adata, color='wt_or_ex', legend_fontsize='5', legend_loc='right margin', size=2, 
#            title='', frameon=True, save='_batch_harmony_wt_or_ex.pdf')

sc.tl.leiden(adata, resolution=0.5)
adata.obs['leiden'].value_counts()
# 0     7501
# 1     6942
# 2     6304
# 3     4821
# 4     4702
# 5     3477
# 6     3017
# 7     1831
# 8     1619
# 9     1350
# 10    1342
# 11    1231
# 12     619
# 13     402
# 14     172
# 15     171

sc.pl.umap(adata, color='leiden', legend_fontsize='small', legend_loc='on data',
           title='', frameon=True, save='_leiden.pdf')

adata.write('res.h5ad')


marker_genes = ['TAGLN', 'MYH11', 'ACTA2', 'MYL9', 'TPM2',
                'DCN', 'GSN', 'DPT', 'COL1A2', 'IGFBP6',
                'PECAM1', 'AQP1', 'VWF', 'FABP4', 'EGFL7',
                'LYZ2', 'C1QA', 'C1QB', 'C1QC', 'H2-AA',
                'MPZ', 'MBP', 'PLP1', 'NCMAP', 'KCNA1',
                'IL7R', 'CD3G', 'CD52', 'CCL5', 'IL2RB',
                'CD74', 'IGKC', 'CD79A', 'IGLC2', 'MS4A1',
                'GPM6B', 'DBI', 'PRMP', 'SCN7A',
                'BPIFA1', 'SCGB1A1', 'WFDC2', 'REG3G', 'SCGB3A1',
                'CXCR2', 'IL1R2', 'CD177', 'MMP9', 'S100A9',
                'MKI67', 'TOP2A',
                'GYPA', 'ALAS2', 'HBA-A1', 'HBA-A2', 'HBB']
# KeyError: "Could not find keys '['ALAS2', 'BPIFA1', 'GYPA', 'H2-AA', 'HBA-A1', 'HBA-A2', 'LYZ2', 'PRMP', 'REG3G']' in columns of `adata.obs` or in adata.raw.var_names."

marker_genes = ['TAGLN', 'MYH11', 'ACTA2', 'MYL9', 'TPM2',
                'DCN', 'GSN', 'DPT', 'COL1A2', 'IGFBP6',
                'PECAM1', 'AQP1', 'VWF', 'FABP4', 'EGFL7',
                'C1QA', 'C1QB', 'C1QC',
                'MPZ', 'MBP', 'PLP1', 'NCMAP', 'KCNA1',
                'IL7R', 'CD3G', 'CD52', 'CCL5', 'IL2RB',
                'CD74', 'IGKC', 'CD79A', 'IGLC2', 'MS4A1',
                'GPM6B', 'DBI', 'SCN7A',
                'SCGB1A1', 'WFDC2', 'SCGB3A1',
                'CXCR2', 'IL1R2', 'CD177', 'MMP9', 'S100A9',
                'MKI67', 'TOP2A',
                'HBB']

sc.pl.dotplot(adata, marker_genes, groupby='leiden', dendrogram=True, standard_scale='var',
              save='marker_genes.pdf')

sc.tl.rank_genes_groups(adata, 'leiden')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20).to_csv('cluster_marker_genes.txt', index=False, sep='\t')

cell_anno_dict = {'0':'MoMaphDC', '1':'EC', '2':'MoMaphDC', '3':'Fibroblast', '4':'SMC', '5':'MoMaphDC',
                  '6':'T cell', '7':'MoMaphDC', '8':'MoMaphDC', '9':'EC', '10':'Neutrophil', '11':'SMC',
                  '12':'Fibroblast', '13':'Mast cell', '14':'MoMaphDC', '15':'EC'}
adata.obs['cell_anno'] = adata.obs['leiden']
adata.obs['cell_anno'].replace(cell_anno_dict, inplace=True)

# rcParams["figure.figsize"] = (2, 2)
sc.pl.umap(adata, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin',
           title='', frameon=True, save='_cell_anno.pdf')

marker_genes = ['TAGLN', 'MYH11', 'ACTA2', 'MYL9', 'TPM2',
                'DCN', 'GSN', 'DPT', 'COL1A2', 'IGFBP6',
                'PECAM1', 'AQP1', 'VWF', 'FABP4', 'EGFL7',
                'C1QA', 'C1QB', 'C1QC',
                'IL7R', 'CD3G', 'CD52', 'CCL5', 'IL2RB',
                'CXCR2', 'IL1R2', 'CD177', 'MMP9', 'S100A9',
                'KIT', 'CPA3', 'TPSB2']  # add marker genes of mast cell
adata.obs['cell_anno'] = pd.Categorical(adata.obs['cell_anno'], categories=['SMC', 'Fibroblast', 'EC', 'MoMaphDC', 'T cell',
                                                                            'Neutrophil', 'Mast cell'])
sc.pl.matrixplot(adata, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='cell_anno_heatmap.pdf')

adata[adata.obs['wt_or_ex']=='wt'].obs['cell_anno'].value_counts()
# MoMaphDC      7857
# EC            6010
# SMC           4685
# Fibroblast    3799
# T cell        1651
# Mast cell      202
# Neutrophil      89

adata[adata.obs['wt_or_ex']=='ex'].obs['cell_anno'].value_counts()
# MoMaphDC      13047
# EC             2453
# Fibroblast     1641
# T cell         1366
# Neutrophil     1253
# SMC            1248
# Mast cell       200

adata.write('res.h5ad')

## plot marker genes in UMAP
from matplotlib.colors import LinearSegmentedColormap

sc.pl.umap(adata, color=['TAGLN', 'DCN', 'PECAM1', 'C1QA', 'IL7R', 'CXCR2', 'KIT'],
           cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_loc="on data", frameon=True, ncols=3, save='_marker_genes.pdf')

#### plot bar for cell proportion
wt_df = (adata[adata.obs['wt_or_ex']=='wt'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
wt_df['wt_or_ex'] = 'wt'
ex_df = (adata[adata.obs['wt_or_ex']=='ex'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ex_df['wt_or_ex'] = 'ex'
df = pd.concat([wt_df, ex_df])
df.reset_index(drop=True)
df.to_csv('cell_proportion_wt_ex.txt', index=False, sep='\t')

import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('cell_proportion_wt_ex.txt')
df['cell_anno'] = pd.Categorical(df['cell_anno'], categories=['SMC', 'Fibroblast', 'EC', 'MoMaphDC', 'T cell',
                                                              'Neutrophil', 'Mast cell'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='wt_or_ex', y='proportion', fill='cell_anno')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') + coord_flip() +\
                                                                      scale_y_continuous(limits=[0, 1], breaks=np.arange(0, 1+0.1, 0.2)) + theme_bw()  
# scale_fill_brewer(type='qualitative',palette='Set1') +\
p.save(filename='cell_proportion_wt_ex.pdf', dpi=600, height=3, width=6)


#### Fibroblast
adata = sc.read_h5ad('res.h5ad')
fibro = adata[adata.obs['cell_anno']=='Fibroblast'].copy()   # 5440 × 3022

# sc.pl.umap(fibro, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1,
#            title='', frameon=True, save='_fibro_cell_anno.pdf')
# sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='wt'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1,
#            title='', frameon=True, save='_fibro_cell_anno_wt.pdf')
# sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='ex'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1,
#            title='', frameon=True, save='_fibro_cell_anno_ex.pdf')

## CTHRC1 not in the top marker gene list
# sc.tl.leiden(fibro, resolution=0.1)
# fibro.obs['leiden'].value_counts()
# # 0    3732
# # 1    1089
# # 2     619

# fibro[fibro.obs['wt_or_ex']=='wt'].obs['leiden'].value_counts()
# # 0    2693
# # 1     863
# # 2     243

# fibro[fibro.obs['wt_or_ex']=='ex'].obs['leiden'].value_counts()
# # 0    1039
# # 2     376
# # 1     226

# sc.pl.umap(fibro, color='leiden', legend_fontsize='small', legend_loc='on data', size=2,
#            title='', frameon=True, save='_fibro_leiden.pdf')

# sc.tl.rank_genes_groups(fibro, 'leiden')
# pd.DataFrame(fibro.uns['rank_genes_groups']['names']).head(10)
# #           0         1       2
# # 0     FBLN1     ADIRF   CDH19
# # 1    CCDC80    IGFBP2    APOD
# # 2        C7     MFGE8   TGFBI
# # 3     SFRP2       FN1    NGFR
# # 4  SERPINF1      DSTN   CLDN1
# # 5       C1R     MYH11    NRP2
# # 6       DCN  PPP1R14A  PHLDA1
# # 7       C1S       BGN  CYP1B1
# # 8        C3     TAGLN     VIT
# # 9       CFD     CRIP1    PLK2

# del fibro.uns['dendrogram_leiden']
# sc.pl.rank_genes_groups_dotplot(fibro, groupby="leiden", n_genes=10, cmap='bwr', values_to_plot='logfoldchanges',
#                                 vmin=-4, vmax=4, save='fibro_marker_gene.pdf') # min_logfoldchange=2

# fibro.obs['cell_anno'] = fibro.obs['leiden']
# fibro.obs['cell_anno'].replace({'0':'fibro_1', '1':'fibro_2', '2':'fibro_3'}, inplace=True)  # to highlight cluster 1
# fibro.obs['cell_anno'].value_counts()
# # fibro_1    5810
# # fibro_2    2239
# # fibro_3     994

# marker_genes = ['Cthrc1', 'Postn', 'Cxcl2', 'Col8a1', 'Thbs2',
#                 'Fmo2', 'Serping1', 'Cst3', 'Clec3b', 'Plpp3',
#                 'Igkc', 'Ddx5', 'Cd74', 'Rps27', 'Eif1']
# fibro.obs['cell_anno'] = pd.Categorical(fibro.obs['cell_anno'], categories=['fibro_1', 'fibro_2', 'fibro_3'])
# plt.rcParams['pdf.fonttype'] = 42
# sc.pl.matrixplot(fibro, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='fibro_cell_anno_heatmap.pdf')
# sc.pl.matrixplot(fibro, ['CTHRC1'], groupby='leiden', cmap='coolwarm', standard_scale='var', save='fibro_cell_anno_heatmap.pdf')

# ## wt & ex
# rcParams["figure.figsize"] = (3, 3)
# sc.pl.umap(fibro, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
#            title='', frameon=True)
# plt.xlim(1, 12)
# plt.ylim(-8, 3) 
# plt.savefig('./figures/umap_fibro_cell_anno.pdf')
# plt.close()

# ## wt
# sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='wt'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
#            title='', frameon=True)
# plt.xlim(1, 12)
# plt.ylim(-8, 3) 
# plt.savefig('./figures/umap_fibro_cell_anno_wt.pdf')
# plt.close()

# ## ex
# sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='ex'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
#            title='', frameon=True)
# plt.xlim(1, 12)
# plt.ylim(-8, 3) 
# plt.savefig('./figures/umap_fibro_cell_anno_ex.pdf')
# plt.close()


## plot umap for CTHRC1
from matplotlib.colors import LinearSegmentedColormap

sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='wt'], color=['CTHRC1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=10, vmin=0, vmax=5, title='', frameon=True)
plt.xlim(-4, 7)
plt.ylim(-7, 1) 
plt.savefig('./figures/umap_fibro_cthrc1_wt.pdf')
plt.close()

sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='ex'], color=['CTHRC1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=10, vmin=0, vmax=5, title='', frameon=True)
plt.xlim(-4, 7)
plt.ylim(-7, 1) 
plt.savefig('./figures/umap_fibro_cthrc1_ex.pdf')
plt.close()

fibro.write('res_fibro.h5ad')


#### smc subclusters
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
sc.settings.set_figure_params(dpi=600)

adata = sc.read_h5ad('res.h5ad')
smc = adata[adata.obs['cell_anno']=='SMC'].copy()  # 5933 × 3022

sc.tl.leiden(smc, resolution=0.5)
sc.pl.umap(smc, color='leiden', legend_fontsize='medium', legend_loc='right margin', size=10, 
           title='', frameon=True, save='_smc_leiden.pdf')

smc.obs['leiden'].value_counts()
# 0    1528
# 1    1277
# 2    1037
# 3     742
# 4     690
# 5     541
# 6     118

marker_genes = ['ACTA2', 'ACTG1', 'TAGLN2', 'MYL12A', 'TUBB2A',
                'CCNL2', 'CCNT2', 'MALAT1', 'CCN1',
                'COL1A2','LUM', 'DCN', 'GSN', 'CYGB',
                'NOTCH3', 'ATF3', 'RGS4', 'CACNA1H', 'PLN']  

del smc.uns['dendrogram_leiden']
sc.pl.dotplot(smc, marker_genes, groupby='leiden', standard_scale='var', dendrogram=True,
              save='smc_marker_genes.pdf')

cell_anno_dict = {'0':'Fibromyocyte', '1':'Contractile SMC', '2':'Contractile SMC',
                  '3':'Stressed SMC', '4':'Proliferating SMC', '5':'Contractile SMC',
                  '6':'Contractile SMC'}
smc.obs['cell_anno'] = smc.obs['leiden']
smc.obs['cell_anno'].replace(cell_anno_dict, inplace=True)

# sc.pl.umap(smc, color='cell_anno', legend_fontsize='5', legend_loc='right margin', size=10,
#            title='', frameon=True, save='_smc_cell_anno.pdf')

sc.pl.umap(smc, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=20,
           title='', frameon=True)
plt.xlim(-8, 2)
plt.ylim(-6, 3) 
plt.savefig('./figures/umap_smc_cell_anno.pdf')
plt.close()

sc.pl.umap(smc[smc.obs['wt_or_ex']=='wt'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=20,
           title='', frameon=True)
plt.xlim(-8, 2)
plt.ylim(-6, 3) 
plt.savefig('./figures/umap_smc_cell_anno_wt.pdf')
plt.close()

sc.pl.umap(smc[smc.obs['wt_or_ex']=='ex'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=20,
           title='', frameon=True)
plt.xlim(-8, 2)
plt.ylim(-6, 3) 
plt.savefig('./figures/umap_smc_cell_anno_ex.pdf')
plt.close()

smc.obs['cell_anno'] = pd.Categorical(smc.obs['cell_anno'],
                                      categories=['Contractile SMC', 'Proliferating SMC', 'Fibromyocyte', 'Stressed SMC'])
marker_genes = ['ACTA2', 'MYL12A',
                'CCNL2', 'MALAT1',
                'DCN', 'GSN',
                'ATF3', 'CACNA1H']             
sc.pl.matrixplot(smc, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='smc_cell_anno_heatmap.pdf')

## output marker genes for each cluster
sc.tl.rank_genes_groups(smc, 'leiden')
pd.DataFrame(smc.uns['rank_genes_groups']['names']).head(20).to_csv('smc_cluster_marker_genes.txt', index=False, sep='\t')

smc.write('./smc_res.h5ad')

## smc proportion
wt_smc_cnt = sum(smc.obs['wt_or_ex']=='wt')   # 4685
ex_smc_cnt = sum(smc.obs['wt_or_ex']=='ex')   # 1248

wt_smc_con_cnt = sum((smc.obs['cell_anno']=='Contractile SMC') & (smc.obs['wt_or_ex']=='wt'))     # 2479
ex_smc_con_cnt = sum((smc.obs['cell_anno']=='Contractile SMC') & (smc.obs['wt_or_ex']=='ex'))     # 494

wt_smc_pro_cnt = sum((smc.obs['cell_anno']=='Proliferating SMC') & (smc.obs['wt_or_ex']=='wt'))   # 303
ex_smc_pro_cnt = sum((smc.obs['cell_anno']=='Proliferating SMC') & (smc.obs['wt_or_ex']=='ex'))   # 387

wt_smc_str_cnt = sum((smc.obs['cell_anno']=='Stressed SMC') & (smc.obs['wt_or_ex']=='wt'))   # 494
ex_smc_str_cnt = sum((smc.obs['cell_anno']=='Stressed SMC') & (smc.obs['wt_or_ex']=='ex'))   # 248

wt_smc_fib_cnt = sum((smc.obs['cell_anno']=='Fibromyocyte') & (smc.obs['wt_or_ex']=='wt'))   # 1409
ex_smc_fib_cnt = sum((smc.obs['cell_anno']=='Fibromyocyte') & (smc.obs['wt_or_ex']=='ex'))   # 119

smc_df = pd.DataFrame({'smc': ['Contractile SMC']*2+['Proliferating SMC']*2+['Stressed SMC']*2+['Fibromyocyte']*2,
                       'trt': ['wt', 'ex']*4,
                       'val': 0})
smc_df.loc[(smc_df['smc']=='Contractile SMC') & (smc_df['trt']=='wt'), 'val'] = wt_smc_con_cnt/wt_smc_cnt
smc_df.loc[(smc_df['smc']=='Contractile SMC') & (smc_df['trt']=='ex'), 'val'] = ex_smc_con_cnt/ex_smc_cnt
smc_df.loc[(smc_df['smc']=='Proliferating SMC') & (smc_df['trt']=='wt'), 'val'] = wt_smc_pro_cnt/wt_smc_cnt
smc_df.loc[(smc_df['smc']=='Proliferating SMC') & (smc_df['trt']=='ex'), 'val'] = ex_smc_pro_cnt/ex_smc_cnt
smc_df.loc[(smc_df['smc']=='Stressed SMC') & (smc_df['trt']=='wt'), 'val'] = wt_smc_str_cnt/wt_smc_cnt
smc_df.loc[(smc_df['smc']=='Stressed SMC') & (smc_df['trt']=='ex'), 'val'] = ex_smc_str_cnt/ex_smc_cnt
smc_df.loc[(smc_df['smc']=='Fibromyocyte') & (smc_df['trt']=='wt'), 'val'] = wt_smc_fib_cnt/wt_smc_cnt
smc_df.loc[(smc_df['smc']=='Fibromyocyte') & (smc_df['trt']=='ex'), 'val'] = ex_smc_fib_cnt/ex_smc_cnt
smc_df.to_csv('smc_subcluster_proportion.txt', index=False, sep='\t')

smc.obs['cell_anno'].value_counts()
# Contractile SMC      2973
# Fibromyocyte         1528
# Stressed SMC          742
# Proliferating SMC     690

#### plot bar for smc subcluster proportion
import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('smc_subcluster_proportion.txt')
smc = ['Contractile SMC', 'Proliferating SMC', 'Stressed SMC', 'Fibromyocyte']
df['smc'] = pd.Categorical(df['smc'], categories=smc)
df['trt'] = pd.Categorical(df['trt'], categories=['wt', 'ex'])

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='smc', y='val', fill='trt')) + geom_bar(stat='identity', position='dodge', width=0.75) + xlab('') + ylab('Fraction') +\
                                                    scale_y_continuous(limits=[0, 0.6], breaks=np.arange(0, 0.6+0.1, 0.2)) + theme_bw()
p.save(filename='smc_subcluster_proportion.pdf', dpi=600, height=3, width=6)


#### find differential genes
## Contractile SMC
import scanpy as sc
import pandas as pd

smc = sc.read_h5ad('smc_res.h5ad')

smc_con = smc[smc.obs['cell_anno']=='Contractile SMC'].copy()  # 2973 × 3022
sc.tl.rank_genes_groups(smc_con, 'wt_or_ex')
result = smc_con.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_con_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_con_up_cnt = sum((smc_con_res['ex_logfoldchanges']>1) & (smc_con_res['ex_pvals_adj']<0.05))   # 2102
smc_con_dw_cnt = sum((smc_con_res['ex_logfoldchanges']<-1) & (smc_con_res['ex_pvals_adj']<0.05))  # 1418
smc_con_df_cnt = smc_con_up_cnt + smc_con_dw_cnt                                                  # 3520
smc_con_res[(smc_con_res['ex_logfoldchanges']<-1) & (smc_con_res['ex_pvals_adj']<0.05)]['ex_names'].to_csv('smc_con_down_genes.txt', index=False, header=False)
smc_con_res[(smc_con_res['ex_logfoldchanges']>1) & (smc_con_res['ex_pvals_adj']<0.05)]['ex_names'].to_csv('smc_con_up_genes.txt', index=False, header=False)

## Proliferating SMC
smc_pro = smc[smc.obs['cell_anno']=='Proliferating SMC'].copy()  # 690 × 3022
sc.tl.rank_genes_groups(smc_pro, 'wt_or_ex')
result = smc_pro.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_pro_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_pro_up_cnt = sum((smc_pro_res['ex_logfoldchanges']>1) & (smc_pro_res['ex_pvals_adj']<0.05))   # 2357
smc_pro_dw_cnt = sum((smc_pro_res['ex_logfoldchanges']<-1) & (smc_pro_res['ex_pvals_adj']<0.05))  # 558
smc_pro_df_cnt = smc_pro_up_cnt + smc_pro_dw_cnt                                                  # 2915
smc_pro_res[(smc_pro_res['ex_logfoldchanges']<-1) & (smc_pro_res['ex_pvals_adj']<0.05)]['ex_names'].to_csv('smc_pro_down_genes.txt', index=False, header=False)
smc_pro_res[(smc_pro_res['ex_logfoldchanges']>1) & (smc_pro_res['ex_pvals_adj']<0.05)]['ex_names'].to_csv('smc_pro_up_genes.txt', index=False, header=False)

## Stressed SMC
smc_str = smc[smc.obs['cell_anno']=='Stressed SMC'].copy()  # 742 × 3022
sc.tl.rank_genes_groups(smc_str, 'wt_or_ex')
result = smc_str.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_str_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_str_up_cnt = sum((smc_str_res['ex_logfoldchanges']>1) & (smc_str_res['ex_pvals_adj']<0.05))   # 1303
smc_str_dw_cnt = sum((smc_str_res['ex_logfoldchanges']<-1) & (smc_str_res['ex_pvals_adj']<0.05))  # 1051
smc_str_df_cnt = smc_str_up_cnt + smc_str_dw_cnt                                                  # 2354
smc_str_res[(smc_str_res['ex_logfoldchanges']<-1) & (smc_str_res['ex_pvals_adj']<0.05)]['ex_names'].to_csv('smc_str_down_genes.txt', index=False, header=False)
smc_str_res[(smc_str_res['ex_logfoldchanges']>1) & (smc_str_res['ex_pvals_adj']<0.05)]['ex_names'].to_csv('smc_str_up_genes.txt', index=False, header=False)

## Fibromyocyte
smc_fib = smc[smc.obs['cell_anno']=='Fibromyocyte'].copy()  # 1528 × 3022
sc.tl.rank_genes_groups(smc_fib, 'wt_or_ex')
result = smc_fib.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_fib_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_fib_up_cnt = sum((smc_fib_res['ex_logfoldchanges']>1) & (smc_fib_res['ex_pvals_adj']<0.05))   # 788
smc_fib_dw_cnt = sum((smc_fib_res['ex_logfoldchanges']<-1) & (smc_fib_res['ex_pvals_adj']<0.05))  # 2013
smc_fib_df_cnt = smc_fib_up_cnt + smc_fib_dw_cnt                                                  # 2801
smc_fib_res[(smc_fib_res['ex_logfoldchanges']<-1) & (smc_fib_res['ex_pvals_adj']<0.05)]['ex_names'].to_csv('smc_fib_down_genes.txt', index=False, header=False)
smc_fib_res[(smc_fib_res['ex_logfoldchanges']>1) & (smc_fib_res['ex_pvals_adj']<0.05)]['ex_names'].to_csv('smc_fib_up_genes.txt', index=False, header=False)

smc_df = pd.DataFrame({'smc': ['Contractile SMC']*2+['Proliferating SMC']*2+['Stressed SMC']*2+['Fibromyocyte']*2,
                       'chg': ['up', 'down']*4,
                       'val': 0})
smc_df.loc[(smc_df['smc']=='Contractile SMC') & (smc_df['chg']=='up'), 'val'] = smc_con_up_cnt
smc_df.loc[(smc_df['smc']=='Contractile SMC') & (smc_df['chg']=='down'), 'val'] = smc_con_dw_cnt
smc_df.loc[(smc_df['smc']=='Proliferating SMC') & (smc_df['chg']=='up'), 'val'] = smc_pro_up_cnt
smc_df.loc[(smc_df['smc']=='Proliferating SMC') & (smc_df['chg']=='down'), 'val'] = smc_pro_dw_cnt
smc_df.loc[(smc_df['smc']=='Stressed SMC') & (smc_df['chg']=='up'), 'val'] = smc_str_up_cnt
smc_df.loc[(smc_df['smc']=='Stressed SMC') & (smc_df['chg']=='down'), 'val'] = smc_str_dw_cnt
smc_df.loc[(smc_df['smc']=='Fibromyocyte') & (smc_df['chg']=='up'), 'val'] = smc_fib_up_cnt
smc_df.loc[(smc_df['smc']=='Fibromyocyte') & (smc_df['chg']=='down'), 'val'] = smc_fib_dw_cnt
smc_df.to_csv('smc_subcluster_deg_cnt.txt', index=False, sep='\t')


#### plot bar for smc subcluster DEG count
import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('smc_subcluster_deg_cnt.txt')
smc = ['Contractile SMC', 'Proliferating SMC', 'Stressed SMC', 'Fibromyocyte']
df['smc'] = pd.Categorical(df['smc'], categories=smc)
chg = ['up', 'down']
df['chg'] = pd.Categorical(df['chg'], categories=chg)

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='smc', y='val', fill='chg')) + geom_bar(stat='identity', position=position_dodge(), width=0.75) +\
                                                    xlab('') + ylab('DEG Count') +\
                                                    scale_y_continuous(limits=[0, 2500], breaks=np.arange(0, 2500+1, 500)) + theme_bw()
p.save(filename='smc_subcluster_deg_cnt.pdf', dpi=600, height=3, width=5)


#### GO enrichment analysis for DEGs
library(clusterProfiler)
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

## smc_con_up
smc_con_up <- read.table('smc_con_up_genes.txt')[, 1]
smc_con_up_ego <- enrichGO(gene=smc_con_up, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_con_up_genes_go.pdf')
dotplot(smc_con_up_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_con_up_ego), file="smc_con_up_genes_go.csv")

## smc_con_down
smc_con_down <- read.table('smc_con_down_genes.txt')[, 1]
smc_con_down_ego <- enrichGO(gene=smc_con_down, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_con_down_genes_go.pdf')
dotplot(smc_con_down_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_con_down_ego), file="smc_con_down_genes_go.csv")

terms <- c('GO:0003012', 'GO:0006936', 'GO:0006937', 'GO:0006941', 'GO:0090257',
           'GO:0030239', 'GO:0030048', 'GO:0070252', 'GO:1903115', 'GO:0006940')
pdf('smc_con_down_genes_go_selected.pdf')
dotplot(smc_con_down_ego, showCategory=smc_con_down_ego@result$Description[smc_con_down_ego@result$ID %in% terms], font.size=10)
dev.off()

## smc_pro_up
smc_pro_up <- read.table('smc_pro_up_genes.txt')[, 1]
smc_pro_up_ego <- enrichGO(gene=smc_pro_up, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_pro_up_genes_go.pdf')
dotplot(smc_pro_up_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_pro_up_ego), file="smc_pro_up_genes_go.csv")

## smc_pro_down
smc_pro_down <- read.table('smc_pro_down_genes.txt')[, 1]
smc_pro_down_ego <- enrichGO(gene=smc_pro_down, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_pro_down_genes_go.pdf')
dotplot(smc_pro_down_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_pro_down_ego), file="smc_pro_down_genes_go.csv")

## smc_str_up
smc_str_up <- read.table('smc_str_up_genes.txt')[, 1]
smc_str_up_ego <- enrichGO(gene=smc_str_up, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_str_up_genes_go.pdf')
dotplot(smc_str_up_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_str_up_ego), file="smc_str_up_genes_go.csv")

## smc_str_down
smc_str_down <- read.table('smc_str_down_genes.txt')[, 1]
smc_str_down_ego <- enrichGO(gene=smc_str_down, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_str_down_genes_go.pdf')
dotplot(smc_str_down_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_str_down_ego), file="smc_str_down_genes_go.csv")

## smc_fib_up
smc_fib_up <- read.table('smc_fib_up_genes.txt')[, 1]
smc_fib_up_ego <- enrichGO(gene=smc_fib_up, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_fib_up_genes_go.pdf')
dotplot(smc_fib_up_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_fib_up_ego), file="smc_fib_up_genes_go.csv")

terms <- c('GO:0097193', 'GO:0032103', 'GO:0052547', 'GO:0001819', 'GO:0050729',
           'GO:0045862', 'GO:0030595', 'GO:0001649', 'GO:0070371', 'GO:0048661')
pdf('smc_fib_up_genes_go_selected.pdf')
dotplot(smc_fib_up_ego, showCategory=smc_fib_up_ego@result$Description[smc_fib_up_ego@result$ID %in% terms], font.size=10)
dev.off()

## smc_fib_down
smc_fib_down <- read.table('smc_fib_down_genes.txt')[, 1]
smc_fib_down_ego <- enrichGO(gene=smc_fib_down, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_fib_down_genes_go.pdf')
dotplot(smc_fib_down_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_fib_down_ego), file="smc_fib_down_genes_go.csv")

# terms <- c('GO:0032103', 'GO:0043254', 'GO:0045862', 'GO:0031589', 'GO:0045860',
#            'GO:0001503', 'GO:0007160', 'GO:0014909', 'GO:0030198', 'GO:0051403')
# pdf('smc_fib_down_genes_go_selected.pdf')
# dotplot(smc_fib_down_ego, showCategory=smc_fib_down_ego@result$Description[smc_fib_down_ego@result$ID %in% terms], font.size=10)
# dev.off()


#### plot box for specific geneset
from plotnine import *
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

smc = sc.read_h5ad('smc_res.h5ad')  # 5933 × 3022

## contractile genes
con_genes = ['RAMP1', 'CNN1', 'ACTG2', 'TAGLN', 'ACTA2', 'TPM2', 'MYL6', 'MYL9', 'MYH11', 'MYLK']  # delete Mylk2 due to low expression
dat = pd.DataFrame()
for gene in con_genes:
    dat_wt = pd.DataFrame({'gene_id':gene,
                           'idx':'wt',
                           'val':smc[smc.obs['wt_or_ex']=='wt'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat_ex = pd.DataFrame({'gene_id':gene,
                           'idx':'ex',
                           'val':smc[smc.obs['wt_or_ex']=='ex'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat = pd.concat([dat, dat_wt, dat_ex])

dat['idx'] = pd.Categorical(dat['idx'], categories=['wt', 'ex'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(dat, aes(x='gene_id', y='val', fill='idx')) + geom_boxplot(width=0.5, outlier_shape='', show_legend=False) +\
                                                         scale_y_continuous(limits=[0, 7], breaks=np.arange(0, 7+0.1, 1.0)) + theme_bw()
p.save(filename='smc_contraction_seperately.pdf', dpi=600, height=4, width=5)

for gene in con_genes:
    print(gene, 'pvalue:', stats.ttest_ind(dat[(dat['idx']=='wt') & (dat['gene_id']==gene)]['val'],
                                           dat[(dat['idx']=='ex') & (dat['gene_id']==gene)]['val'])[1])
# RAMP1 pvalue: 3.9442636121081895e-78
# CNN1 pvalue: 2.0017443201907703e-233
# ACTG2 pvalue: 1.8810419558685794e-47
# TAGLN pvalue: 1.794799943065285e-244
# ACTA2 pvalue: 0.0
# TPM2 pvalue: 1.6206710770337024e-113
# MYL6 pvalue: 4.200715012764694e-59
# MYL9 pvalue: 0.0
# MYH11 pvalue: 0.0
# MYLK pvalue: 2.2644993961937348e-35

## mmp genes (median is all zero!!!)
#mmp_genes = ['Mmp14', 'Mmp2']  # delete Mmp24 & Mmp8 & Mmp9 & Mmp15 & Mmp13 due to low expression
mmp_genes = ['MMP14', 'MMP2', 'MMP24', 'MMP8', 'MMP9', 'MMP15', 'MMP13']
dat = pd.DataFrame()
for gene in mmp_genes:
    dat_wt = pd.DataFrame({'gene_id':gene,
                           'idx':'wt',
                           'val':smc[smc.obs['wt_or_ex']=='wt'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat_ex = pd.DataFrame({'gene_id':gene,
                           'idx':'ex',
                           'val':smc[smc.obs['wt_or_ex']=='ex'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat = pd.concat([dat, dat_wt, dat_ex])

dat['idx'] = pd.Categorical(dat['idx'], categories=['wt', 'ex'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(dat, aes(x='gene_id', y='val', fill='idx')) + geom_boxplot(width=0.5, outlier_shape='', show_legend=False) +\
                                                         scale_y_continuous(limits=[0, 4], breaks=np.arange(0, 4+0.1, 1)) + theme_bw()
p.save(filename='smc_mmp_seperately.pdf', dpi=600, height=4, width=5)

for gene in mmp_genes:
    print(gene, 'pvalue:', stats.ttest_ind(dat[(dat['idx']=='wt') & (dat['gene_id']==gene)]['val'],
                                           dat[(dat['idx']=='ex') & (dat['gene_id']==gene)]['val'])[1])
# MMP14 pvalue: 0.013659190456015793
# MMP2 pvalue: 0.19792162168259925
# MMP24 pvalue: 0.041169598125871336
# MMP8 pvalue: 0.05267359089000932
# MMP9 pvalue: 0.0001019132944538918
# MMP15 pvalue: 0.7634441949198965
# MMP13 pvalue: nan


## CTHRC1 in fibroblasts
from plotnine import *
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
plt.rcParams['pdf.fonttype'] = 42

adata = sc.read_h5ad('res.h5ad')  
fib = adata[adata.obs['cell_anno']=='Fibroblast', :].copy()  # 5440 × 3022

dat_wt_1 = pd.DataFrame({'gene_id':'CTHRC1', 'idx':'wt-1',
                       'val':fib[fib.obs['sample_groups']=='wt-1'].raw.X[:, np.argwhere(fib.raw.var.index=='CTHRC1').flatten()].toarray().flatten()})
dat_wt_2 = pd.DataFrame({'gene_id':'CTHRC1', 'idx':'wt-2',
                       'val':fib[fib.obs['sample_groups']=='wt-2'].raw.X[:, np.argwhere(fib.raw.var.index=='CTHRC1').flatten()].toarray().flatten()})
dat_ex_1 = pd.DataFrame({'gene_id':'CTHRC1', 'idx':'ex-1',
                       'val':fib[fib.obs['sample_groups']=='ex-1'].raw.X[:, np.argwhere(fib.raw.var.index=='CTHRC1').flatten()].toarray().flatten()})
dat_ex_2 = pd.DataFrame({'gene_id':'CTHRC1', 'idx':'ex-2',
                       'val':fib[fib.obs['sample_groups']=='ex-2'].raw.X[:, np.argwhere(fib.raw.var.index=='CTHRC1').flatten()].toarray().flatten()})
dat = pd.concat([dat_wt_1, dat_wt_2, dat_ex_1, dat_ex_2])

dat = dat[dat['idx']!='ex-1']  # remove 'ex-1'
dat['idx'].replace({'wt-1':'wt', 'wt-2':'wt', 'ex-2':'ex'}, inplace=True)
dat['idx'] = pd.Categorical(dat['idx'], categories=['wt', 'ex'])

# all control cells
group_means= dat.groupby('idx')['val'].mean().reset_index()
p = ggplot(dat, aes(x='idx', y='val')) + geom_jitter(aes(color='idx'), width=0.3, size=0.75, alpha=0.75) + geom_hline(aes(yintercept='val', color='idx'), data=group_means) +\
                                         scale_y_continuous(limits=[0, 4], breaks=np.arange(0, 4+0.1, 1)) + theme_bw()
p.save(filename='fib_cthrc1_seperately.pdf', dpi=600, height=4, width=5)


stats.ttest_ind(dat[(dat['idx']=='wt') & (dat['gene_id']=='CTHRC1')]['val'], 
                dat[(dat['idx']=='ex') & (dat['gene_id']=='CTHRC1')]['val'])[1]
# 2.525439866778745e-05

# randomly selected control cells
dat.reset_index(drop=True, inplace=True)

np.random.seed(0)
df = pd.concat([dat.loc[np.random.choice(dat[dat['idx']=='wt'].index, sum(dat['idx']=='ex'), replace=False)], 
                dat[dat['idx']=='ex']])

group_means= df.groupby('idx')['val'].mean().reset_index()
p = ggplot(df, aes(x='idx', y='val')) + geom_jitter(aes(color='idx'), width=0.3, size=0.75, alpha=0.75) + geom_hline(aes(yintercept='val', color='idx'), data=group_means) +\
                                         scale_y_continuous(limits=[0, 4], breaks=np.arange(0, 4+0.1, 1)) + theme_bw()
p.save(filename='fib_cthrc1_seperately_random.pdf', dpi=600, height=4, width=5)

stats.ttest_ind(df[(df['idx']=='wt') & (df['gene_id']=='CTHRC1')]['val'], 
                df[(df['idx']=='ex') & (df['gene_id']=='CTHRC1')]['val'])[1]
# 0.006160499280139001
