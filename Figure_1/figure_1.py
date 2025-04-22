## workdir: /fs/home/jiluzhang/xiaoqing20231204/Science/Figure_1_new

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import harmonypy as hm
from matplotlib import rcParams

sc.settings.set_figure_params(dpi=600)
plt.rcParams['pdf.fonttype'] = 42

wt_1 = sc.read_10x_mtx('./raw_data/wt_1')   # 26642 × 32662
wt_2 = sc.read_10x_mtx('./raw_data/wt_2')   # 25602 × 31969
ex_1 = sc.read_10x_mtx('./raw_data/ex_1')   # 16381 × 29180
ex_2 = sc.read_10x_mtx('./raw_data/ex_2')   # 25547 × 32460

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
gene_info.drop_duplicates(inplace=True)  # 39006*2
gene_info.index = gene_info['gene_name'].values

adata = ad.concat([wt_1, wt_2, ex_1, ex_2], join='outer')   # 94172 × 39006
adata.var['gene_ids'] = gene_info.loc[adata.var.index]['gene_ids'].values

sc.pp.filter_cells(adata, min_genes=500)     # 86738 × 39006
sc.pp.filter_genes(adata, min_cells=10)      # 86738 × 26020

adata.obs['sample_groups'] = [i.split('_')[0] for i in adata.obs.index]
adata.obs['wt_or_ex'] = [i.split('-')[0] for i in adata.obs['sample_groups']]

## qc for mt genes
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt<10, :]  # 82326 × 26020

## normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

## hvg identification
sc.pp.highly_variable_genes(adata, batch_key='sample_groups')
adata.raw = adata
adata = adata[:, adata.var.highly_variable]  # 82326 × 2660

## plot qc metrics
# n_genes_by_counts: the number of genes expressed in the count matrix
sc.pl.violin(adata, keys='n_genes_by_counts', groupby='sample_groups', rotation=90, stripplot=False,
             order=['wt-1', 'wt-2', 'ex-1', 'ex-2'], save='_n_genes.pdf')

# total_counts: the total counts per cell
sc.pl.violin(adata, keys='total_counts', groupby='sample_groups', rotation=90, stripplot=False,
             order=['wt-1', 'wt-2', 'ex-1', 'ex-2'], save='_nCounts.pdf')

# pct_counts_mt: the percentage of counts in mitochondrial genes
sc.pl.violin(adata, keys='pct_counts_mt', groupby='sample_groups', rotation=90, stripplot=False,
             order=['wt-1', 'wt-2', 'ex-1', 'ex-2'], save='_pct_mt.pdf')

# actb expression level
sc.pl.violin(adata, keys='Actb', groupby='sample_groups', rotation=90, stripplot=False,
             order=['wt-1', 'wt-2', 'ex-1', 'ex-2'], save='_actb_exp.pdf')

## PCA
sc.tl.pca(adata)

## Clustering
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'sample_groups')  # Converged after 3 iterations
adata.obsm['X_pca_harmony'] = ho.Z_corr.T
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.pl.umap(adata, color='sample_groups', legend_fontsize='5', legend_loc='right margin', size=2, 
           title='', frameon=True, save='_batch_harmony_sample_groups.pdf')
sc.pl.umap(adata, color='wt_or_ex', legend_fontsize='5', legend_loc='right margin', size=2, 
           title='', frameon=True, save='_batch_harmony_wt_or_ex.pdf')

sc.tl.leiden(adata, resolution=0.5)
adata.obs['leiden'].value_counts()
# 0     19717
# 1     15229
# 2     13014
# 3     10223
# 4      8324
# 5      7068
# 6      1980
# 7      1950
# 8      1299
# 9      1100
# 10      719
# 11      573
# 12      444
# 13      394
# 14      174
# 15      118

sc.pl.umap(adata, color='leiden', legend_fontsize='small', legend_loc='on data',
           title='', frameon=True, save='_leiden.pdf')

adata.write('res.h5ad')


marker_genes = ['Tagln', 'Myh11', 'Acta2', 'Myl9', 'Tpm2',
                'Dcn', 'Gsn', 'Dpt', 'Col1a2', 'Igfbp6',
                'Pecam1', 'Aqp1', 'Vwf', 'Fabp4', 'Egfl7', 
                'Lyz2', 'C1qa', 'C1qb', 'C1qc', 'H2-Aa',
                'Mpz', 'Mbp', 'Plp1', 'Ncmap', 'Kcna1',
                'Il7r', 'Cd3g', 'Cd52', 'Ccl5', 'Il2rb', 
                'Cd74', 'Igkc', 'Cd79a', 'Iglc2', 'Ms4a1', 
                'Gpm6b', 'Dbi', 'Prmp', 'Scn7a', 
                'Bpifa1', 'Scgb1a1', 'Wfdc2', 'Reg3g', 'Scgb3a1',
                'Cxcr2', 'Il1r2', 'Cd177', 'Mmp9', 'S100a9',
                'Mki67', 'Top2a',
                'Gypa', 'Alas2', 'Hba-a1', 'Hba-a2', 'Hbb']
# KeyError: "Could not find keys '['Bpifa1', 'Prmp', 'Reg3g', 'Hbb']' in columns of `adata.obs` or in adata.raw.var_names."

marker_genes = ['Tagln', 'Myh11', 'Acta2', 'Myl9', 'Tpm2',
                'Dcn', 'Gsn', 'Dpt', 'Col1a2', 'Igfbp6',
                'Pecam1', 'Aqp1', 'Vwf', 'Fabp4', 'Egfl7', 
                'Lyz2', 'C1qa', 'C1qb', 'C1qc', 'H2-Aa',
                'Mpz', 'Mbp', 'Plp1', 'Ncmap', 'Kcna1',
                'Il7r', 'Cd3g', 'Cd52', 'Ccl5', 'Il2rb', 
                'Cd74', 'Igkc', 'Cd79a', 'Iglc2', 'Ms4a1', 
                'Gpm6b', 'Dbi', 'Scn7a', 
                'Scgb1a1', 'Wfdc2', 'Scgb3a1',
                'Cxcr2', 'Il1r2', 'Cd177', 'Mmp9', 'S100a9',
                'Mki67', 'Top2a',
                'Gypa', 'Alas2', 'Hba-a1', 'Hba-a2']

sc.pl.dotplot(adata, marker_genes, groupby='leiden', dendrogram=True, standard_scale='var',
              save='marker_genes.pdf')

sc.tl.rank_genes_groups(adata, 'leiden')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20).to_csv('cluster_marker_genes.txt', index=False, sep='\t')

cell_anno_dict = {'0':'B cell', '1':'T cell', '2':'Neutrophil', '3':'MoMaphDC', '4':'Fibroblast', '5':'SMC',
                  '6':'MoMaphDC', '7':'EC', '8':'B cell', '9':'SMC', '10':'Fibroblast', '11':'B cell',
                  '12':'T cell', '13':'Schwann', '14':'Neutrophil', '15':'Erythrocyte'}
adata.obs['cell_anno'] = adata.obs['leiden']
adata.obs['cell_anno'].replace(cell_anno_dict, inplace=True)

# rcParams["figure.figsize"] = (2, 2)
sc.pl.umap(adata, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin',
           title='', frameon=True, save='_cell_anno.png')

marker_genes = ['Tagln', 'Myh11', 'Acta2', 'Myl9', 'Tpm2',
                'Dcn', 'Dpt', 'Col1a2', 'Igfbp6',
                'Pecam1', 'Aqp1', 'Vwf', 'Fabp4', 'Egfl7', 
                'Lyz2', 'C1qa', 'C1qb', 'C1qc',
                'Mpz', 'Mbp', 'Plp1', 'Ncmap', 'Kcna1',
                'Il7r', 'Cd3g', 'Ccl5', 'Il2rb', 
                'Cd74', 'Igkc', 'Cd79a', 'Iglc2', 'Ms4a1', 
                'Cxcr2', 'Il1r2', 'Cd177', 'Mmp9', 'S100a9',
                'Gypa', 'Alas2', 'Hba-a1', 'Hba-a2']
adata.obs['cell_anno'] = pd.Categorical(adata.obs['cell_anno'], categories=['SMC', 'Fibroblast', 'EC', 'MoMaphDC', 'Schwann', 'T cell',
                                                                            'B cell', 'Neutrophil', 'Erythrocyte'])
sc.pl.matrixplot(adata, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='cell_anno_heatmap.pdf')

adata[adata.obs['wt_or_ex']=='wt'].obs['cell_anno'].value_counts()
# B cell         19611
# T cell         12165
# SMC             7146
# Fibroblast      3545
# MoMaphDC        1757
# EC              1318
# Schwann          365
# Erythrocyte      117
# Neutrophil        82

adata[adata.obs['wt_or_ex']=='ex'].obs['cell_anno'].value_counts()
# Neutrophil     13106
# MoMaphDC       10446
# Fibroblast      5498
# T cell          3508
# B cell          1978
# SMC             1022
# EC               632
# Schwann           29
# Erythrocyte        1

adata.write('res.h5ad')

## plot marker genes in UMAP
from matplotlib.colors import LinearSegmentedColormap

sc.pl.umap(adata, color=['Tagln', 'Dcn', 'Pecam1', 'Lyz2', 'Mpz', 'Il7r', 'Cd74', 'Cxcr2', 'Gypa'],
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
df['cell_anno'] = pd.Categorical(df['cell_anno'], categories=['SMC', 'Fibroblast', 'EC', 'MoMaphDC', 'Schwann', 'T cell',
                                                              'B cell', 'Neutrophil', 'Erythrocyte'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='wt_or_ex', y='proportion', fill='cell_anno')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') + coord_flip() +\
                                                                      scale_y_continuous(limits=[0, 1], breaks=np.arange(0, 1+0.1, 0.2)) + theme_bw()  
# scale_fill_brewer(type='qualitative',palette='Set1') +\
p.save(filename='cell_proportion_wt_ex.pdf', dpi=600, height=3, width=6)

## SMC & Fibroblast & EC
import scanpy as sc
import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

adata = sc.read_h5ad('res.h5ad')
wt_df = (adata[(adata.obs['wt_or_ex']=='wt') & (adata.obs['cell_anno'].isin(['SMC', 'Fibroblast', 'EC']))].obs['cell_anno'].value_counts(normalize=True)).reset_index()
wt_df['wt_or_ex'] = 'wt'
ex_df = (adata[(adata.obs['wt_or_ex']=='ex') & (adata.obs['cell_anno'].isin(['SMC', 'Fibroblast', 'EC']))].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ex_df['wt_or_ex'] = 'ex'
df = pd.concat([wt_df, ex_df])
df.reset_index(drop=True, inplace=True)
df['cell_anno'] = pd.Categorical(df['cell_anno'], categories=['SMC', 'Fibroblast', 'EC'])
df['wt_or_ex'] = pd.Categorical(df['wt_or_ex'], categories=['wt', 'ex'])

plt.rcParams['pdf.fonttype'] = 42
# without coord_flip
p = ggplot(df, aes(x='wt_or_ex', y='proportion', fill='cell_anno')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') +\
                                                                      scale_y_continuous(limits=[0, 1], breaks=np.arange(0, 1+0.1, 0.2)) + theme_bw()  
p.save(filename='cell_proportion_wt_ex_smc_fibro_ec.pdf', dpi=600, height=4, width=3)


#### Fibroblast
adata = sc.read_h5ad('res.h5ad')
fibro = adata[adata.obs['cell_anno']=='Fibroblast'].copy()   # 9043 × 2660

# rcParams["figure.figsize"] = (5, 5)
# sc.pl.umap(fibro, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1,
#            title='', frameon=True, save='_fibro_cell_anno.pdf')
# sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='wt'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1,
#            title='', frameon=True, save='_fibro_cell_anno_wt.pdf')
# sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='ex'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1,
#            title='', frameon=True, save='_fibro_cell_anno_ex.pdf')

sc.tl.leiden(fibro, resolution=0.1)
fibro.obs['leiden'].value_counts()
# 0    5810
# 1    2239
# 2     994

fibro[fibro.obs['wt_or_ex']=='wt'].obs['leiden'].value_counts()
# 1    2198
# 2     980
# 0     367

fibro[fibro.obs['wt_or_ex']=='ex'].obs['leiden'].value_counts()
# 0    5443
# 1      41
# 2      14

sc.pl.umap(fibro, color='leiden', legend_fontsize='small', legend_loc='on data', size=2,
           title='', frameon=True, save='_fibro_leiden.pdf')

sc.tl.rank_genes_groups(fibro, 'leiden')
pd.DataFrame(fibro.uns['rank_genes_groups']['names']).head(10)
#          0         1      2
# 0   Cthrc1       Gsn    Gsn
# 1    Postn      Fmo2   Igkc
# 2    Cxcl2     Gstm1    Dcn
# 3   Col8a1  Serping1   Ddx5
# 4    Thbs2       Dcn    Ubb
# 5     Lyz2     Itm2b   Cd74
# 6  Col11a1      Cst3  Rps27
# 7   Col5a1    Clec3b  Gstm1
# 8    Mmp14     Plpp3   Eif1
# 9   Col5a2    Rnase4   Junb

pd.DataFrame(fibro.uns['rank_genes_groups']['names']).head(100).to_csv('fibro_cluster_marker_genes.txt', index=False, sep='\t')

del fibro.uns['dendrogram_leiden']
sc.pl.rank_genes_groups_dotplot(fibro, groupby="leiden", n_genes=10, cmap='bwr', values_to_plot='logfoldchanges',
                                vmin=-4, vmax=4, save='fibro_marker_gene.pdf') # min_logfoldchange=2

fibro.obs['cell_anno'] = fibro.obs['leiden']
fibro.obs['cell_anno'].replace({'0':'fibro_1', '1':'fibro_2', '2':'fibro_3'}, inplace=True)  # to highlight cluster 1
fibro.obs['cell_anno'].value_counts()
# fibro_1    5810
# fibro_2    2239
# fibro_3     994

marker_genes = ['Cthrc1', 'Postn', 'Cxcl2', 'Col8a1', 'Thbs2',
                'Fmo2', 'Serping1', 'Cst3', 'Clec3b', 'Plpp3',
                'Igkc', 'Ddx5', 'Cd74', 'Rps27', 'Eif1']
fibro.obs['cell_anno'] = pd.Categorical(fibro.obs['cell_anno'], categories=['fibro_1', 'fibro_2', 'fibro_3'])
plt.rcParams['pdf.fonttype'] = 42
sc.pl.matrixplot(fibro, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='fibro_cell_anno_heatmap.pdf')

## wt & ex
rcParams["figure.figsize"] = (3, 3)
sc.pl.umap(fibro, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True)
plt.xlim(1, 12)
plt.ylim(-8, 3) 
plt.savefig('./figures/umap_fibro_cell_anno.pdf')
plt.close()

## wt
sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='wt'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True)
plt.xlim(1, 12)
plt.ylim(-8, 3) 
plt.savefig('./figures/umap_fibro_cell_anno_wt.pdf')
plt.close()

## ex
sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='ex'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True)
plt.xlim(1, 12)
plt.ylim(-8, 3) 
plt.savefig('./figures/umap_fibro_cell_anno_ex.pdf')
plt.close()


## plot umap for Cthrc1
from matplotlib.colors import LinearSegmentedColormap

sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='wt'], color=['Cthrc1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=10, vmin=0, vmax=5, title='', frameon=True)
plt.xlim(1, 12)
plt.ylim(-8, 3) 
plt.savefig('./figures/umap_fibro_cthrc1_wt.pdf')
plt.close()

sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='ex'], color=['Cthrc1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=10, vmin=0, vmax=5, title='', frameon=True)
plt.xlim(1, 12)
plt.ylim(-8, 3) 
plt.savefig('./figures/umap_fibro_cthrc1_ex.pdf')
plt.close()

fibro.write('res_fibro.h5ad')


#### 2 clusters for fibroblasts
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42

adata = sc.read_h5ad('res.h5ad')
fibro = adata[adata.obs['cell_anno']=='Fibroblast'].copy()   # 9043 × 2660

sc.tl.leiden(fibro, resolution=0.1)
fibro.obs['leiden'].value_counts()
# 0    5810
# 1    2239
# 2     994

fibro[fibro.obs['wt_or_ex']=='wt'].obs['leiden'].value_counts()
# 1    2198
# 2     980
# 0     367

fibro[fibro.obs['wt_or_ex']=='ex'].obs['leiden'].value_counts()
# 0    5443
# 1      41
# 2      14

# sc.pl.umap(fibro, color='leiden', legend_fontsize='small', legend_loc='on data', size=2,
#            title='', frameon=True, save='_fibro_leiden.pdf')

fibro.obs['cell_anno'] = fibro.obs['leiden']
fibro.obs['cell_anno'].replace({'0':'fibro_1', '1':'fibro_2', '2':'fibro_2'}, inplace=True)
fibro.obs['cell_anno'].value_counts()
# fibro_1    5810
# fibro_2    3233

sc.tl.rank_genes_groups(fibro, 'cell_anno')
pd.DataFrame(fibro.uns['rank_genes_groups']['names']).head(10)
# fibro_1  fibro_2
# 0   Cthrc1      Gsn
# 1    Postn    Gstm1
# 2    Cxcl2     Igkc
# 3   Col8a1      Dcn
# 4    Thbs2  Rarres2
# 5     Lyz2    Itm2b
# 6  Col11a1     Cd34
# 7   Col5a1     Fmo2
# 8    Mmp14      Cfh
# 9   Col5a2    Pmp22

pd.DataFrame(fibro.uns['rank_genes_groups']['names']).head(100).to_csv('fibro_cluster_marker_genes_2.txt', index=False, sep='\t')

marker_genes = ['Cthrc1', 'Postn', 'Cxcl2', 'Col8a1', 'Thbs2',
                'Gsn', 'Gstm1', 'Igkc', 'Dcn', 'Rarres2']
fibro.obs['cell_anno'] = pd.Categorical(fibro.obs['cell_anno'], categories=['fibro_1', 'fibro_2'])

sc.pl.matrixplot(fibro, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='fibro_cell_anno_heatmap_2.pdf')

## wt & ex
# rcParams["figure.figsize"] = (3, 3)
sc.pl.umap(fibro, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True)
plt.xlim(1, 12)
plt.ylim(-8, 3) 
plt.savefig('./figures/umap_fibro_cell_anno_2.pdf')
plt.close()

## wt
sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='wt'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True)
plt.xlim(1, 12)
plt.ylim(-8, 3) 
plt.savefig('./figures/umap_fibro_cell_anno_wt_2.pdf')
plt.close()

## ex
sc.pl.umap(fibro[fibro.obs['wt_or_ex']=='ex'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True)
plt.xlim(1, 12)
plt.ylim(-8, 3) 
plt.savefig('./figures/umap_fibro_cell_anno_ex_2.pdf')
plt.close()

fibro.write('res_fibro_2.h5ad')


#### Microarray dataset scatter plot
from plotnine import *
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_table('A_pos_A_neg.txt')
df['pos_neg'] = df['A_pos']/df['A_neg']

def assign_label(x):
    if x>(3/2):
        return '1'
    elif x<(2/3):
        return '-1'
    else:
        return '0'

df['up_down'] = df['pos_neg'].apply(assign_label)
df['up_down'] = pd.Categorical(df['up_down'], categories=['1', '-1'])

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='A_neg', y='A_pos', color='up_down')) + geom_point(size=0.5) +\
    scale_color_manual(values=['red', 'blue']) +\
    scale_x_continuous(limits=[0, 30000], breaks=range(0, 30000+1, 5000)) +\
    scale_y_continuous(limits=[0, 30000], breaks=range(0, 30000+1, 5000)) +\
    theme_bw()
p.save(filename='A_pos_A_neg.pdf', dpi=600, height=4, width=5)


#### GSE155468 (single cell)
#### public human single cell dataset
#### GSE155468
## https://explore.data.humancellatlas.org/projects/07073c12-8006-4710-a00b-23abdb814904/project-matrices

### extract fibroblasts
import pandas as pd

dat = pd.read_csv('AorticTissue_cell_types_annotations.csv')
dat['celltype'].value_counts()
# Tcell         18950
# MonoMaphDC    13475
# SMC1           4775
# Fibroblast     3652
# NK             2499
# SMC2           1666
# MSC             996
# EC              891
# Plasma          654
# Mastcell        265
# Bcell           259

fibro = dat[dat['celltype']=='Fibroblast'].iloc[:, [0,1,4]]
fibro.to_csv('fibroblast_info.txt', sep='\t', index=False, header=False)

## can not know sample info from loom files
## use raw umi files from GEO
import pandas as pd
import scanpy as sc
import numpy as np
import glob
import harmonypy as hm
import anndata as ad

def gen_adata(sample):
    raw = pd.read_csv(glob.glob('*'+sample+'*')[0], sep='\t')
    adata = sc.AnnData(raw.T.values, obs=raw.columns.values, var=raw.index.values)
    adata.obs.index = sample+'-'+adata.obs[0].values
    adata.obs['sample'] = sample
    
    if sample[:3]=='Con':
        adata.obs['sample_type'] = 'Control'
    else:
        adata.obs['sample_type'] = 'TAA'
    
    del adata.obs[0]
    adata.var.index = adata.var[0].values
    del adata.var[0]
    
    return adata

con_4 = gen_adata('Con4')
con_6 = gen_adata('Con6')
taa_1 = gen_adata('TAA1')
taa_2 = gen_adata('TAA2')
taa_3 = gen_adata('TAA3')
taa_4 = gen_adata('TAA4')
taa_5 = gen_adata('TAA5')
taa_6 = gen_adata('TAA6')
taa_7 = gen_adata('TAA7')
taa_8 = gen_adata('TAA8')

adata = ad.concat([con_4, con_6, taa_1, taa_2, taa_3, taa_4, taa_5, taa_6, taa_7, taa_8])  # 44221 × 12429
sc.pp.filter_cells(adata, min_genes=500)     # 41466 × 12429
sc.pp.filter_genes(adata, min_cells=10)      # 41466 × 12429

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt<10, :]  # 41466 × 12429
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='sample')  # set batch_key is important!!!
adata.raw = adata
adata = adata[:, adata.var.highly_variable]  # 41466 × 1448
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'], n_jobs=40)
# sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata)
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'sample')  # Converged after 3 iterations
adata.obsm['X_pca_harmony'] = ho.Z_corr.T
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)

sc.pl.umap(adata, color='sample_type', legend_fontsize='5', legend_loc='right margin', size=2, 
           title='', frameon=True, save='_batch_harmony.pdf')
sc.pl.umap(adata, color='sample_type', groups=['Control'], legend_fontsize='5', legend_loc='right margin', size=2, 
           title='', frameon=True, save='_batch_harmony_control.pdf')

fibro = pd.read_csv('fibroblast_info.txt', sep='\t', header=None)

def mod_idx(x):
    if x[:3]=='Con':
        return x.split('_')[0]
    else:
        return x.split('_')[0][1:]

fibro_idx = (fibro[1]+'-'+fibro[0]).apply(mod_idx)

adata_fibro = adata[[i for i in fibro_idx.values if i in adata.obs.index]].copy()  # 2910 × 1448
adata_fibro.obs['sample_type'].value_counts()
# TAA        2267
# Control     643

sc.pl.umap(adata_fibro, color='sample_type', legend_fontsize='5', legend_loc='right margin', size=5, 
           title='', frameon=True, save='_fibro.pdf')

## plot umap for CTHRC1
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sc.pl.umap(adata_fibro[adata_fibro.obs['sample_type']=='Control'], color=['CTHRC1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=30, vmin=0, vmax=5, title='', frameon=True)
plt.xlim(-5, 3)
plt.ylim(5, 16) 
plt.savefig('./figures/umap_fibro_cthrc1_control.pdf')
plt.close()

sc.pl.umap(adata_fibro[adata_fibro.obs['sample_type']=='TAA'], color=['CTHRC1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=30, vmin=0, vmax=5, title='', frameon=True)
plt.xlim(-5, 3)
plt.ylim(5, 16) 
plt.savefig('./figures/umap_fibro_cthrc1_taa.pdf')
plt.close()
