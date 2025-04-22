## workdir: /fs/home/jiluzhang/xiaoqing20231204/Science/Figure_2

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import harmonypy as hm
from matplotlib import rcParams

sc.settings.set_figure_params(dpi=600)
plt.rcParams['pdf.fonttype'] = 42

np.random.seed(0)

wt_1_raw = sc.read_10x_mtx('./raw_data/WT-T-80')   # 20785 × 53715
wt_1 = wt_1_raw[np.random.choice(wt_1_raw.n_obs, (wt_1_raw.n_obs)//2, replace=False), :].copy()  # 10392 × 53715

wt_2_raw = sc.read_10x_mtx('./raw_data/WT-T-82')   # 19151 × 53715
wt_2 = wt_2_raw[np.random.choice(wt_2_raw.n_obs, (wt_2_raw.n_obs)//2, replace=False), :].copy()  # 9575 × 53715

ko_1 = sc.read_10x_mtx('./raw_data/KO-T-76')   # 16290 × 53715
# ko_2 = sc.read_10x_mtx('./raw_data/KO-T-77') # discarded

wt_1.obs_names = ['wt-1_'+i for i in wt_1.obs_names]
wt_2.obs_names = ['wt-2_'+i for i in wt_2.obs_names]
ko_1.obs_names = ['ko-1_'+i for i in ko_1.obs_names]
# ko_2.obs_names = ['ko-2_'+i for i in ko_2.obs_names]

wt_1_gene_info = wt_1.var.copy()
wt_1_gene_info['gene_name'] = wt_1_gene_info.index.values
wt_1_gene_info.reset_index(drop=True, inplace=True)

wt_2_gene_info = wt_2.var.copy()
wt_2_gene_info['gene_name'] = wt_2_gene_info.index.values
wt_2_gene_info.reset_index(drop=True, inplace=True)

ko_1_gene_info = ko_1.var.copy()
ko_1_gene_info['gene_name'] = ko_1_gene_info.index.values
ko_1_gene_info.reset_index(drop=True, inplace=True)

# ko_2_gene_info = ko_2.var.copy()
# ko_2_gene_info['gene_name'] = ko_2_gene_info.index.values
# ko_2_gene_info.reset_index(drop=True, inplace=True)

gene_info = pd.concat([wt_1_gene_info, wt_2_gene_info, ko_1_gene_info])
gene_info.drop_duplicates(inplace=True)  # 53715*2
gene_info.index = gene_info['gene_name'].values

adata = ad.concat([wt_1, wt_2, ko_1], join='outer')   # 36257 × 53715
adata.var['gene_ids'] = gene_info.loc[adata.var.index]['gene_ids'].values

sc.pp.filter_cells(adata, min_genes=500)     # 28700 × 53715
sc.pp.filter_genes(adata, min_cells=10)      # 28700 × 21573

adata.obs['sample_groups'] = [i.split('_')[0] for i in adata.obs.index]
adata.obs['wt_or_ko'] = [i.split('-')[0] for i in adata.obs['sample_groups']]

adata.obs['wt_or_ko'].value_counts()
# wt    16541
# ko    12159

## qc for mt genes
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt<10, :]  # 27917 × 21573

## normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

## hvg identification
sc.pp.highly_variable_genes(adata, batch_key='wt_or_ko')
adata.raw = adata
adata = adata[:, adata.var.highly_variable]  # 27917 × 2346

## plot qc metrics
# n_genes_by_counts: the number of genes expressed in the count matrix
sc.pl.violin(adata, keys='n_genes_by_counts', groupby='wt_or_ko', rotation=90, stripplot=False,
             order=['wt', 'ko'], save='_n_genes.pdf')

# total_counts: the total counts per cell
sc.pl.violin(adata, keys='total_counts', groupby='wt_or_ko', rotation=90, stripplot=False,
             order=['wt', 'ko'], save='_nCounts.pdf')

# pct_counts_mt: the percentage of counts in mitochondrial genes
sc.pl.violin(adata, keys='pct_counts_mt', groupby='wt_or_ko', rotation=90, stripplot=False,
             order=['wt', 'ko'], save='_pct_mt.pdf')

# actb expression level
sc.pl.violin(adata, keys='Actb', groupby='wt_or_ko', rotation=90, stripplot=False,
             order=['wt', 'ko'], save='_actb_exp.pdf')

## PCA
sc.tl.pca(adata)

## Clustering
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'wt_or_ko')  # Converged after 2 iterations
adata.obsm['X_pca_harmony'] = ho.Z_corr.T
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
# sc.pl.umap(adata, color='sample_groups', legend_fontsize='5', legend_loc='right margin', size=2, 
#            title='', frameon=True, save='_batch_harmony_sample_groups.pdf')
sc.pl.umap(adata, color='wt_or_ko', legend_fontsize='5', legend_loc='right margin', size=2, 
           title='', frameon=True, save='_batch_harmony_wt_or_ko.pdf')

sc.tl.leiden(adata, resolution=0.5)
adata.obs['leiden'].value_counts()
# 0     5722
# 1     4802
# 2     4092
# 3     3392
# 4     2543
# 5     1962
# 6     1569
# 7     1146
# 8     1116
# 9      515
# 10     376
# 11     355
# 12     171
# 13     156

sc.pl.umap(adata, color='leiden', legend_fontsize='small', legend_loc='on data', size=1.5,
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
# KeyError: "Could not find keys '['Alas2', 'Bpifa1', 'Gypa', 'Hbb', 'Prmp', 'Reg3g']' in columns of `adata.obs` or in adata.raw.var_names."

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
                'Hba-a1', 'Hba-a2']

sc.pl.dotplot(adata, marker_genes, groupby='leiden', dendrogram=True, standard_scale='var',
              save='marker_genes.pdf')

sc.tl.rank_genes_groups(adata, 'leiden')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20).to_csv('cluster_marker_genes.txt', index=False, sep='\t')

cell_anno_dict = {'0':'Fibroblast', '1':'Fibroblast', '2':'MoMaphDC', '3':'Fibroblast', '4':'MoMaphDC', '5':'SMC',
                  '6':'T cell', '7':'Fibroblast', '8':'EC', '9':'Neutrophil', '10':'B cell', '11':'Schwann',
                  '12':'MoMaphDC', '13':'Fibroblast'}
adata.obs['cell_anno'] = adata.obs['leiden']
adata.obs['cell_anno'].replace(cell_anno_dict, inplace=True)

rcParams["figure.figsize"] = (2, 2)
sc.pl.umap(adata, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1.5,
           title='', frameon=True, save='_cell_anno.png')

marker_genes = ['Tagln', 'Myh11', 'Acta2', 'Myl9', 'Tpm2',
                'Dcn', 'Dpt', 'Col1a2', 'Igfbp6',
                'Pecam1', 'Aqp1', 'Vwf', 'Fabp4', 'Egfl7', 
                'Lyz2', 'C1qa', 'C1qb', 'C1qc',
                'Mpz', 'Mbp', 'Plp1', 'Ncmap', 'Kcna1',
                'Il7r', 'Cd3g', 'Ccl5', 'Il2rb', 
                'Igkc', 'Cd79a', 'Iglc2', 'Ms4a1', 
                'Cxcr2', 'Il1r2', 'Cd177', 'Mmp9', 'S100a9']
adata.obs['cell_anno'] = pd.Categorical(adata.obs['cell_anno'], categories=['SMC', 'Fibroblast', 'EC', 'MoMaphDC', 'Schwann', 'T cell',
                                                                            'B cell', 'Neutrophil'])
sc.pl.matrixplot(adata, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='cell_anno_heatmap.pdf')

adata[adata.obs['wt_or_ko']=='wt'].obs['cell_anno'].value_counts()
# Fibroblast    6223
# MoMaphDC      5622
# T cell        1464
# SMC            875
# EC             668
# Neutrophil     465
# B cell         320
# Schwann        300

adata[adata.obs['wt_or_ko']=='ko'].obs['cell_anno'].value_counts()
# Fibroblast    8995
# MoMaphDC      1184
# SMC           1087
# EC             448
# T cell         105
# B cell          56
# Schwann         55
# Neutrophil      50

adata.write('res.h5ad')

## plot marker genes in UMAP
from matplotlib.colors import LinearSegmentedColormap

sc.pl.umap(adata, color=['Tagln', 'Dcn', 'Pecam1', 'Lyz2', 'Mpz', 'Il7r', 'Igkc', 'Cxcr2'],
           cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_loc="on data", frameon=True, ncols=4, save='_marker_genes.pdf')

#### plot bar for cell proportion
wt_df = (adata[adata.obs['wt_or_ko']=='wt'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
wt_df['wt_or_ko'] = 'wt'
ko_df = (adata[adata.obs['wt_or_ko']=='ko'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ko_df['wt_or_ko'] = 'ko'
df = pd.concat([wt_df, ko_df])
df.reset_index(drop=True)
df.to_csv('cell_proportion_wt_ko.txt', index=False, sep='\t')

import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('cell_proportion_wt_ko.txt')
df['cell_anno'] = pd.Categorical(df['cell_anno'], categories=['SMC', 'Fibroblast', 'EC', 'MoMaphDC', 'Schwann', 'T cell',
                                                              'B cell', 'Neutrophil'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='wt_or_ko', y='proportion', fill='cell_anno')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') + coord_flip() +\
                                                                      scale_y_continuous(limits=[0, 1], breaks=np.arange(0, 1+0.1, 0.2)) + theme_bw()  
# scale_fill_brewer(type='qualitative',palette='Set1') +\
p.save(filename='cell_proportion_wt_ko.pdf', dpi=600, height=3, width=6)

## SMC & Fibroblast & EC
import scanpy as sc
import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

adata = sc.read_h5ad('res.h5ad')
wt_df = (adata[(adata.obs['wt_or_ko']=='wt') & (adata.obs['cell_anno'].isin(['SMC', 'Fibroblast', 'EC']))].obs['cell_anno'].value_counts(normalize=True)).reset_index()
wt_df['wt_or_ko'] = 'wt'
ko_df = (adata[(adata.obs['wt_or_ko']=='ko') & (adata.obs['cell_anno'].isin(['SMC', 'Fibroblast', 'EC']))].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ko_df['wt_or_ko'] = 'ko'
df = pd.concat([wt_df, ko_df])
df.reset_index(drop=True, inplace=True)
df['cell_anno'] = pd.Categorical(df['cell_anno'], categories=['SMC', 'Fibroblast', 'EC'])
df['wt_or_ko'] = pd.Categorical(df['wt_or_ko'], categories=['wt', 'ko'])

plt.rcParams['pdf.fonttype'] = 42
# without coord_flip
p = ggplot(df, aes(x='wt_or_ko', y='proportion', fill='cell_anno')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') +\
                                                                      scale_y_continuous(limits=[0, 1], breaks=np.arange(0, 1+0.1, 0.2)) + theme_bw()  
p.save(filename='cell_proportion_wt_ko_smc_fibro_ec.pdf', dpi=600, height=4, width=3)


## plot umap for wt & ko
# sc.pl.umap(adata[adata.obs['wt_or_ko']=='wt'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1.5,
#            title='', frameon=True, save='_cell_anno_wt.png')
# sc.pl.umap(adata[adata.obs['wt_or_ko']=='ko'], color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1.5,
#            title='', frameon=True, save='_cell_anno_ko.png')

## randomly show 10000 cells
adata.obs['cell_anno'] = pd.Categorical(adata.obs['cell_anno'], categories=['Fibroblast', 'MoMaphDC', 'SMC', 'T cell', 'EC', 'Neutrophil', 'B cell', 'Schwann'])
np.random.seed(0)
adata_wt = adata[adata.obs['wt_or_ko']=='wt']
sc.pl.umap(adata_wt[np.random.choice(adata_wt.n_obs, 10000, replace=False)],
           color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1.5,
           title='', frameon=True, save='_cell_anno_wt_10000.png')
adata_ko = adata[adata.obs['wt_or_ko']=='ko']
sc.pl.umap(adata_ko[np.random.choice(adata_ko.n_obs, 10000, replace=False)],
           color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1.5,
           title='', frameon=True, save='_cell_anno_ko_10000.png')


## plot cthrc1 exp in fibroblast
fibro = adata[adata.obs['cell_anno']=='Fibroblast'].copy()  # 15218 × 2346

sc.pl.umap(fibro[fibro.obs['wt_or_ko']=='wt'], color=['Cthrc1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=5, vmin=0, vmax=6, title='', frameon=True)
plt.xlim(2, 18)
plt.ylim(-5, 8) 
plt.savefig('./figures/umap_fibro_cthrc1_wt.pdf')
plt.close()

sc.pl.umap(fibro[fibro.obs['wt_or_ko']=='ko'], color=['Cthrc1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=5, vmin=0, vmax=6, title='', frameon=True)
plt.xlim(2, 18)
plt.ylim(-5, 8) 
plt.savefig('./figures/umap_fibro_cthrc1_ko.pdf')
plt.close()

sc.pl.umap(fibro[fibro.obs['wt_or_ko']=='wt'], color=['Cthrc1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=15, vmin=0, vmax=7, title='', frameon=True)
plt.xlim(2, 18)
plt.ylim(-5, 8) 
plt.savefig('./figures/umap_fibro_cthrc1_wt_2.pdf')
plt.close()

sc.pl.umap(fibro[fibro.obs['wt_or_ko']=='ko'], color=['Cthrc1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=15, vmin=0, vmax=7, title='', frameon=True)
plt.xlim(2, 18)
plt.ylim(-5, 8) 
plt.savefig('./figures/umap_fibro_cthrc1_ko_2.pdf')
plt.close()


## fibroblast subcluster & proportion
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

adata = sc.read_h5ad('res.h5ad')
fibro = adata[adata.obs['cell_anno']=='Fibroblast'].copy()  # 15218 × 2346

# fibro_ref = sc.read_h5ad('/fs/home/jiluzhang/xiaoqing20231204/Science/Figure_1_new/res_fibro.h5ad')
# inter_genes = np.intersect1d(fibro.var.index, fibro_ref.var.index)
# fibro_inter = fibro[:, inter_genes].copy()
# fibro_ref_inter = fibro_ref[:, inter_genes].copy()
# sc.pp.neighbors(fibro_ref_inter, use_rep='X_pca_harmony')
# sc.tl.umap(fibro_ref_inter)
# sc.tl.ingest(fibro_inter, fibro_ref_inter, obs="leiden")
# sc.pl.umap(fibro_inter, color='leiden', save='_fibro_leiden_map.pdf')
# sc.pl.umap(fibro_ref_inter, color='leiden', save='_fibro_ref_leiden_map.pdf')
# sc.pl.umap(fibro_ref, color='leiden', save='_fibro_ref_leiden.pdf')

sc.tl.leiden(fibro, resolution=0.20)
fibro.obs['leiden'].value_counts().reset_index(drop=True).sort_index()
# 0    5624
# 1    4981
# 2    4426
# 3     153
# 4      34

fibro[fibro.obs['wt_or_ko']=='wt'].obs['leiden'].value_counts(normalize=True).sort_index()
# 0    0.287000
# 1    0.509079
# 2    0.191869
# 3    0.007392
# 4    0.004660

fibro[fibro.obs['wt_or_ko']=='ko'].obs['leiden'].value_counts(normalize=True).sort_index()
# 0    0.426681
# 1    0.201556
# 2    0.359311
# 3    0.011895
# 4    0.000556

marker_genes = ['Cthrc1', 'Postn', 'Cxcl2', 'Col8a1', 'Thbs2',
                'Fmo2', 'Serping1', 'Cst3', 'Clec3b', 'Plpp3',
                'Igkc', 'Ddx5', 'Cd74', 'Rps27', 'Eif1']
sc.pl.matrixplot(fibro, marker_genes, groupby='leiden', cmap='coolwarm', standard_scale='var', save='fibro_leiden_heatmap.pdf')

# fibro_1_marker_genes = ['Cthrc1', 'Postn', 'Cxcl2', 'Col8a1', 'Thbs2']
# fibro_2_marker_genes = ['Fmo2', 'Serping1', 'Cst3', 'Clec3b', 'Plpp3']
# fibro_3_marker_genes = ['Igkc', 'Ddx5', 'Cd74', 'Rps27', 'Eif1']
# sc.tl.score_genes(fibro, gene_list=fibro_1_marker_genes, score_name='fibro_1_score')
# sc.tl.score_genes(fibro, gene_list=fibro_2_marker_genes, score_name='fibro_2_score')
# sc.tl.score_genes(fibro, gene_list=fibro_3_marker_genes, score_name='fibro_3_score')
# fibro.obs['cell_anno'] = fibro.obs[['fibro_1_score', 'fibro_2_score', 'fibro_3_score']].idxmax(axis=1).apply(lambda x: x[:7])

fibro.obs['cell_anno'] = fibro.obs['leiden']
fibro.obs['cell_anno'].replace({'0':'fibro_1', '1':'fibro_1', '2':'fibro_2', '3':'fibro_2', '4':'fibro_3'}, inplace=True)

sc.pl.umap(fibro, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True, save='_fibro_cell_anno.pdf')

fibro.write('res_fibro.h5ad')

## plot heatmap for marker genes
marker_genes = ['Cthrc1', 'Postn', 'Thbs2',
                'Fmo2', 'Cst3', 'Clec3b',
                'Igkc', 'Ddx5', 'Rps27']
fibro.obs['cell_anno'] = pd.Categorical(fibro.obs['cell_anno'], categories=['fibro_1', 'fibro_2', 'fibro_3'])
plt.rcParams['pdf.fonttype'] = 42
sc.pl.matrixplot(fibro, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='fibro_cell_anno_heatmap.pdf')


## plot bar for proportion
wt_df = (fibro[fibro.obs['wt_or_ko']=='wt'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
wt_df['wt_or_ko'] = 'wt'
ko_df = (fibro[fibro.obs['wt_or_ko']=='ko'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ko_df['wt_or_ko'] = 'ko'
df = pd.concat([wt_df, ko_df])
df.reset_index(drop=True)
df.to_csv('fibro_cell_proportion_wt_ko.txt', index=False, sep='\t')

import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('fibro_cell_proportion_wt_ko.txt')
df['cell_anno'] = pd.Categorical(df['cell_anno'], categories=['fibro_1', 'fibro_2', 'fibro_3'])

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='wt_or_ko', y='proportion', fill='cell_anno')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') + coord_flip() +\
                                                                      scale_y_continuous(limits=[0, 1], breaks=np.arange(0, 1+0.1, 0.2)) + theme_bw()  
# scale_fill_brewer(type='qualitative',palette='Set1') +\
p.save(filename='fibro_cell_proportion_wt_ko.pdf', dpi=600, height=3, width=6)


## plot umap for fibroblast subcluster
plt.rcParams['pdf.fonttype'] = 42

sc.pl.umap(fibro[fibro.obs['wt_or_ko']=='wt'], color=['cell_anno'],
           legend_fontsize='xx-small', legend_loc='right margin', size=5, title='', frameon=True)
plt.xlim(2, 18)
plt.ylim(-5, 8) 
plt.savefig('./figures/umap_fibro_subcluster_wt.pdf')
plt.close()

sc.pl.umap(fibro[fibro.obs['wt_or_ko']=='ko'], color=['cell_anno'],
           legend_fontsize='xx-small', legend_loc='right margin', size=5, title='', frameon=True)
plt.xlim(2, 18)
plt.ylim(-5, 8) 
plt.savefig('./figures/umap_fibro_subcluster_ko.pdf')
plt.close()

fibro.obs['wt_or_ko'].value_counts()
# wt    6223
# ko    8995


## plot cthrc1 exp in fibroblast-1
# from matplotlib.colors import LinearSegmentedColormap

# sc.pl.umap(fibro[(fibro.obs['wt_or_ko']=='wt') & (fibro.obs['cell_anno']=='fibro_1')], color=['Cthrc1'], 
#            cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
#            legend_fontsize='xx-small', legend_loc='right margin', size=40, vmin=0, vmax=7, title='', frameon=True)
# plt.xlim(4, 13)
# plt.ylim(-2, 8) 
# plt.savefig('./figures/umap_fibro_1_cthrc1_wt.pdf')
# plt.close()

# sc.pl.umap(fibro[(fibro.obs['wt_or_ko']=='ko') & (fibro.obs['cell_anno']=='fibro_1')], color=['Cthrc1'], 
#            cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
#            legend_fontsize='xx-small', legend_loc='right margin', size=40, vmin=0, vmax=7, title='', frameon=True)
# plt.xlim(4, 13)
# plt.ylim(-2, 8) 
# plt.savefig('./figures/umap_fibro_1_cthrc1_ko.pdf')
# plt.close()


## fibroblast subcluster & proportion (for 2 clusters)
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['pdf.fonttype'] = 42

adata = sc.read_h5ad('res.h5ad')
fibro = adata[adata.obs['cell_anno']=='Fibroblast'].copy()  # 15218 × 2346

sc.tl.leiden(fibro, resolution=0.20)
fibro.obs['leiden'].value_counts().reset_index(drop=True).sort_index()
# 0    5624
# 1    4981
# 2    4426
# 3     153
# 4      34

fibro[fibro.obs['wt_or_ko']=='wt'].obs['leiden'].value_counts(normalize=True).sort_index()
# 0    0.287000
# 1    0.509079
# 2    0.191869
# 3    0.007392
# 4    0.004660

fibro[fibro.obs['wt_or_ko']=='ko'].obs['leiden'].value_counts(normalize=True).sort_index()
# 0    0.426681
# 1    0.201556
# 2    0.359311
# 3    0.011895
# 4    0.000556

sc.pl.umap(fibro, color='leiden', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True, save='_fibro_leiden_2.pdf')

fibro.obs['cell_anno'] = fibro.obs['leiden']
fibro.obs['cell_anno'].replace({'0':'fibro_1', '1':'fibro_1', '2':'fibro_2', '3':'fibro_2', '4':'fibro_2'}, inplace=True)

sc.pl.umap(fibro, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True, save='_fibro_cell_anno_2.pdf')

fibro.write('res_fibro_2.h5ad')

## plot heatmap for marker genes
marker_genes = ['Cthrc1', 'Postn', 'Cxcl2', 'Col8a1', 'Thbs2',
                'Gsn', 'Gstm1', 'Dcn', 'Rarres2']
fibro.obs['cell_anno'] = pd.Categorical(fibro.obs['cell_anno'], categories=['fibro_1', 'fibro_2'])

sc.pl.matrixplot(fibro, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='fibro_cell_anno_heatmap_2.pdf')


## plot bar for proportion
wt_df = (fibro[fibro.obs['wt_or_ko']=='wt'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
wt_df['wt_or_ko'] = 'wt'
ko_df = (fibro[fibro.obs['wt_or_ko']=='ko'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ko_df['wt_or_ko'] = 'ko'
df = pd.concat([wt_df, ko_df])
df.reset_index(drop=True)
df.to_csv('fibro_cell_proportion_wt_ko_2.txt', index=False, sep='\t')

import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42

df = pd.read_table('fibro_cell_proportion_wt_ko_2.txt')
df['cell_anno'] = pd.Categorical(df['cell_anno'], categories=['fibro_1', 'fibro_2'])

p = ggplot(df, aes(x='wt_or_ko', y='proportion', fill='cell_anno')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') + coord_flip() +\
                                                                      scale_y_continuous(limits=[0, 1], breaks=np.arange(0, 1+0.1, 0.2)) + theme_bw()  
# scale_fill_brewer(type='qualitative',palette='Set1') +\
p.save(filename='fibro_cell_proportion_wt_ko_2.pdf', dpi=600, height=3, width=6)

## plot umap for fibroblast subcluster
sc.pl.umap(fibro[fibro.obs['wt_or_ko']=='wt'], color=['cell_anno'],
           legend_fontsize='xx-small', legend_loc='right margin', size=5, title='', frameon=True)
plt.xlim(2, 18)
plt.ylim(-5, 8) 
plt.savefig('./figures/umap_fibro_subcluster_wt_2.pdf')
plt.close()

sc.pl.umap(fibro[fibro.obs['wt_or_ko']=='ko'], color=['cell_anno'],
           legend_fontsize='xx-small', legend_loc='right margin', size=5, title='', frameon=True)
plt.xlim(2, 18)
plt.ylim(-5, 8) 
plt.savefig('./figures/umap_fibro_subcluster_ko_2.pdf')
plt.close()
