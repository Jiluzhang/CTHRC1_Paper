#### Antibody treatment 
## workdir: /fs/home/jiluzhang/xiaoqing20231204/Science/Figure_4/wt_ab

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import harmonypy as hm
from matplotlib import rcParams

plt.rcParams['pdf.fonttype'] = 42

ig_1 = sc.read_10x_mtx('./raw_data/Ig_1')   # 13653 × 32589
ig_2 = sc.read_10x_mtx('./raw_data/Ig_2')   # 11683 × 32589
ab_1 = sc.read_10x_mtx('./raw_data/Ab_1')   # 14065 × 32589
ab_2 = sc.read_10x_mtx('./raw_data/Ab_2')   # 13812 × 32589

ig_1.obs_names = ['ig-1_'+i for i in ig_1.obs_names]
ig_2.obs_names = ['ig-2_'+i for i in ig_2.obs_names]
ab_1.obs_names = ['ab-1_'+i for i in ab_1.obs_names]
ab_2.obs_names = ['ab-2_'+i for i in ab_2.obs_names]

ig_1_gene_info = ig_1.var.copy()
ig_1_gene_info['gene_name'] = ig_1_gene_info.index.values
ig_1_gene_info.reset_index(drop=True, inplace=True)

ig_2_gene_info = ig_2.var.copy()
ig_2_gene_info['gene_name'] = ig_2_gene_info.index.values
ig_2_gene_info.reset_index(drop=True, inplace=True)

ab_1_gene_info = ab_1.var.copy()
ab_1_gene_info['gene_name'] = ab_1_gene_info.index.values
ab_1_gene_info.reset_index(drop=True, inplace=True)

ab_2_gene_info = ab_2.var.copy()
ab_2_gene_info['gene_name'] = ab_2_gene_info.index.values
ab_2_gene_info.reset_index(drop=True, inplace=True)

gene_info = pd.concat([ig_1_gene_info, ig_2_gene_info, ab_1_gene_info, ab_2_gene_info])
gene_info.drop_duplicates(inplace=True)  # 32589*2
gene_info.index = gene_info['gene_name'].values

adata = ad.concat([ig_1, ig_2, ab_1, ab_2], join='outer')   # 53213 × 32589
adata.var['gene_ids'] = gene_info.loc[adata.var.index]['gene_ids'].values

sc.pp.filter_cells(adata, min_genes=500)     # 43370 × 32589
sc.pp.filter_genes(adata, min_cells=10)      # 43370 × 20081

adata.obs['sample_groups'] = [i.split('_')[0] for i in adata.obs.index]
adata.obs['ig_or_ab'] = [i.split('-')[0] for i in adata.obs['sample_groups']]

## qc for mt genes
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt<10, :]  # 42632 × 20081

## normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

## hvg identification
sc.pp.highly_variable_genes(adata, batch_key='sample_groups')
adata.raw = adata
adata = adata[:, adata.var.highly_variable]  # 42632 × 2466

## plot qc metrics
# n_genes_by_counts: the number of genes expressed in the count matrix
sc.pl.violin(adata, keys='n_genes_by_counts', groupby='sample_groups', rotation=90, stripplot=False,
             order=['ig-1', 'ig-2', 'ab-1', 'ab-2'], save='_n_genes.pdf')

# total_counts: the total counts per cell
sc.pl.violin(adata, keys='total_counts', groupby='sample_groups', rotation=90, stripplot=False,
             order=['ig-1', 'ig-2', 'ab-1', 'ab-2'], save='_nCounts.pdf')

# pct_counts_mt: the percentage of counts in mitochondrial genes
sc.pl.violin(adata, keys='pct_counts_mt', groupby='sample_groups', rotation=90, stripplot=False,
             order=['ig-1', 'ig-2', 'ab-1', 'ab-2'], save='_pct_mt.pdf')

# actb expression level
sc.pl.violin(adata, keys='Actb', groupby='sample_groups', rotation=90, stripplot=False,
             order=['ig-1', 'ig-2', 'ab-1', 'ab-2'], save='_actb_exp.pdf')

## PCA
sc.tl.pca(adata)

## Clustering
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'sample_groups')  # Converged after 3 iterations
adata.obsm['X_pca_harmony'] = ho.Z_corr.T
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.pl.umap(adata, color='sample_groups', legend_fontsize='5', legend_loc='right margin', size=2, 
           title='', frameon=True, save='_batch_harmony_sample_groups.pdf')
sc.pl.umap(adata, color='ig_or_ab', legend_fontsize='5', legend_loc='right margin', size=2, 
           title='', frameon=True, save='_batch_harmony_ig_or_ab.pdf')

sc.tl.leiden(adata, resolution=0.5)
adata.obs['leiden'].value_counts()
# 0     10605
# 1      6165
# 2      5686
# 3      5557
# 4      4242
# 5      1959
# 6      1931
# 7      1757
# 8      1626
# 9      1515
# 10      590
# 11      354
# 12      328
# 13      253
# 14       64

sc.pl.umap(adata, color='leiden', legend_fontsize='small', legend_loc='right margin',
           title='', frameon=True, show=False)
plt.gcf().set_size_inches(1.5, 1.5)
plt.savefig('./figures/umap_leiden.png', dpi=600)
plt.close()

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
# KeyError: "Could not find keys '['Bpifa1', 'Hbb', 'Igkc', 'Prmp', 'Reg3g']' in columns of `adata.obs` or in adata.raw.var_names."

marker_genes = ['Tagln', 'Myh11', 'Acta2', 'Myl9', 'Tpm2',
                'Dcn', 'Gsn', 'Dpt', 'Col1a2', 'Igfbp6',
                'Pecam1', 'Aqp1', 'Vwf', 'Fabp4', 'Egfl7', 
                'Lyz2', 'C1qa', 'C1qb', 'C1qc', 'H2-Aa',
                'Mpz', 'Mbp', 'Plp1', 'Ncmap', 'Kcna1',
                'Il7r', 'Cd3g', 'Cd52', 'Ccl5', 'Il2rb', 
                'Cd74', 'Cd79a', 'Iglc2', 'Ms4a1', 
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

cell_anno_dict = {'0':'SMC', '1':'Fibroblast', '2':'Fibroblast', '3':'MoMaphDC', '4':'Fibroblast', '5':'Neutrophil',
                  '6':'EC', '7':'EC', '8':'MoMaphDC', '9':'T cell', '10':'Schwann', '11':'SMC',
                  '12':'B cell', '13':'MoMaphDC', '14':'Fibroblast'}
adata.obs['cell_anno'] = adata.obs['leiden']
adata.obs['cell_anno'].replace(cell_anno_dict, inplace=True)

sc.pl.umap(adata, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin',
           title='', frameon=True, show=False)
plt.gcf().set_size_inches(1.5, 1.5)
plt.savefig('./figures/umap_cell_anno.png', dpi=600)
plt.close()

sc.pl.umap(adata, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin',
           title='', frameon=True, save='_cell_anno.pdf')


marker_genes = ['Tagln', 'Myh11', 'Acta2', 'Myl9', 'Tpm2',
                'Dcn', 'Dpt', 'Col1a2', 'Igfbp6',
                'Pecam1', 'Vwf', 'Fabp4', 'Egfl7', 
                'Lyz2', 'C1qa', 'C1qb', 'C1qc',
                'Mpz', 'Mbp', 'Plp1', 'Ncmap', 'Kcna1',
                'Il7r', 'Cd3g', 'Ccl5', 'Il2rb', 
                'Cd74', 'Cd79a', 'Iglc2', 'Ms4a1',
                'Cxcr2', 'Il1r2', 'Cd177', 'Mmp9', 'S100a9']
adata.obs['cell_anno'] = pd.Categorical(adata.obs['cell_anno'], categories=['SMC', 'Fibroblast', 'EC', 'MoMaphDC', 'Schwann', 'T cell',
                                                                            'B cell', 'Neutrophil'])
sc.pl.matrixplot(adata, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='cell_anno_heatmap.pdf')

adata[adata.obs['ig_or_ab']=='ig'].obs['cell_anno'].value_counts()
# MoMaphDC      5691
# Fibroblast    5444
# SMC           2819
# EC            2001
# Neutrophil    1858
# T cell        1366
# Schwann        539
# B cell         288

adata[adata.obs['ig_or_ab']=='ab'].obs['cell_anno'].value_counts()
# Fibroblast    10713
# SMC            8140
# MoMaphDC       1745
# EC             1687
# T cell          149
# Neutrophil      101
# Schwann          51
# B cell           40

adata.write('res.h5ad')


## plot marker genes in UMAP
from matplotlib.colors import LinearSegmentedColormap

sc.settings.set_figure_params(dpi=600)
sc.pl.umap(adata, color=['Tagln', 'Dcn', 'Pecam1', 'Lyz2', 'Mpz', 'Il7r', 'Cd74', 'Cxcr2'],
           cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_loc="on data", frameon=True, ncols=4, save='_marker_genes.pdf')


#### plot bar for cell proportion
ig_df = (adata[adata.obs['ig_or_ab']=='ig'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ig_df['ig_or_ab'] = 'ig'
ab_df = (adata[adata.obs['ig_or_ab']=='ab'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ab_df['ig_or_ab'] = 'ab'
df = pd.concat([ig_df, ab_df])
df.reset_index(drop=True)
df.to_csv('cell_proportion_ig_ab.txt', index=False, sep='\t')

import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('cell_proportion_ig_ab.txt')
df['cell_anno'] = pd.Categorical(df['cell_anno'], categories=['SMC', 'Fibroblast', 'EC', 'MoMaphDC', 'Schwann', 'T cell',
                                                              'B cell', 'Neutrophil'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='ig_or_ab', y='proportion', fill='cell_anno')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') + coord_flip() +\
                                                                      scale_y_continuous(limits=[0, 1], breaks=np.arange(0, 1+0.1, 0.2)) + theme_bw()  
# scale_fill_brewer(type='qualitative',palette='Set1') +\
p.save(filename='cell_proportion_ig_ab.pdf', dpi=600, height=3, width=6)


## plot umap for ig & ab
## randomly show 20000 cells
adata.obs['cell_anno'] = pd.Categorical(adata.obs['cell_anno'], categories=['SMC', 'Fibroblast', 'MoMaphDC', 'Neutrophil', 'EC', 'T cell', 'Schwann', 'B cell'])
np.random.seed(0)
adata_ig = adata[adata.obs['ig_or_ab']=='ig']
sc.pl.umap(adata_ig[np.random.choice(adata_ig.n_obs, 20000, replace=False)],
           color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1.5,
           title='', frameon=True, save='_cell_anno_ig_20000.png')
adata_ab = adata[adata.obs['ig_or_ab']=='ab']
sc.pl.umap(adata_ab[np.random.choice(adata_ab.n_obs, 20000, replace=False)],
           color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=1.5,
           title='', frameon=True, save='_cell_anno_ab_20000.png')


## plot cthrc1 exp in fibroblast
fibro = adata[adata.obs['cell_anno']=='Fibroblast'].copy()  # 16157 × 2466

sc.pl.umap(fibro[fibro.obs['ig_or_ab']=='ig'], color=['Cthrc1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=5, vmin=0, vmax=10, title='', frameon=True)   # 5444 × 2466
plt.xlim(-9, 12)
plt.ylim(0, 18) 
plt.savefig('./figures/umap_fibro_cthrc1_ig.pdf')
plt.close()

sc.pl.umap(fibro[fibro.obs['ig_or_ab']=='ab'], color=['Cthrc1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
           legend_fontsize='xx-small', legend_loc='right margin', size=5, vmin=0, vmax=10, title='', frameon=True)   # 10713 × 2466
plt.xlim(-9, 12)
plt.ylim(0, 18) 
plt.savefig('./figures/umap_fibro_cthrc1_ab.pdf')
plt.close()


#### plot smc proportion vesus others
import scanpy as sc
import pandas as pd
from plotnine import *
import matplotlib.pyplot as plt
import numpy as np

adata = sc.read_h5ad('./res.h5ad')

ig_cnt = sum(adata.obs['ig_or_ab']=='ig')   # 20006
ab_cnt = sum(adata.obs['ig_or_ab']=='ab')   # 22626

ig_cnt_smc = sum((adata.obs['ig_or_ab']=='ig') & (adata.obs['cell_anno']=='SMC'))   # 2819
ab_cnt_smc = sum((adata.obs['ig_or_ab']=='ab') & (adata.obs['cell_anno']=='SMC'))   # 8140

df = pd.DataFrame({'cls': ['SMC']*2+['Others']*2,
                   'trt': ['ig', 'ab']*2,
                   'val': 0})
df.loc[(df['cls']=='SMC') & (df['trt']=='ig'), 'val'] = ig_cnt_smc/ig_cnt
df.loc[(df['cls']=='SMC') & (df['trt']=='ab'), 'val'] = ab_cnt_smc/ab_cnt
df.loc[(df['cls']=='Others') & (df['trt']=='ig'), 'val'] = 1-ig_cnt_smc/ig_cnt
df.loc[(df['cls']=='Others') & (df['trt']=='ab'), 'val'] = 1-ab_cnt_smc/ab_cnt
df.to_csv('smc_vs_others_proportion.txt', index=False, sep='\t')

df = pd.read_table('smc_vs_others_proportion.txt')
df['cls'] = pd.Categorical(df['cls'], categories=['Others', 'SMC'])
df['trt'] = pd.Categorical(df['trt'], categories=['ig', 'ab'])

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='trt', y='val', fill='cls')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') +\
                                                    scale_y_continuous(limits=[0, 1.0], breaks=np.arange(0, 1.0+0.1, 0.2)) + theme_bw()
p.save(filename='smc_vs_others_proportion.pdf', dpi=600, height=3, width=3)


#### smc & fibro & ec
import scanpy as sc
import pandas as pd
from plotnine import *
import matplotlib.pyplot as plt
import numpy as np

adata = sc.read_h5ad('./res.h5ad')

ig_cnt = sum(adata.obs['ig_or_ab']=='ig')   # 20006
ab_cnt = sum(adata.obs['ig_or_ab']=='ab')   # 22626

ig_cnt_smc = sum((adata.obs['ig_or_ab']=='ig') & (adata.obs['cell_anno']=='SMC'))   # 2819
ab_cnt_smc = sum((adata.obs['ig_or_ab']=='ab') & (adata.obs['cell_anno']=='SMC'))   # 8140
ig_cnt_fibro = sum((adata.obs['ig_or_ab']=='ig') & (adata.obs['cell_anno']=='Fibroblast'))   # 5444
ab_cnt_fibro = sum((adata.obs['ig_or_ab']=='ab') & (adata.obs['cell_anno']=='Fibroblast'))   # 10713
ig_cnt_ec = sum((adata.obs['ig_or_ab']=='ig') & (adata.obs['cell_anno']=='EC'))   # 2001
ab_cnt_ec = sum((adata.obs['ig_or_ab']=='ab') & (adata.obs['cell_anno']=='EC'))   # 1687

df = pd.DataFrame({'cls': ['SMC']*2+['Fibroblast']*2+['EC']*2,
                   'trt': ['ig', 'ab']*3,
                   'val': 0})
df.loc[(df['cls']=='SMC') & (df['trt']=='ig'), 'val'] = ig_cnt_smc/ig_cnt
df.loc[(df['cls']=='SMC') & (df['trt']=='ab'), 'val'] = ab_cnt_smc/ab_cnt
df.loc[(df['cls']=='Fibroblast') & (df['trt']=='ig'), 'val'] = ig_cnt_fibro/ig_cnt
df.loc[(df['cls']=='Fibroblast') & (df['trt']=='ab'), 'val'] = ab_cnt_fibro/ab_cnt
df.loc[(df['cls']=='EC') & (df['trt']=='ig'), 'val'] = ig_cnt_ec/ig_cnt
df.loc[(df['cls']=='EC') & (df['trt']=='ab'), 'val'] = ab_cnt_ec/ab_cnt

df['cls'] = pd.Categorical(df['cls'], categories=['SMC', 'Fibroblast', 'EC'])
df['trt'] = pd.Categorical(df['trt'], categories=['ig', 'ab'])

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='cls', y='val', fill='trt')) + geom_bar(stat='identity', position='dodge', width=0.75) + xlab('') + ylab('Fraction') +\
                                                    scale_y_continuous(limits=[0, 0.5], breaks=np.arange(0, 0.5+0.1, 0.1)) + theme_bw()
p.save(filename='smc_fibro_ec_proportion.pdf', dpi=600, height=3, width=3)

#### smc subclusters
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
sc.settings.set_figure_params(dpi=600)

adata = sc.read_h5ad('res.h5ad')
smc = adata[adata.obs['cell_anno']=='SMC'].copy()  # 10959 × 2466

sc.tl.leiden(smc, resolution=0.8)
sc.pl.umap(smc, color='leiden', legend_fontsize='medium', legend_loc='right margin', size=20, 
           title='', frameon=True, save='_smc_leiden.pdf')

smc.obs['leiden'].value_counts()
# 0    1847
# 1    1587
# 2    1559
# 3    1533
# 4    1447
# 5    1267
# 6    1119
# 7     468
# 8     132

marker_genes = ['Acta2', 'Actg1', 'Tagln2', 'Myl12a', 'Tubb2a',
                'Ccnl2', 'Ccnt2', 'Malat1', 'Ccn1', 'mt-Rnr1',
                'Col1a2', 'Lum', 'Dcn', 'Gsn', 'Cygb',
                'Notch3', 'Atf3', 'Rgs4', 'Cacna1h', 'Pln']
# "Could not find keys '['mt-Rnr1']' in columns of `adata.obs` or in adata.raw.var_names."

marker_genes = ['Acta2', 'Actg1', 'Tagln2', 'Myl12a', 'Tubb2a',
                'Ccnl2', 'Ccnt2', 'Malat1', 'Ccn1',
                'Col1a2', 'Lum', 'Dcn', 'Gsn', 'Cygb',
                'Notch3', 'Atf3', 'Rgs4', 'Cacna1h', 'Pln']

del smc.uns['dendrogram_leiden']
sc.pl.dotplot(smc, marker_genes, groupby='leiden', standard_scale='var', dendrogram=True,
              save='smc_marker_genes.pdf')

cell_anno_dict = {'0':'Proliferating SMC', '1':'Fibromyocyte', '2':'Stressed SMC',
                  '3':'Proliferating SMC', '4':'Contractile SMC', '5':'Contractile SMC',
                  '6':'Fibromyocyte', '7':'Stressed SMC', '8':'Fibromyocyte'}
smc.obs['cell_anno'] = smc.obs['leiden']
smc.obs['cell_anno'].replace(cell_anno_dict, inplace=True)

sc.pl.umap(smc, color='cell_anno', legend_fontsize='5', legend_loc='right margin', size=10,
           title='', frameon=True, show=False)
plt.gcf().set_size_inches(2, 2)
plt.savefig('./figures/umap_smc_cell_anno.pdf', dpi=600)
plt.close()

smc.obs['cell_anno'] = pd.Categorical(smc.obs['cell_anno'],
                                      categories=['Contractile SMC', 'Proliferating SMC', 'Fibromyocyte', 'Stressed SMC'])
marker_genes = ['Actg1', 'Tagln2', 'Myl12a',
                'Ccnl2', 'Malat1',
                'Lum', 'Dcn',
                'Notch3', 'Atf3', 'Cacna1h', 'Pln']
sc.pl.matrixplot(smc, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='smc_cell_anno_heatmap.pdf')

## output marker genes for each cluster
sc.tl.rank_genes_groups(smc, 'leiden')
pd.DataFrame(smc.uns['rank_genes_groups']['names']).head(20).to_csv('smc_cluster_marker_genes.txt', index=False, sep='\t')

smc.write('./smc_res.h5ad')

## smc proportion
ig_smc_cnt = sum(smc.obs['ig_or_ab']=='ig')   # 2819
ab_smc_cnt = sum(smc.obs['ig_or_ab']=='ab')   # 8140

ig_smc_con_cnt = sum((smc.obs['cell_anno']=='Contractile SMC') & (smc.obs['ig_or_ab']=='ig'))     # 612
ab_smc_con_cnt = sum((smc.obs['cell_anno']=='Contractile SMC') & (smc.obs['ig_or_ab']=='ab'))     # 2102

ig_smc_pro_cnt = sum((smc.obs['cell_anno']=='Proliferating SMC') & (smc.obs['ig_or_ab']=='ig'))   # 771
ab_smc_pro_cnt = sum((smc.obs['cell_anno']=='Proliferating SMC') & (smc.obs['ig_or_ab']=='ab'))   # 2609

ig_smc_str_cnt = sum((smc.obs['cell_anno']=='Stressed SMC') & (smc.obs['ig_or_ab']=='ig'))   # 699
ab_smc_str_cnt = sum((smc.obs['cell_anno']=='Stressed SMC') & (smc.obs['ig_or_ab']=='ab'))   # 1328

ig_smc_fib_cnt = sum((smc.obs['cell_anno']=='Fibromyocyte') & (smc.obs['ig_or_ab']=='ig'))   # 737
ab_smc_fib_cnt = sum((smc.obs['cell_anno']=='Fibromyocyte') & (smc.obs['ig_or_ab']=='ab'))   # 2101

smc_df = pd.DataFrame({'smc': ['Contractile SMC']*2+['Proliferating SMC']*2+['Stressed SMC']*2+['Fibromyocyte']*2,
                       'trt': ['ig', 'ab']*4,
                       'val': 0})
smc_df.loc[(smc_df['smc']=='Contractile SMC') & (smc_df['trt']=='ig'), 'val'] = ig_smc_con_cnt/ig_smc_cnt
smc_df.loc[(smc_df['smc']=='Contractile SMC') & (smc_df['trt']=='ab'), 'val'] = ab_smc_con_cnt/ab_smc_cnt
smc_df.loc[(smc_df['smc']=='Proliferating SMC') & (smc_df['trt']=='ig'), 'val'] = ig_smc_pro_cnt/ig_smc_cnt
smc_df.loc[(smc_df['smc']=='Proliferating SMC') & (smc_df['trt']=='ab'), 'val'] = ab_smc_pro_cnt/ab_smc_cnt
smc_df.loc[(smc_df['smc']=='Stressed SMC') & (smc_df['trt']=='ig'), 'val'] = ig_smc_str_cnt/ig_smc_cnt
smc_df.loc[(smc_df['smc']=='Stressed SMC') & (smc_df['trt']=='ab'), 'val'] = ab_smc_str_cnt/ab_smc_cnt
smc_df.loc[(smc_df['smc']=='Fibromyocyte') & (smc_df['trt']=='ig'), 'val'] = ig_smc_fib_cnt/ig_smc_cnt
smc_df.loc[(smc_df['smc']=='Fibromyocyte') & (smc_df['trt']=='ab'), 'val'] = ab_smc_fib_cnt/ab_smc_cnt
smc_df.to_csv('smc_subcluster_proportion.txt', index=False, sep='\t')

smc.obs['cell_anno'].value_counts()
# Proliferating SMC    3380
# Fibromyocyte         2838
# Contractile SMC      2714
# Stressed SMC         2027

#### plot bar for smc subcluster proportion
import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('smc_subcluster_proportion.txt')
smc = ['Contractile SMC', 'Proliferating SMC', 'Stressed SMC', 'Fibromyocyte']
df['smc'] = pd.Categorical(df['smc'], categories=smc)
df['trt'] = pd.Categorical(df['trt'], categories=['ig', 'ab'])

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='smc', y='val', fill='trt')) + geom_bar(stat='identity', position='dodge', width=0.75) + xlab('') + ylab('Fraction') +\
                                                    scale_y_continuous(limits=[0, 0.4], breaks=np.arange(0, 0.4+0.1, 0.1)) + theme_bw()
p.save(filename='smc_subcluster_proportion.pdf', dpi=600, height=3, width=6)


#### find differential genes
## Contractile SMC
smc_con = smc[smc.obs['cell_anno']=='Contractile SMC'].copy()  # 2714 × 2466
sc.tl.rank_genes_groups(smc_con, 'ig_or_ab')
result = smc_con.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_con_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_con_up_cnt = sum((smc_con_res['ab_logfoldchanges']>1) & (smc_con_res['ab_pvals_adj']<0.05))   # 253
smc_con_dw_cnt = sum((smc_con_res['ab_logfoldchanges']<-1) & (smc_con_res['ab_pvals_adj']<0.05))  # 409
smc_con_df_cnt = smc_con_up_cnt + smc_con_dw_cnt                                                  # 662
smc_con_res[(smc_con_res['ab_logfoldchanges']<-1) & (smc_con_res['ab_pvals_adj']<0.05)]['ab_names'].to_csv('smc_con_down_genes.txt', index=False, header=False)
smc_con_res[(smc_con_res['ab_logfoldchanges']>1) & (smc_con_res['ab_pvals_adj']<0.05)]['ab_names'].to_csv('smc_con_up_genes.txt', index=False, header=False)

## Proliferating SMC
smc_pro = smc[smc.obs['cell_anno']=='Proliferating SMC'].copy()  # 3380 × 2466
sc.tl.rank_genes_groups(smc_pro, 'ig_or_ab')
result = smc_pro.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_pro_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_pro_up_cnt = sum((smc_pro_res['ab_logfoldchanges']>1) & (smc_pro_res['ab_pvals_adj']<0.05))   # 312
smc_pro_dw_cnt = sum((smc_pro_res['ab_logfoldchanges']<-1) & (smc_pro_res['ab_pvals_adj']<0.05))  # 475
smc_pro_df_cnt = smc_pro_up_cnt + smc_pro_dw_cnt                                                  # 787
smc_pro_res[(smc_pro_res['ab_logfoldchanges']<-1) & (smc_pro_res['ab_pvals_adj']<0.05)]['ab_names'].to_csv('smc_pro_down_genes.txt', index=False, header=False)
smc_pro_res[(smc_pro_res['ab_logfoldchanges']>1) & (smc_pro_res['ab_pvals_adj']<0.05)]['ab_names'].to_csv('smc_pro_up_genes.txt', index=False, header=False)

## Stressed SMC
smc_str = smc[smc.obs['cell_anno']=='Stressed SMC'].copy()  # 2027 × 2466
sc.tl.rank_genes_groups(smc_str, 'ig_or_ab')
result = smc_str.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_str_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_str_up_cnt = sum((smc_str_res['ab_logfoldchanges']>1) & (smc_str_res['ab_pvals_adj']<0.05))   # 161
smc_str_dw_cnt = sum((smc_str_res['ab_logfoldchanges']<-1) & (smc_str_res['ab_pvals_adj']<0.05))  # 828
smc_str_df_cnt = smc_str_up_cnt + smc_str_dw_cnt                                                  # 989
smc_str_res[(smc_str_res['ab_logfoldchanges']<-1) & (smc_str_res['ab_pvals_adj']<0.05)]['ab_names'].to_csv('smc_str_down_genes.txt', index=False, header=False)
smc_str_res[(smc_str_res['ab_logfoldchanges']>1) & (smc_str_res['ab_pvals_adj']<0.05)]['ab_names'].to_csv('smc_str_up_genes.txt', index=False, header=False)

## Fibromyocyte
smc_fib = smc[smc.obs['cell_anno']=='Fibromyocyte'].copy()  # 2838 × 2466
sc.tl.rank_genes_groups(smc_fib, 'ig_or_ab')
result = smc_fib.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_fib_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_fib_up_cnt = sum((smc_fib_res['ab_logfoldchanges']>1) & (smc_fib_res['ab_pvals_adj']<0.05))   # 252
smc_fib_dw_cnt = sum((smc_fib_res['ab_logfoldchanges']<-1) & (smc_fib_res['ab_pvals_adj']<0.05))  # 553
smc_fib_df_cnt = smc_fib_up_cnt + smc_fib_dw_cnt                                                  # 805
smc_fib_res[(smc_fib_res['ab_logfoldchanges']<-1) & (smc_fib_res['ab_pvals_adj']<0.05)]['ab_names'].to_csv('smc_fib_down_genes.txt', index=False, header=False)
smc_fib_res[(smc_fib_res['ab_logfoldchanges']>1) & (smc_fib_res['ab_pvals_adj']<0.05)]['ab_names'].to_csv('smc_fib_up_genes.txt', index=False, header=False)

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
                                                    scale_y_continuous(limits=[0, 900], breaks=np.arange(0, 900+1, 300)) + theme_bw()
p.save(filename='smc_subcluster_deg_cnt.pdf', dpi=600, height=3, width=5)


#### GO enrichment analysis for DEGs
library(clusterProfiler)
library(org.Mm.eg.db)

## smc_con_up
smc_con_up <- read.table('smc_con_up_genes.txt')[, 1]
smc_con_up_ego <- enrichGO(gene=smc_con_up, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_con_up_genes_go.pdf')
dotplot(smc_con_up_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_con_up_ego), file="smc_con_up_genes_go.csv")

## smc_con_down
smc_con_down <- read.table('smc_con_down_genes.txt')[, 1]
smc_con_down_ego <- enrichGO(gene=smc_con_down, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_con_down_genes_go.pdf')
dotplot(smc_con_down_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_con_down_ego), file="smc_con_down_genes_go.csv")

## smc_pro_up
smc_pro_up <- read.table('smc_pro_up_genes.txt')[, 1]
smc_pro_up_ego <- enrichGO(gene=smc_pro_up, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_pro_up_genes_go.pdf')
dotplot(smc_pro_up_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_pro_up_ego), file="smc_pro_up_genes_go.csv")

## smc_pro_down
smc_pro_down <- read.table('smc_pro_down_genes.txt')[, 1]
smc_pro_down_ego <- enrichGO(gene=smc_pro_down, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_pro_down_genes_go.pdf')
dotplot(smc_pro_down_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_pro_down_ego), file="smc_pro_down_genes_go.csv")

## smc_str_up
smc_str_up <- read.table('smc_str_up_genes.txt')[, 1]
smc_str_up_ego <- enrichGO(gene=smc_str_up, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_str_up_genes_go.pdf')
dotplot(smc_str_up_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_str_up_ego), file="smc_str_up_genes_go.csv")

## smc_str_down
smc_str_down <- read.table('smc_str_down_genes.txt')[, 1]
smc_str_down_ego <- enrichGO(gene=smc_str_down, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_str_down_genes_go.pdf')
dotplot(smc_str_down_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_str_down_ego), file="smc_str_down_genes_go.csv")

## smc_fib_up
smc_fib_up <- read.table('smc_fib_up_genes.txt')[, 1]
smc_fib_up_ego <- enrichGO(gene=smc_fib_up, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_fib_up_genes_go.pdf')
dotplot(smc_fib_up_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_fib_up_ego), file="smc_fib_up_genes_go.csv")

## smc_fib_down
smc_fib_down <- read.table('smc_fib_down_genes.txt')[, 1]
smc_fib_down_ego <- enrichGO(gene=smc_fib_down, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_fib_down_genes_go.pdf')
dotplot(smc_fib_down_ego, showCategory=10, font.size=10)
dev.off()
write.csv(as.data.frame(smc_fib_down_ego), file="smc_fib_down_genes_go.csv")

terms <- c('GO:0032103', 'GO:0043254', 'GO:0045862', 'GO:0031589', 'GO:0045860',
           'GO:0001503', 'GO:0007160', 'GO:0014909', 'GO:0030198', 'GO:0051403')
pdf('smc_fib_down_genes_go_selected.pdf')
dotplot(smc_fib_down_ego, showCategory=smc_fib_down_ego@result$Description[smc_fib_down_ego@result$ID %in% terms], font.size=10)
dev.off()


#### plot box for specific geneset
from plotnine import *
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

smc = sc.read_h5ad('smc_res.h5ad')  # 10959 × 2466

## contractile genes
con_genes = ['Ramp1', 'Cnn1', 'Actg2', 'Tagln', 'Acta2', 'Tpm2', 'Myl6', 'Myl9', 'Myh11', 'Mylk']  # delete Mylk2 due to low expression
dat = pd.DataFrame()
for gene in con_genes:
    dat_ig = pd.DataFrame({'gene_id':gene,
                           'idx':'ig',
                           'val':smc[smc.obs['ig_or_ab']=='ig'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat_ab = pd.DataFrame({'gene_id':gene,
                           'idx':'ab',
                           'val':smc[smc.obs['ig_or_ab']=='ab'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat = pd.concat([dat, dat_ig, dat_ab])

dat['idx'] = pd.Categorical(dat['idx'], categories=['ig', 'ab'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(dat, aes(x='gene_id', y='val', fill='idx')) + geom_boxplot(width=0.5, outlier_shape='', show_legend=False) +\
                                                         scale_y_continuous(limits=[0, 7], breaks=np.arange(0, 7+0.1, 1.0)) + theme_bw()
p.save(filename='smc_contraction_seperately.pdf', dpi=600, height=4, width=5)

for gene in con_genes:
    print(gene, 'pvalue:', stats.ttest_ind(dat[(dat['idx']=='ig') & (dat['gene_id']==gene)]['val'],
                                           dat[(dat['idx']=='ab') & (dat['gene_id']==gene)]['val'])[1])
# Ramp1 pvalue: 1.1904422611762816e-186
# Cnn1 pvalue: 6.938803871008987e-163
# Actg2 pvalue: 5.2397676936744384e-08
# Tagln pvalue: 4.998066175766804e-118
# Acta2 pvalue: 2.0620044115931753e-59
# Tpm2 pvalue: 2.426914441161974e-109
# Myl6 pvalue: 5.3287015261517904e-213
# Myl9 pvalue: 3.5357452996145175e-302
# Myh11 pvalue: 1.477387477690079e-20
# Mylk pvalue: 0.1821665710202039

## mmp genes
mmp_genes = ['Mmp14', 'Mmp2']  # delete Mmp24 & Mmp8 & Mmp9 & Mmp15 & Mmp13 due to low expression
dat = pd.DataFrame()
for gene in mmp_genes:
    dat_ig = pd.DataFrame({'gene_id':gene,
                           'idx':'ig',
                           'val':smc[smc.obs['ig_or_ab']=='ig'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat_ab = pd.DataFrame({'gene_id':gene,
                           'idx':'ab',
                           'val':smc[smc.obs['ig_or_ab']=='ab'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat = pd.concat([dat, dat_ig, dat_ab])

dat['idx'] = pd.Categorical(dat['idx'], categories=['ig', 'ab'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(dat, aes(x='gene_id', y='val', fill='idx')) + geom_boxplot(width=0.5, outlier_shape='', show_legend=False) +\
                                                         scale_y_continuous(limits=[0, 4], breaks=np.arange(0, 4+0.1, 1)) + theme_bw()
p.save(filename='smc_mmp_seperately.pdf', dpi=600, height=4, width=5)

for gene in mmp_genes:
    print(gene, 'pvalue:', stats.ttest_ind(dat[(dat['idx']=='ig') & (dat['gene_id']==gene)]['val'],
                                           dat[(dat['idx']=='ab') & (dat['gene_id']==gene)]['val'])[1])
# Mmp14 pvalue: 5.231537047848063e-12
# Mmp2 pvalue: 7.142945897893802e-10


## fibroblast subcluster & proportion
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

adata = sc.read_h5ad('res.h5ad')
fibro = adata[adata.obs['cell_anno']=='Fibroblast'].copy()  # 16157 × 2466

sc.tl.leiden(fibro, resolution=0.20)
fibro.obs['leiden'].value_counts()
# 0    5963
# 1    5841
# 2    4289
# 3      64

fibro[fibro.obs['ig_or_ab']=='ig'].obs['leiden'].value_counts(normalize=True).sort_index()
# 0    0.298677
# 1    0.313740
# 2    0.376010
# 3    0.011572

fibro[fibro.obs['ig_or_ab']=='ab'].obs['leiden'].value_counts(normalize=True).sort_index()
# 0    0.404835
# 1    0.385793
# 2    0.209278
# 3    0.000093

marker_genes = ['Cthrc1', 'Postn', 'Cxcl2', 'Col8a1', 'Thbs2',
                'Fmo2', 'Serping1', 'Cst3', 'Clec3b', 'Plpp3',
                'Ddx5', 'Cd74', 'Rps27', 'Eif1']
sc.pl.matrixplot(fibro, marker_genes, groupby='leiden', cmap='coolwarm', standard_scale='var', save='fibro_leiden_heatmap.pdf')

fibro.obs['cell_anno'] = fibro.obs['leiden']
fibro.obs['cell_anno'].replace({'0':'fibro_2', '1':'fibro_1', '2':'fibro_1', '3':'fibro_3'}, inplace=True)

sc.pl.umap(fibro, color='cell_anno', legend_fontsize='xx-small', legend_loc='right margin', size=10,
           title='', frameon=True, save='_fibro_cell_anno.pdf')

fibro.write('res_fibro.h5ad')

## plot heatmap for marker genes
marker_genes = ['Cthrc1', 'Postn', 'Thbs2',
                'Fmo2', 'Cst3', 'Clec3b',
                'Igkc', 'Ddx5', 'Rps27']

marker_genes = ['Cthrc1', 'Postn', 'Col8a1', 'Thbs2',
                'Fmo2', 'Cst3',
                'Ddx5', 'Cd74', 'Rps27']
fibro.obs['cell_anno'] = pd.Categorical(fibro.obs['cell_anno'], categories=['fibro_1', 'fibro_2', 'fibro_3'])
plt.rcParams['pdf.fonttype'] = 42
sc.pl.matrixplot(fibro, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='fibro_cell_anno_heatmap.pdf')


## plot bar for proportion
ig_df = (fibro[fibro.obs['ig_or_ab']=='ig'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ig_df['ig_or_ab'] = 'ig'
ab_df = (fibro[fibro.obs['ig_or_ab']=='ab'].obs['cell_anno'].value_counts(normalize=True)).reset_index()
ab_df['ig_or_ab'] = 'ab'
df = pd.concat([ig_df, ab_df])
df.reset_index(drop=True)
df.to_csv('fibro_cell_proportion_ig_ab.txt', index=False, sep='\t')

import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('fibro_cell_proportion_ig_ab.txt')
df['cell_anno'] = pd.Categorical(df['cell_anno'], categories=['fibro_1', 'fibro_2', 'fibro_3'])

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='ig_or_ab', y='proportion', fill='cell_anno')) + geom_bar(stat='identity', width=0.75) + xlab('') + ylab('Fraction') + coord_flip() +\
                                                                      scale_y_continuous(limits=[0, 1], breaks=np.arange(0, 1+0.1, 0.2)) + theme_bw()  
# scale_fill_brewer(type='qualitative',palette='Set1') +\
p.save(filename='fibro_cell_proportion_ig_ab.pdf', dpi=600, height=3, width=6)

## plot umap for fibroblast subcluster
plt.rcParams['pdf.fonttype'] = 42

sc.pl.umap(fibro[fibro.obs['ig_or_ab']=='ig'], color=['cell_anno'],
           legend_fontsize='xx-small', legend_loc='right margin', size=5, title='', frameon=True)
plt.xlim(-7, 12)
plt.ylim(-2, 18) 
plt.savefig('./figures/umap_fibro_subcluster_ig.pdf')
plt.close()

sc.pl.umap(fibro[fibro.obs['ig_or_ab']=='ab'], color=['cell_anno'],
           legend_fontsize='xx-small', legend_loc='right margin', size=5, title='', frameon=True)
plt.xlim(-7, 12)
plt.ylim(-2, 18) 
plt.savefig('./figures/umap_fibro_subcluster_ab.pdf')
plt.close()

fibro.obs['ig_or_ab'].value_counts()
# ig     5444
# ab    10713
