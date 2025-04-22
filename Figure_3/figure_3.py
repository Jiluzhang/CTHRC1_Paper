#### cell-cell interaction
## https://github.com/jinworks/CellChat

# install.packages('devtools')

# Warning messages:
# 1: In install.packages("devtools") :
#   installation of package ‘textshaping’ had non-zero exit status
# 2: In install.packages("devtools") :
#   installation of package ‘ragg’ had non-zero exit status
# 3: In install.packages("devtools") :
#   installation of package ‘pkgdown’ had non-zero exit status
# 4: In install.packages("devtools") :
#   installation of package ‘devtools’ had non-zero exit status

# conda install conda-forge::r-devtools  # it works!!!

install.packages('ggsci')
install.packages('ggsignif')
install.packages('ggalluvial')
BiocManager::install("ComplexHeatmap")
install.packages('ggpubr')
BiocManager::install("BiocNeighbors")
install.packages('ggrepel')
install.packages('htmltools')
install.packages('anndata')
devtools::install_github("jinworks/CellChat")
devtools::install_github('immunogenomics/presto')


## Comparison analysis of multiple datasets using CellChat
## https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
## Inference and analysis of cell-cell communication using CellChat
## https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
## Update CellChatDB by integrating new ligand-receptor pairs from other resources or utilizing a custom ligand-receptor interaction database
## https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Update-CellChatDB.html

## workdir: /fs/home/jiluzhang/xiaoqing20231204/Science/Figure_4/cci
## cp /fs/home/jiluzhang/xiaoqing20231204/Science/Figure_1_new/res.h5ad .

# import scanpy as sc
# dat = sc.read_h5ad('res.h5ad')  # 82326 × 2660
# wt = dat[(dat.obs['sample_groups'].isin(['wt-1', 'wt-2']))&(dat.obs['cell_anno']!='Erythrocyte'), :].copy()  # Erythrocyte count is few
# sc.AnnData(wt.raw.X, obs=wt.obs, var=wt.raw.var).write('wt.h5ad')  # 45989 × 26020
# ex = dat[(dat.obs['sample_groups'].isin(['ex-1', 'ex-2']))&(dat.obs['cell_anno']!='Erythrocyte'), :].copy()
# sc.AnnData(ex.raw.X, obs=ex.obs, var=ex.raw.var).write('ex.h5ad')  # 36219 × 26020

## Load the required libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors=FALSE)
options(future.globals.maxSize=1024^10)
library(anndata)

## Data input & processing and initialization of CellChat object
## gene expression & cell label
## gene × cell (normalized count)
## Starting from an Anndata object

future::plan("multisession", workers=4)

run_cellchat <- function(sample='wt'){
    ad <- read_h5ad(paste0(sample, '.h5ad'))  # Would you like to create a default Python environment for the reticulate package? (Yes/no/cancel) cancel
    counts <- t(as.matrix(ad$X))
    data.input <- as(counts, "dgCMatrix")
    rownames(data.input) <- rownames(ad$var)
    colnames(data.input) <- rownames(ad$obs)
    meta <- ad$obs 
    meta$labels <- meta[["cell_anno"]]
    
    ## Create a CellChat object
    cellchat <- createCellChat(object=data.input, meta=meta, group.by="labels")
    
    ## Set the ligand-receptor interaction database
    CellChatDB <- CellChatDB.mouse
    CellChatDB.use <- subsetDB(CellChatDB)
    cellchat@DB <- CellChatDB.use
    
    ## Preprocessing the expression data for cell-cell communication analysis
    cellchat <- subsetData(cellchat)
    # future::plan("multisession", workers=8)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    # The number of highly variable ligand-receptor pairs used for signaling inference is *
    
    ## Compute the communication probability and infer cellular communication network
    cellchat <- computeCommunProb(cellchat, type="triMean")
    cellchat <- filterCommunication(cellchat, min.cells=10)
    
    ## Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    ## Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    # groupSize <- as.numeric(table(cellchat@idents))
    # pdf('test.pdf')
    # par(mfrow=c(1,2), xpd=TRUE)
    # netVisual_circle(cellchat@net$count, vertex.weight=groupSize, weight.scale=TRUE, label.edge=FALSE, title.name="Number of interactions")
    # netVisual_circle(cellchat@net$weight, vertex.weight=groupSize, weight.scale=TRUE, label.edge=FALSE, title.name="Interaction weights/strength")
    # dev.off()
    
    saveRDS(cellchat, file = paste0(sample, '_cellchat.rds'))
}

run_cellchat('wt')
run_cellchat('ex')


#### Comparison analysis of multiple datasets using CellChat
library(CellChat)
library(patchwork)

## Load CellChat object of each dataset and merge them together
cellchat.wt <- readRDS('wt_cellchat.rds')
cellchat.ex <- readRDS('ex_cellchat.rds')
object.list <- list(wt=cellchat.wt, ex=cellchat.ex)
cellchat <- mergeCellChat(object.list, add.names=names(object.list))

# save(object.list, file="cellchat_object.list_wt_ex.RData")
# save(cellchat, file="cellchat_merged_wt_ex.RData")

## Compare the total number of interactions and interaction strength
pdf('cci_cnt_strength_total.pdf')
gg1 <- compareInteractions(cellchat, show.legend=FALSE, group=c(1,2))
gg2 <- compareInteractions(cellchat, show.legend=FALSE, group=c(1,2), measure="weight")
gg1 + gg2
dev.off()

## Compare the number of interactions and interaction strength among different cell populations
pdf('cci_cnt_strength_cell_types.pdf')
par(mfrow=c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale=TRUE)
netVisual_diffInteraction(cellchat, weight.scale=TRUE, measure="weight")
dev.off()

## Identify dysfunctional signaling by comparing the communication probabities
## up-regulated
pdf('cci_up_fibro_smc.pdf', width=5, height=20)
gg1 <- netVisual_bubble(cellchat, sources.use='Fibroblast', targets.use='SMC', comparison=c(1, 2), max.dataset=2, title.name="Increased signaling in AD", angle.x=45, remove.isolate=TRUE)
gg1
dev.off()
## down-regulated
pdf('cci_down_fibro_smc.pdf', width=5, height=10)
gg2 <- netVisual_bubble(cellchat, sources.use='Fibroblast', targets.use='SMC', comparison=c(1, 2), max.dataset=1, title.name="Decreased signaling in AD", angle.x=45, remove.isolate=TRUE)
gg2
dev.off()


# #### Access the ligand-receptor interaction information in CellChatDB
# CellChatDB <- CellChatDB.mouse
# write.table(CellChatDB$interaction, file="interaction_input_CellChatDB.csv")
# write.table(CellChatDB$complex, file="complex_input_CellChatDB.csv")
# write.table(CellChatDB$cofactor, file="cofactor_input_CellChatDB.csv")
# write.table(CellChatDB$geneInfo, file="geneInfo_input_CellChatDB.csv")
# # CTHRC1-ADAM9 interaction not found


#### workdir: /fs/home/jiluzhang/xiaoqing20231204/Science/Figure_4/wt_ko
#### plot smc proportion vesus others
import scanpy as sc
import pandas as pd
from plotnine import *
import matplotlib.pyplot as plt
import numpy as np

adata = sc.read_h5ad('./res.h5ad')

wt_cnt = sum(adata.obs['wt_or_ko']=='wt')   # 15937
ko_cnt = sum(adata.obs['wt_or_ko']=='ko')   # 11980

wt_cnt_smc = sum((adata.obs['wt_or_ko']=='wt') & (adata.obs['cell_anno']=='SMC'))   # 875
ko_cnt_smc = sum((adata.obs['wt_or_ko']=='ko') & (adata.obs['cell_anno']=='SMC'))   # 1087

df = pd.DataFrame({'cls': ['SMC']*2+['Others']*2,
                   'trt': ['wt', 'ko']*2,
                   'val': 0})
df.loc[(df['cls']=='SMC') & (df['trt']=='wt'), 'val'] = wt_cnt_smc/wt_cnt
df.loc[(df['cls']=='SMC') & (df['trt']=='ko'), 'val'] = ko_cnt_smc/ko_cnt
df.loc[(df['cls']=='Others') & (df['trt']=='wt'), 'val'] = 1-wt_cnt_smc/wt_cnt
df.loc[(df['cls']=='Others') & (df['trt']=='ko'), 'val'] = 1-ko_cnt_smc/ko_cnt
df.to_csv('smc_vs_others_proportion.txt', index=False, sep='\t')

df = pd.read_table('smc_vs_others_proportion.txt')
df['cls'] = pd.Categorical(df['cls'], categories=['Others', 'SMC'])
df['trt'] = pd.Categorical(df['trt'], categories=['wt', 'ko'])

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

wt_cnt = sum(adata.obs['wt_or_ko']=='wt')   # 15937
ko_cnt = sum(adata.obs['wt_or_ko']=='ko')   # 11980

wt_cnt_smc = sum((adata.obs['wt_or_ko']=='wt') & (adata.obs['cell_anno']=='SMC'))   # 875
ko_cnt_smc = sum((adata.obs['wt_or_ko']=='ko') & (adata.obs['cell_anno']=='SMC'))   # 1087
wt_cnt_fibro = sum((adata.obs['wt_or_ko']=='wt') & (adata.obs['cell_anno']=='Fibroblast'))   # 6223
ko_cnt_fibro = sum((adata.obs['wt_or_ko']=='ko') & (adata.obs['cell_anno']=='Fibroblast'))   # 8995
wt_cnt_ec = sum((adata.obs['wt_or_ko']=='wt') & (adata.obs['cell_anno']=='EC'))   # 668
ko_cnt_ec = sum((adata.obs['wt_or_ko']=='ko') & (adata.obs['cell_anno']=='EC'))   # 448

df = pd.DataFrame({'cls': ['SMC']*2+['Fibroblast']*2+['EC']*2,
                   'trt': ['wt', 'ko']*3,
                   'val': 0})
df.loc[(df['cls']=='SMC') & (df['trt']=='wt'), 'val'] = wt_cnt_smc/wt_cnt
df.loc[(df['cls']=='SMC') & (df['trt']=='ko'), 'val'] = ko_cnt_smc/ko_cnt
df.loc[(df['cls']=='Fibroblast') & (df['trt']=='wt'), 'val'] = wt_cnt_fibro/wt_cnt
df.loc[(df['cls']=='Fibroblast') & (df['trt']=='ko'), 'val'] = ko_cnt_fibro/ko_cnt
df.loc[(df['cls']=='EC') & (df['trt']=='wt'), 'val'] = wt_cnt_ec/wt_cnt
df.loc[(df['cls']=='EC') & (df['trt']=='ko'), 'val'] = ko_cnt_ec/ko_cnt

df['cls'] = pd.Categorical(df['cls'], categories=['SMC', 'Fibroblast', 'EC'])
df['trt'] = pd.Categorical(df['trt'], categories=['wt', 'ko'])

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='cls', y='val', fill='trt')) + geom_bar(stat='identity', position='dodge', width=0.75) + xlab('') + ylab('Fraction') +\
                                                    scale_y_continuous(limits=[0, 0.8], breaks=np.arange(0, 0.8+0.1, 0.2)) + theme_bw()
p.save(filename='smc_fibro_ec_proportion.pdf', dpi=600, height=3, width=3)


#### smc subclusters
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42

adata = sc.read_h5ad('res.h5ad')
smc = adata[adata.obs['cell_anno']=='SMC'].copy()  # 1962 × 2346

sc.tl.leiden(smc, resolution=0.45)
sc.pl.umap(smc, color='leiden', legend_fontsize='medium', legend_loc='right margin', size=20, 
           title='', frameon=True, save='_smc_leiden.pdf')

smc.obs['leiden'].value_counts()
# 0    626
# 1    608
# 2    295
# 3    229
# 4    120
# 5     56
# 6     28

marker_genes = ['Acta2', 'Actg1', 'Tagln2', 'Myl12a', 'Tubb2a',
                'Ccnl2', 'Ccnt2', 'Malat1', 'Ccn1', 'mt-Rnr1',
                'Col1a2', 'Lum', 'Dcn', 'Gsn', 'Cygb',
                'Notch3', 'Atf3', 'Rgs4', 'Cacna1h', 'Pln']
# "Could not find keys '['Ccn1']' in columns of `adata.obs` or in adata.raw.var_names."

marker_genes = ['Acta2', 'Actg1', 'Tagln2', 'Myl12a', 'Tubb2a',
                'Ccnl2', 'Ccnt2', 'Malat1', 'mt-Rnr1',
                'Col1a2', 'Lum', 'Dcn', 'Gsn', 'Cygb',
                'Notch3', 'Atf3', 'Rgs4', 'Cacna1h', 'Pln']

del smc.uns['dendrogram_leiden']
sc.pl.dotplot(smc, marker_genes, groupby='leiden', standard_scale='var', dendrogram=True,
              save='smc_marker_genes.pdf')

cell_anno_dict = {'0':'Fibromyocyte', '1':'Contractile SMC', '2':'Stressed SMC',
                  '3':'Proliferating SMC', '4':'Contractile SMC', '5':'Fibromyocyte',
                  '6':'Stressed SMC'}
smc.obs['cell_anno'] = smc.obs['leiden']
smc.obs['cell_anno'].replace(cell_anno_dict, inplace=True)
sc.pl.umap(smc, color='cell_anno', legend_fontsize='5', legend_loc='right margin', size=20,
           title='', frameon=True, save='_smc_cell_anno.pdf')

smc.obs['cell_anno'] = pd.Categorical(smc.obs['cell_anno'],
                                      categories=['Contractile SMC', 'Proliferating SMC', 'Fibromyocyte', 'Stressed SMC'])
marker_genes = ['Acta2', 'Actg1', 'Tagln2', 'Myl12a',
                'Ccnl2', 'Ccnt2', 'Malat1',
                'Col1a2', 'Lum', 'Col1a1',
                'Atf3', 'Rgs4', 'Cacna1h', 'Pln']
sc.pl.matrixplot(smc, marker_genes, groupby='cell_anno', cmap='coolwarm', standard_scale='var', save='smc_cell_anno_heatmap.pdf')

## output marker genes for each cluster
sc.tl.rank_genes_groups(smc, 'leiden')
pd.DataFrame(smc.uns['rank_genes_groups']['names']).head(20).to_csv('smc_cluster_marker_genes.txt', index=False, sep='\t')

smc.write('./smc_res.h5ad')

## smc proportion
wt_smc_cnt = sum(smc.obs['wt_or_ko']=='wt')   # 875
ko_smc_cnt = sum(smc.obs['wt_or_ko']=='ko')   # 1087

wt_smc_con_cnt = sum((smc.obs['cell_anno']=='Contractile SMC') & (smc.obs['wt_or_ko']=='wt'))     # 285
ko_smc_con_cnt = sum((smc.obs['cell_anno']=='Contractile SMC') & (smc.obs['wt_or_ko']=='ko'))     # 443

wt_smc_pro_cnt = sum((smc.obs['cell_anno']=='Proliferating SMC') & (smc.obs['wt_or_ko']=='wt'))   # 81
ko_smc_pro_cnt = sum((smc.obs['cell_anno']=='Proliferating SMC') & (smc.obs['wt_or_ko']=='ko'))   # 148

wt_smc_str_cnt = sum((smc.obs['cell_anno']=='Stressed SMC') & (smc.obs['wt_or_ko']=='wt'))   # 149
ko_smc_str_cnt = sum((smc.obs['cell_anno']=='Stressed SMC') & (smc.obs['wt_or_ko']=='ko'))   # 174

wt_smc_fib_cnt = sum((smc.obs['cell_anno']=='Fibromyocyte') & (smc.obs['wt_or_ko']=='wt'))   # 13
ko_smc_fib_cnt = sum((smc.obs['cell_anno']=='Fibromyocyte') & (smc.obs['wt_or_ko']=='ko'))   # 43

smc_df = pd.DataFrame({'smc': ['Contractile SMC']*2+['Proliferating SMC']*2+['Stressed SMC']*2+['Fibromyocyte']*2,
                       'trt': ['wt', 'ko']*4,
                       'val': 0})
smc_df.loc[(smc_df['smc']=='Contractile SMC') & (smc_df['trt']=='wt'), 'val'] = wt_smc_con_cnt/wt_smc_cnt
smc_df.loc[(smc_df['smc']=='Contractile SMC') & (smc_df['trt']=='ko'), 'val'] = ko_smc_con_cnt/ko_smc_cnt
smc_df.loc[(smc_df['smc']=='Proliferating SMC') & (smc_df['trt']=='wt'), 'val'] = wt_smc_pro_cnt/wt_smc_cnt
smc_df.loc[(smc_df['smc']=='Proliferating SMC') & (smc_df['trt']=='ko'), 'val'] = ko_smc_pro_cnt/ko_smc_cnt
smc_df.loc[(smc_df['smc']=='Stressed SMC') & (smc_df['trt']=='wt'), 'val'] = wt_smc_str_cnt/wt_smc_cnt
smc_df.loc[(smc_df['smc']=='Stressed SMC') & (smc_df['trt']=='ko'), 'val'] = ko_smc_str_cnt/ko_smc_cnt
smc_df.loc[(smc_df['smc']=='Fibromyocyte') & (smc_df['trt']=='wt'), 'val'] = wt_smc_fib_cnt/wt_smc_cnt
smc_df.loc[(smc_df['smc']=='Fibromyocyte') & (smc_df['trt']=='ko'), 'val'] = ko_smc_fib_cnt/ko_smc_cnt
smc_df.to_csv('smc_subcluster_proportion.txt', index=False, sep='\t')

smc.obs['cell_anno'].value_counts()
# Contractile SMC      728
# Fibromyocyte         682
# Stressed SMC         323
# Proliferating SMC    229

#### plot bar for smc subcluster proportion
import pandas as pd
from plotnine import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('smc_subcluster_proportion.txt')
smc = ['Contractile SMC', 'Proliferating SMC', 'Stressed SMC', 'Fibromyocyte']
df['smc'] = pd.Categorical(df['smc'], categories=smc)
df['trt'] = pd.Categorical(df['trt'], categories=['wt', 'ko'])

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(df, aes(x='smc', y='val', fill='trt')) + geom_bar(stat='identity', position='dodge', width=0.75) + xlab('') + ylab('Fraction') +\
                                                    scale_y_continuous(limits=[0, 0.5], breaks=np.arange(0, 0.5+0.1, 0.1)) + theme_bw()
p.save(filename='smc_subcluster_proportion.pdf', dpi=600, height=3, width=6)


# ## plot umap for specific genes
# from matplotlib.colors import LinearSegmentedColormap

# sc.pl.umap(smc[smc.obs['wt_or_ko']=='wt'], color=['Lmod1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
#            legend_fontsize='xx-small', legend_loc='right margin', size=15, vmin=0, vmax=10, title='', frameon=True)
# plt.xlim(6, 13)
# plt.ylim(-3, 12) 
# plt.savefig('./figures/umap_smc_Lmod1_wt.pdf')
# plt.close()

# sc.pl.umap(smc[smc.obs['wt_or_ko']=='ko'], color=['Lmod1'], cmap=LinearSegmentedColormap.from_list('lightgray_red', ['lightgray', 'red']),
#            legend_fontsize='xx-small', legend_loc='right margin', size=15, vmin=0, vmax=10, title='', frameon=True)
# plt.xlim(6, 13)
# plt.ylim(-3, 12) 
# plt.savefig('./figures/umap_smc_Lmod1_ko.pdf')
# plt.close()



#### find differential genes
## Contractile SMC
smc_con = smc[smc.obs['cell_anno']=='Contractile SMC'].copy()  # 728 × 2346
sc.tl.rank_genes_groups(smc_con, 'wt_or_ko')
result = smc_con.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_con_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_con_up_cnt = sum((smc_con_res['ko_logfoldchanges']>1) & (smc_con_res['ko_pvals_adj']<0.05))   # 178
smc_con_dw_cnt = sum((smc_con_res['ko_logfoldchanges']<-1) & (smc_con_res['ko_pvals_adj']<0.05))  # 253
smc_con_df_cnt = smc_con_up_cnt + smc_con_dw_cnt                                                  # 431
smc_con_res[(smc_con_res['ko_logfoldchanges']<-1) & (smc_con_res['ko_pvals_adj']<0.05)]['ko_names'].to_csv('smc_con_down_genes.txt', index=False, header=False)
smc_con_res[(smc_con_res['ko_logfoldchanges']>1) & (smc_con_res['ko_pvals_adj']<0.05)]['ko_names'].to_csv('smc_con_up_genes.txt', index=False, header=False)

## Proliferating SMC
smc_pro = smc[smc.obs['cell_anno']=='Proliferating SMC'].copy()
sc.tl.rank_genes_groups(smc_pro, 'wt_or_ko')
result = smc_pro.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_pro_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_pro_up_cnt = sum((smc_pro_res['ko_logfoldchanges']>1) & (smc_pro_res['ko_pvals_adj']<0.05))   # 179
smc_pro_dw_cnt = sum((smc_pro_res['ko_logfoldchanges']<-1) & (smc_pro_res['ko_pvals_adj']<0.05))  # 88
smc_pro_df_cnt = smc_pro_up_cnt + smc_pro_dw_cnt                                                  # 267
smc_pro_res[(smc_pro_res['ko_logfoldchanges']<-1) & (smc_pro_res['ko_pvals_adj']<0.05)]['ko_names'].to_csv('smc_pro_down_genes.txt', index=False, header=False)
smc_pro_res[(smc_pro_res['ko_logfoldchanges']>1) & (smc_pro_res['ko_pvals_adj']<0.05)]['ko_names'].to_csv('smc_pro_up_genes.txt', index=False, header=False)

## Stressed SMC
smc_str = smc[smc.obs['cell_anno']=='Stressed SMC'].copy()
sc.tl.rank_genes_groups(smc_str, 'wt_or_ko')
result = smc_str.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_str_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_str_up_cnt = sum((smc_str_res['ko_logfoldchanges']>1) & (smc_str_res['ko_pvals_adj']<0.05))   # 28
smc_str_dw_cnt = sum((smc_str_res['ko_logfoldchanges']<-1) & (smc_str_res['ko_pvals_adj']<0.05))  # 42
smc_str_df_cnt = smc_str_up_cnt + smc_str_dw_cnt                                                  # 70
smc_str_res[(smc_str_res['ko_logfoldchanges']<-1) & (smc_str_res['ko_pvals_adj']<0.05)]['ko_names'].to_csv('smc_str_down_genes.txt', index=False, header=False)
smc_str_res[(smc_str_res['ko_logfoldchanges']>1) & (smc_str_res['ko_pvals_adj']<0.05)]['ko_names'].to_csv('smc_str_up_genes.txt', index=False, header=False)

## Fibromyocyte
smc_fib = smc[smc.obs['cell_anno']=='Fibromyocyte'].copy()
sc.tl.rank_genes_groups(smc_fib, 'wt_or_ko')
result = smc_fib.uns["rank_genes_groups"]
groups = result["names"].dtype.names
smc_fib_res = pd.DataFrame({group + '_' + key: result[key][group] for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges']})
smc_fib_up_cnt = sum((smc_fib_res['ko_logfoldchanges']>1) & (smc_fib_res['ko_pvals_adj']<0.05))   # 116
smc_fib_dw_cnt = sum((smc_fib_res['ko_logfoldchanges']<-1) & (smc_fib_res['ko_pvals_adj']<0.05))  # 1465
smc_fib_df_cnt = smc_fib_up_cnt + smc_fib_dw_cnt                                                  # 1581
smc_fib_res[(smc_fib_res['ko_logfoldchanges']<-1) & (smc_fib_res['ko_pvals_adj']<0.05)]['ko_names'].to_csv('smc_fib_down_genes.txt', index=False, header=False)
smc_fib_res[(smc_fib_res['ko_logfoldchanges']>1) & (smc_fib_res['ko_pvals_adj']<0.05)]['ko_names'].to_csv('smc_fib_up_genes.txt', index=False, header=False)

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
                                                    scale_y_continuous(limits=[0, 1500], breaks=np.arange(0, 1500+1, 500)) + theme_bw()
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
write.table(as.data.frame(smc_con_up_ego), file="smc_con_up_genes_go.txt")

write.csv(as.data.frame(smc_con_up_ego), file="smc_con_up_genes_go.csv")
terms <- c('GO:0090101', 'GO:0006941', 'GO:0055088', 'GO:0048705', 'GO:0086003',
           'GO:0003009', 'GO:0003012', 'GO:0014888', 'GO:0070252', 'GO:0043500')
pdf('smc_con_up_genes_go_selected.pdf')
dotplot(smc_con_up_ego, showCategory=smc_con_up_ego@result$Description[smc_con_up_ego@result$ID %in% terms], font.size=10)
dev.off()

## smc_con_down
smc_con_down <- read.table('smc_con_down_genes.txt')[, 1]
smc_con_down_ego <- enrichGO(gene=smc_con_down, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_con_down_genes_go.pdf')
dotplot(smc_con_down_ego, showCategory=10, font.size=10)
dev.off()
write.table(as.data.frame(smc_con_down_ego), file="smc_con_down_genes_go.txt")

## smc_pro_up
smc_pro_up <- read.table('smc_pro_up_genes.txt')[, 1]
smc_pro_up_ego <- enrichGO(gene=smc_pro_up, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_pro_up_genes_go.pdf')
dotplot(smc_pro_up_ego, showCategory=10, font.size=10)
dev.off()
write.table(as.data.frame(smc_pro_up_ego), file="smc_pro_up_genes_go.txt")

## smc_pro_down
smc_pro_down <- read.table('smc_pro_down_genes.txt')[, 1]
smc_pro_down_ego <- enrichGO(gene=smc_pro_down, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_pro_down_genes_go.pdf')
dotplot(smc_pro_down_ego, showCategory=10, font.size=10)
dev.off()
write.table(as.data.frame(smc_pro_down_ego), file="smc_pro_down_genes_go.txt")

## smc_str_up
smc_str_up <- read.table('smc_str_up_genes.txt')[, 1]
smc_str_up_ego <- enrichGO(gene=smc_str_up, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_str_up_genes_go.pdf')
dotplot(smc_str_up_ego, showCategory=10, font.size=10)
dev.off()
write.table(as.data.frame(smc_str_up_ego), file="smc_str_up_genes_go.txt")

## smc_str_down
smc_str_down <- read.table('smc_str_down_genes.txt')[, 1]
smc_str_down_ego <- enrichGO(gene=smc_str_down, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_str_down_genes_go.pdf')
dotplot(smc_str_down_ego, showCategory=10, font.size=10)
dev.off()
write.table(as.data.frame(smc_str_down_ego), file="smc_str_down_genes_go.txt")

## smc_fib_up
smc_fib_up <- read.table('smc_fib_up_genes.txt')[, 1]
smc_fib_up_ego <- enrichGO(gene=smc_fib_up, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                           pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_fib_up_genes_go.pdf')
dotplot(smc_fib_up_ego, showCategory=10, font.size=10)
dev.off()
write.table(as.data.frame(smc_fib_up_ego), file="smc_fib_up_genes_go.txt")

## smc_fib_down
smc_fib_down <- read.table('smc_fib_down_genes.txt')[, 1]
smc_fib_down_ego <- enrichGO(gene=smc_fib_down, OrgDb=org.Mm.eg.db, keyType='SYMBOL', ont='BP',
                             pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE)
pdf('smc_fib_down_genes_go.pdf')
dotplot(smc_fib_down_ego, showCategory=10, font.size=10)
dev.off()
write.table(as.data.frame(smc_fib_down_ego), file="smc_fib_down_genes_go.txt")

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

smc = sc.read_h5ad('smc_res.h5ad')

## contractile genes
dat = pd.DataFrame({'gene_id':['Mylk2', 'Ramp1', 'Cnn1', 'Actg2', 'Tagln', 'Acta2', 'Tpm2', 'Myl6', 'Myl9', 'Myh11', 'Mylk']*2,
                    'idx':['wt']*11+['ko']*11, 'val':[0]*22})
dat['idx'] = pd.Categorical(dat['idx'], categories=['wt', 'ko'])

for i in range(dat.shape[0]):
    dat['val'][i] = smc[smc.obs['wt_or_ko']==dat['idx'][i]].raw.X[:, np.argwhere(smc.raw.var.index==dat['gene_id'][i]).flatten()].toarray().flatten().mean()

p_value = stats.ttest_rel(dat[dat['idx']=='wt']['val'], dat[dat['idx']=='ko']['val'])[1]  # 0.0005793217668344407

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(dat, aes(x='idx', y='val', fill='idx')) + geom_boxplot(width=0.5, outlier_shape='', show_legend=False) +\
                                                     scale_y_continuous(limits=[0, 5], breaks=np.arange(0, 5+0.1, 1.0)) + theme_bw() +\
                                                     annotate("text", x=1.5, y=4, label=f"P = {p_value:.3f}", ha='center')
p.save(filename='smc_contraction.pdf', dpi=600, height=4, width=5)
np.median(dat[dat['idx']=='wt']['val'])   # 2.0280301570892334
np.median(dat[dat['idx']=='ko']['val'])   # 2.632713556289673

## mmp genes
dat = pd.DataFrame({'gene_id':['Mmp14', 'Mmp2', 'Mmp24', 'Mmp8', 'Mmp9', 'Mmp15', 'Mmp13']*2,
                    'idx':['wt']*7+['ko']*7, 'val':[0]*14})
dat['idx'] = pd.Categorical(dat['idx'], categories=['wt', 'ko'])

for i in range(dat.shape[0]):
    dat['val'][i] = smc[smc.obs['wt_or_ko']==dat['idx'][i]].raw.X[:, np.argwhere(smc.raw.var.index==dat['gene_id'][i]).flatten()].toarray().flatten().mean()

p_value = stats.ttest_rel(dat[dat['idx']=='wt']['val'], dat[dat['idx']=='ko']['val'])[1]  # 0.0005793217668344407

plt.rcParams['pdf.fonttype'] = 42
p = ggplot(dat, aes(x='idx', y='val', fill='idx')) + geom_boxplot(width=0.5, outlier_shape='', show_legend=False) +\
                                                     scale_y_continuous(limits=[0, 1.2], breaks=np.arange(0, 1.2+0.1, 0.2)) + theme_bw() +\
                                                     annotate("text", x=1.5, y=1.0, label=f"P = {p_value:.3f}", ha='center')
p.save(filename='smc_mmp.pdf', dpi=600, height=4, width=5)
np.median(dat[dat['idx']=='wt']['val'])   # 0.0064423177391290665
np.median(dat[dat['idx']=='ko']['val'])   # 0.00401050690561533

## contractile genes
con_genes = ['Ramp1', 'Cnn1', 'Actg2', 'Tagln', 'Acta2', 'Tpm2', 'Myl6', 'Myl9', 'Myh11', 'Mylk']  # delete Mylk2 due to low expression
dat = pd.DataFrame()
for gene in con_genes:
    dat_wt = pd.DataFrame({'gene_id':gene,
                           'idx':'wt',
                           'val':smc[smc.obs['wt_or_ko']=='wt'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat_ko = pd.DataFrame({'gene_id':gene,
                           'idx':'ko',
                           'val':smc[smc.obs['wt_or_ko']=='ko'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat = pd.concat([dat, dat_wt, dat_ko])

dat['idx'] = pd.Categorical(dat['idx'], categories=['wt', 'ko'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(dat, aes(x='gene_id', y='val', fill='idx')) + geom_boxplot(width=0.5, outlier_shape='', show_legend=False) +\
                                                         scale_y_continuous(limits=[0, 7], breaks=np.arange(0, 7+0.1, 1.0)) + theme_bw()
p.save(filename='smc_contraction_seperately.pdf', dpi=600, height=4, width=5)

for gene in con_genes:
    print(gene, 'pvalue:', stats.ttest_ind(dat[(dat['idx']=='wt') & (dat['gene_id']==gene)]['val'],
                                           dat[(dat['idx']=='ko') & (dat['gene_id']==gene)]['val'])[1])
# Ramp1 pvalue: 0.006038088164832036
# Cnn1 pvalue: 1.3298991715611946e-23
# Actg2 pvalue: 0.051593355763368876
# Tagln pvalue: 1.323864003077057e-20
# Acta2 pvalue: 6.247875540504227e-26
# Tpm2 pvalue: 1.1167557246656118e-13
# Myl6 pvalue: 1.626148533221483e-18
# Myl9 pvalue: 7.836466964622348e-34
# Myh11 pvalue: 1.3420959603956587e-28
# Mylk pvalue: 0.09060753317076656

## mmp genes
mmp_genes = ['Mmp14', 'Mmp2']  # delete Mmp24 & Mmp8 & Mmp9 & Mmp15 & Mmp13 due to low expression
dat = pd.DataFrame()
for gene in mmp_genes:
    dat_wt = pd.DataFrame({'gene_id':gene,
                           'idx':'wt',
                           'val':smc[smc.obs['wt_or_ko']=='wt'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat_ko = pd.DataFrame({'gene_id':gene,
                           'idx':'ko',
                           'val':smc[smc.obs['wt_or_ko']=='ko'].raw.X[:, np.argwhere(smc.raw.var.index==gene).flatten()].toarray().flatten()})
    dat = pd.concat([dat, dat_wt, dat_ko])

dat['idx'] = pd.Categorical(dat['idx'], categories=['wt', 'ko'])
plt.rcParams['pdf.fonttype'] = 42
p = ggplot(dat, aes(x='gene_id', y='val', fill='idx')) + geom_boxplot(width=0.5, outlier_shape='', show_legend=False) +\
                                                         scale_y_continuous(limits=[0, 4], breaks=np.arange(0, 4+0.1, 1)) + theme_bw()
p.save(filename='smc_mmp_seperately.pdf', dpi=600, height=4, width=5)

for gene in mmp_genes:
    print(gene, 'pvalue:', stats.ttest_ind(dat[(dat['idx']=='wt') & (dat['gene_id']==gene)]['val'],
                                           dat[(dat['idx']=='ko') & (dat['gene_id']==gene)]['val'])[1])
# Mmp14 pvalue: 2.8087396308808363e-24
# Mmp2 pvalue: 1.9754390667289024e-15
