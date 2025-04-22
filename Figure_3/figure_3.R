#### CellChat installation
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

## workdir: /fs/home/jiluzhang/xiaoqing20231204/cci
## cp /fs/home/jiluzhang/xiaoqing20231204/20241005/01/res.h5ad .

# import scanpy as sc
# dat = sc.read_h5ad('res.h5ad')  # 55867 × 2608
# wt = dat[(dat.obs['sample_groups'].isin(['wt-1', 'wt-2']))&(dat.obs['cell_anno']!='Others'), :].copy()
# sc.AnnData(wt.raw.X, obs=wt.obs, var=wt.raw.var).write('wt.h5ad')  # 23575 × 12467
# ex = dat[(dat.obs['sample_groups'].isin(['ex-1', 'ex-2']))&(dat.obs['cell_anno']!='Others'), :].copy()
# sc.AnnData(ex.raw.X, obs=ex.obs, var=ex.raw.var).write('ex.h5ad')  # 19094 × 12467

## Load the required libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors=FALSE)
options(future.globals.maxSize=1024^3)
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


#### Access the ligand-receptor interaction information in CellChatDB
CellChatDB <- CellChatDB.mouse
write.table(CellChatDB$interaction, file="interaction_input_CellChatDB.csv")
write.table(CellChatDB$complex, file="complex_input_CellChatDB.csv")
write.table(CellChatDB$cofactor, file="cofactor_input_CellChatDB.csv")
write.table(CellChatDB$geneInfo, file="geneInfo_input_CellChatDB.csv")
# CTHRC1-ADAM9 interaction not found
