## Rafael Valentin
# InterBIO 2024 Collaboration
# Gene Ontology workup

#------------------------------------------------
## Helper functions
#------------------------------------------------

# -- This function will extract the gene symbols
pullGeneSymbol = function(cpg_vector) {
  cat("Pulling down annotations\n")
  anots = as.data.table(minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
  
  subsetanots = anots[Name %in% cpg_vector, .(Name, chr, Probe_rs, UCSC_RefGene_Accession, UCSC_RefGene_Name)]

  methylationSitesMappedToGenes = do.call(c, lapply(subsetanots$Name, FUN = function(cpg) {
    x = subsetanots[Name == cpg, ]$UCSC_RefGene_Name
    xx = strsplit(x, ";")[[1]]
    
    if(length(xx) == 0) {
      output = cpg } else {
        output = paste0(cpg, "; ", xx[!duplicated(xx)])    
      }

  }))
}

## --- Server down offline version to run the gene2disease processing
# (pulled Aug 2023 - Currated)

offline_gene2disease = function(targs) {
  
  # Load disease gene associations
  DisGeNet_diseases = readRDS("~/Documents/Elysium/DisGeNet_diseases.rds")
  
  dat = DisGeNet_diseases[geneid %in% targs,]
  
  countreferenceList = unique(do.call(c, sapply(unique(dat$disease_class_name), FUN = function(x) {
    gsub("^   ", "", strsplit(x, ";    ")[[1]])
  })))
  
  op = NULL
  
  for(i in 1:nrow(dat)) {
    
    df = as.data.frame(t(data.frame(sum = sapply(countreferenceList, FUN = function(diseaseClass) {
      sum(grepl(diseaseClass, dat[i, disease_class_name]))
    }))))
    
    row.names(df) = i
    
    op[[i]] = df  
    cat(i, " ")
  }
  resSums = as.data.table(melt(do.call(rbind, op)))[, sum(value), by = variable]
  setorder(resSums, V1)
  resSums$variable = factor(resSums$variable, levels = resSums$variable, labels = resSums$variable, ordered = TRUE)
  
  return(resSums)
}


#------------------------------------------------
## Step 1: Load the data
#------------------------------------------------
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)
#library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(topGO)
library(ggplot2)
library(disgenet2r)
library(patchwork)
library(dplyr)


plot_subtitle = "HC2 v HC3 (Sex+GHA+MA corrected)"
data_type = "from volcano"
all.genes = limma_Sex.GA.MA_2v3 ##** CHANGE BASED ON MODEL OF INTEREST **##
all.genes = all.genes[order(all.genes$adj.P.Val), ]
if (any(colnames(all.genes) == "Row.names")) {
  rownames(all.genes) = all.genes$Row.names
  all.genes = all.genes[,2:ncol(all.genes)]
}

all.genes.anots = pullGeneSymbol(cpg_vector = rownames(all.genes))
all.genes.anots = data.frame("CpG" = sapply(strsplit(x = all.genes.anots, "; "), FUN = function(x) x[1]),
                             "gene" = sapply(strsplit(x = all.genes.anots, "; "), FUN = function(x) x[2]))
entID = clusterProfiler::bitr(all.genes.anots$gene, fromType = "ALIAS", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
all.genes.anots = merge(all.genes.anots, entID, by.x = "gene", by.y = "ALIAS")
rm(entID)
all.genes = merge(all.genes, all.genes.anots, by.x = "row.names", by.y = "CpG")
all.genes = all.genes[order(all.genes$adj.P.Val), ]
rm(all.genes.anots)

up.genes = all.genes[all.genes$adj.P.Val <= 0.05 & all.genes$logFC >= log2(1.5), ]
down.genes = all.genes[all.genes$adj.P.Val <= 0.05 & all.genes$logFC <= -log2(1.5), ]
diff.genes = all.genes[all.genes$adj.P.Val <= 0.05 & abs(all.genes$logFC) >= log2(1.5), ]


#------------------------------------------------
## Step 2: Run Gene Ontology analysis (topGO)
#------------------------------------------------

upList = factor(as.integer(all.genes$ENTREZID %in% up.genes$ENTREZID))
names(upList) = all.genes$ENTREZID
head(upList, 30)
table(upList)

dnList = factor(as.integer(all.genes$ENTREZID %in% down.genes$ENTREZID))
names(dnList) = all.genes$ENTREZID
head(dnList, 30)
table(dnList)

diffList = factor(as.integer(all.genes$ENTREZID %in% diff.genes$ENTREZID))
names(diffList) = all.genes$ENTREZID
head(diffList, 30)
table(diffList)


upGOdata <- new("topGOdata", ontology = "BP", allGenes = upList, geneSel = function(x)(x == 1), 
                nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ENTREZ")

dnGOdata <- new("topGOdata", ontology = "BP", allGenes = dnList, geneSel = function(x)(x == 1), 
                nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ENTREZ")

diffGOdata <- new("topGOdata", ontology = "BP", allGenes = diffList, geneSel = function(x)(x == 1), 
                  nodeSize = 10, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "ENTREZ")


upRes <- runTest(upGOdata, algorithm = "weight01", statistic = "fisher")
upRes
dnRes <- runTest(dnGOdata, algorithm = "weight01", statistic = "fisher")
dnRes
diffRes <- runTest(diffGOdata, algorithm = "weight01", statistic = "fisher")
diffRes


source("~/Documents/Elysium/plot_enrichment.r")
upGOdf = enrichment_df(upGOdata, upRes, numChar = 50, orderBy = "Scores")
upGOdf = upGOdf[!upGOdf$Term=="biological_process",]
dnGOdf = enrichment_df(dnGOdata, dnRes, numChar = 50, orderBy = "Scores")
dnGOdf = dnGOdf[!dnGOdf$Term=="biological_process",]
diffGOdf = enrichment_df(diffGOdata, diffRes, numChar = 50, orderBy = "Scores")
diffGOdf = diffGOdf[!diffGOdf$Term=="biological_process",]


GO.plt1 = ggplot(upGOdf[1:50,], aes(x = Count, y = Term, fill = Scores)) +
  geom_col() + scale_x_continuous(expand = c(0,0)) +
  scale_fill_continuous(low = "red", high = "blue",
                        name = "weight01\nfisher") +
  theme_bw() + ylab("Enriched GO Terms") + labs(title = "GO-BP of hypermethylated genes",
                                                subtitle = plot_subtitle)

GO.plt2 = ggplot(dnGOdf[1:50,], aes(x = Count, y = Term, fill = Scores)) +
  geom_col() + scale_x_continuous(expand = c(0,0)) +
  scale_fill_continuous(low = "red", high = "blue",
                        name = "weight01\nfisher") +
  theme_bw() + ylab("Enriched GO Terms") + labs(title = "GO-BP of hypomethylated genes",
                                                subtitle = plot_subtitle)

GO.plt2 + GO.plt1 + plot_layout(guides = "collect")




##** GenTable and Sig Genes for all GOs **
##*UP
Full_UP_topGO_results = GenTable(upGOdata, upRes, topNodes = length(upGOdata@graph@nodeData@data))
colnames(Full_UP_topGO_results)[6] = "P.Value"
Full_UP_topGO_results$ENTREZID = sapply(Full_UP_topGO_results$GO.ID, function(x){
  meh1 = genesInTerm(upGOdata, x)
  meh2 = meh1[[1]][meh1[[1]] %in% names(upList[upList==1])]
})
Full_UP_topGO_results$category = rep("BP", nrow(Full_UP_topGO_results))
Full_UP_topGO_results = relocate(Full_UP_topGO_results, category)
Sig_UP_topGO_results = Full_UP_topGO_results[Full_UP_topGO_results$P.Value <= 0.05,]
Sig_UP_topGO_results$adj_pval = p.adjust(Sig_UP_topGO_results$P.Value, method = "fdr")
Sig_UP_topGO_results = relocate(Sig_UP_topGO_results, adj_pval, .after = P.Value)
Sig_UP_topGO_results$Genes = sapply(Sig_UP_topGO_results$GO.ID, function(x){
  meh = annotate::getSYMBOL(unlist(Sig_UP_topGO_results[Sig_UP_topGO_results$GO.ID == x,"ENTREZID"], use.names = F),
                            data='org.Hs.eg')
  meh = unname(meh)
})


##*DOWN
Full_DN_topGO_results = GenTable(dnGOdata, dnRes, topNodes = length(dnGOdata@graph@nodeData@data))
colnames(Full_DN_topGO_results)[6] = "P.Value"
Full_DN_topGO_results$ENTREZID = sapply(Full_DN_topGO_results$GO.ID, function(x){
  meh1 = genesInTerm(dnGOdata, x)
  meh2 = meh1[[1]][meh1[[1]] %in% names(dnList[dnList==1])]
})
Full_DN_topGO_results$category = rep("BP", nrow(Full_DN_topGO_results))
Full_DN_topGO_results = relocate(Full_DN_topGO_results, category)
Sig_DN_topGO_results = Full_DN_topGO_results[Full_DN_topGO_results$P.Value <= 0.05,]
Sig_DN_topGO_results$adj_pval = p.adjust(Sig_DN_topGO_results$P.Value, method = "fdr")
Sig_DN_topGO_results = relocate(Sig_DN_topGO_results, adj_pval, .after = P.Value)
Sig_DN_topGO_results$Genes = sapply(Sig_DN_topGO_results$GO.ID, function(x){
  meh = annotate::getSYMBOL(unlist(Sig_DN_topGO_results[Sig_DN_topGO_results$GO.ID == x,"ENTREZID"], use.names = F),
                            data='org.Hs.eg')
  meh = unname(meh)
})


##*ALL
Full_DIFF_topGO_results = GenTable(diffGOdata, diffRes, topNodes = length(diffGOdata@graph@nodeData@data))
colnames(Full_DIFF_topGO_results)[6] = "P.Value"
Full_DIFF_topGO_results$ENTREZID = sapply(Full_DIFF_topGO_results$GO.ID, function(x){
  meh1 = genesInTerm(diffGOdata, x)
  meh2 = meh1[[1]][meh1[[1]] %in% names(diffList[diffList==1])]
})
Full_DIFF_topGO_results$category = rep("BP", nrow(Full_DIFF_topGO_results))
Full_DIFF_topGO_results = relocate(Full_DIFF_topGO_results, category)
Sig_DIFF_topGO_results = Full_DIFF_topGO_results[Full_DIFF_topGO_results$P.Value <= 0.05,]
Sig_DIFF_topGO_results$adj_pval = p.adjust(Sig_DIFF_topGO_results$P.Value, method = "fdr")
Sig_DIFF_topGO_results = relocate(Sig_DIFF_topGO_results, adj_pval, .after = P.Value)
Sig_DIFF_topGO_results$Genes = sapply(Sig_DIFF_topGO_results$GO.ID, function(x){
  meh = annotate::getSYMBOL(unlist(Sig_DIFF_topGO_results[Sig_DIFF_topGO_results$GO.ID == x,"ENTREZID"], use.names = F),
                            data='org.Hs.eg')
  meh = unname(meh)
})






meh = Sig_DN_topGO_results
for (i in 1:nrow(meh)) {
  meh[1,"Genes"] = paste(unlist(meh[i, "Genes"], use.names = F), collapse = ", ")
}
rm(i)

## Save GenTable output

#Base_model_GO_outputs = list()
#Sex.GA.MA_model_GO_outputs = list()

if (grepl("Sex", plot_subtitle)) {
  
  if (!any(ls() == "Sex.GA.MA_model_GO_outputs")) {
    Sex.GA.MA_model_GO_outputs = list()
  }
  
  Sex.GA.MA_model_GO_outputs[[paste(gsub("\\(.*", "", plot_subtitle), "Hypermethylated")]] = Sig_UP_topGO_results
  Sex.GA.MA_model_GO_outputs[[paste(gsub("\\(.*", "", plot_subtitle), "Hypomethylated")]] = Sig_DN_topGO_results
  
} else {
  
  if (!any(ls() == "Base_model_GO_outputs")) {
    Base_model_GO_outputs = list()
  }
  
  Base_model_GO_outputs[[paste(gsub("\\(.*", "", plot_subtitle), "Hypermethylated")]] = Sig_UP_topGO_results
  Base_model_GO_outputs[[paste(gsub("\\(.*", "", plot_subtitle), "Hypomethylated")]] = Sig_DN_topGO_results
  
}


openxlsx::write.xlsx(Sex.GA.MA_model_GO_outputs, file = "Sex+GA+MA Models GO outputs.xlsx")
openxlsx::write.xlsx(Base_model_GO_outputs, file = "Base Models (Not corrected) GO outputs.xlsx")



#------------------------------------------------
## Step 3: Run Disease Phenotypes
#------------------------------------------------
# -- Mapping disease phenotypes using disgenet2r
#disgenet_api_key = get_disgenet_api_key(
#  email = readline(prompt = paste0("Input account username (email): ")), 
#  password = readline(prompt = paste0("Password: ")) )

# Isolate targets
#targs_up = run_gene2disease(up.genes$ENTREZID[!is.na(up.genes$ENTREZID)])
targs_up = offline_gene2disease(up.genes$ENTREZID[!is.na(up.genes$ENTREZID)])

#targs_dn = run_gene2disease(down.genes$ENTREZID[!is.na(down.genes$ENTREZID)])
targs_dn = offline_gene2disease(down.genes$ENTREZID[!is.na(down.genes$ENTREZID)])


# Hypermethylated CpGs
ggplot(targs_up, aes(x = variable, y = V1)) + geom_bar(stat = "identity") + coord_flip() + 
  labs(y = "Count", x = "Disease Phenotype") +
  theme_bw() + theme(legend.position = "bottom") +
  labs(title = paste("Disease phenotypes of hypermethylated genes\n", plot_subtitle, data_type))

# Hypomethylated CpGs
ggplot(targs_dn, aes(x = variable, y = V1)) + geom_bar(stat = "identity") + coord_flip() + 
  labs(y = "Count", x = "Disease Phenotype") + 
  scale_linetype_manual(values=c("dashed", "solid")) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(title = paste("Disease phenotypes of hypomethylated genes\n", plot_subtitle, data_type))






