library(DESeq2)
#lecture du fichier métadonnées
metadata <-read.csv("metadata.csv", sep=';')
#chercher les fichiers issus de features count
count_files<-list.files(path ="RNAseq_proj/comptage" , pattern ="_counts.txt$", full.names = TRUE )
print(count_files)

#lire les fichiers counts en ignorant l'enete
count_data<-lapply(count_files, function(file) {read.delim(file, skip =1)} )
#Créer une matrice finale en fusionnant les colonnes de comptes (colonne 7)
counts <-data.frame(sapply(count_data, function(x) x[,7]))
#Utiliser les Geneid comme noms de lignes
row.names(counts)<-count_data[[1]]$Geneid
#utiliser les phénotypes pour nommer les colonnes
colnames(counts)<-metadata$sample

head(counts)
#creation de l'objet deseqdataset
dds<-DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition )
#execution de l'analyse différentielle
dds<-DESeq(dds)
#extraction des résultats
results_comparison<-results(dds, contrast= c("condition", "mutant", "wild"))
head(results_comparison)
#tri du tableau par padj
results_ordered<-results_comparison[order(results_comparison$padj), ]
head(results_ordered, 10)
#Filtrage des Gènes Différentiellement Exprimés (DEG)
# On filtre sur les deux critères simultanément : padj < 0.05 ET Log2FC > 1
deseq_significant = subset(results_ordered, padj < 0.05 & abs(log2FoldChange > 1))

nrow(deseq_significant)

results_orderedlog10padj<- -log10(results_ordered$padj)
plot (results_ordered$log2FoldChange,
      results_orderedlog10padj,
      xlab ='log2foldchange',
      ylab= '-log10 padj',
      main= 'volcano plot',
      pch = 20 ,
      col = 'red'
      )
