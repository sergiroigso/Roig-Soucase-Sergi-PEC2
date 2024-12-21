######### Paquetes

# Instalar BiocManager si no está disponible
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Lista de paquetes de Bioconductor
bioconductor_packages <- c(
  "clusterProfiler", "org.Mm.eg.db", "GO.db", "limma", "oligo",
  "Biobase", "arrayQualityMetrics", "genefilter", "mouse4302.db",
  "AnnotationDbi"
)

# Instalar paquetes de Bioconductor
for (pkg in bioconductor_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Lista de paquetes de CRAN
cran_packages <- c("R.utils", "pheatmap", "pvca", "ggplot2")

# Instalar paquetes de CRAN
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(clusterProfiler)
library(org.Mm.eg.db)
library(clusterProfiler)
library(GO.db)
library(limma)
library(R.utils)
library(oligo)
library(Biobase)
library(arrayQualityMetrics)
library(genefilter)
library(mouse4302.db)
library(AnnotationDbi)
library(pheatmap)
library(pvca)
library(ggplot2)

######### Preparación de los datos

# creamos un working directory
workingDir <- getwd()
# ruta donde se encuentran los datos
dataDir <- file.path(workingDir, "datos")
# ruta donde se depositarán los resultados
resultsDir <- file.path(workingDir, "resultados")

# función selectSamples.R:
filter_microarray <- function(allTargets, seed = 123) {
  # Configurar la semilla aleatoria
  set.seed(21012126)
  
  # Filtrar las filas donde 'time' no sea 'hour 2'
  filtered <- subset(allTargets, time != "hour 2")
  
  # Dividir el dataset por grupos únicos de 'infection' + 'agent'
  filtered$group <- interaction(filtered$infection, filtered$agent)
  
  # Seleccionar 4 muestras al azar de cada grupo
  selected <- do.call(rbind, lapply(split(filtered, filtered$group), function(group_data) {
    if (nrow(group_data) > 4) {
      group_data[sample(1:nrow(group_data), 4), ]
    } else {
      group_data
    }
  }))
  
  # Obtener los índices originales como nombres de las filas seleccionadas
  original_indices <- match(selected$sample, allTargets$sample)
  
  # Modificar los rownames usando 'sample' y los índices originales
  rownames(selected) <- paste0(selected$sample, ".", original_indices)
  
  # Eliminar la columna 'group' y devolver el resultado
  selected$group <- NULL
  return(selected)
}
# Simular el dataset basado en la descripción proporcionada
allTargets <- data.frame(
  sample = c("GSM944831", "GSM944838", "GSM944845", "GSM944852", "GSM944859",
             "GSM944833", "GSM944840", "GSM944847", "GSM944854", "GSM944861",
             "GSM944834", "GSM944841", "GSM944848", "GSM944855", "GSM944862",
             "GSM944832", "GSM944839", "GSM944846", "GSM944853", "GSM944860",
             "GSM944835", "GSM944842", "GSM944849", "GSM944856", "GSM944863",
             "GSM944836", "GSM944843", "GSM944850", "GSM944857", "GSM944864",
             "GSM944837", "GSM944844", "GSM944851", "GSM944858", "GSM944865"),
  infection = c(rep("uninfected", 15), rep("S. aureus USA300", 20)),
  time = c(rep("hour 0", 15), rep("hour 2", 5), rep("hour 24", 15)),
  agent = c(rep("untreated", 5), rep("linezolid", 5), rep("vancomycin", 5),
            rep("untreated", 5), rep("untreated", 5), rep("linezolid", 5), rep("vancomycin", 5))
)
# Aplicar la función (cambiar 123 por vuestro ID de la UOC u otro número que podáis escribir en el documento)
result <- filter_microarray(allTargets, seed=21012126)

cambio <- T
library(R.utils)
if (cambio) {
  # RENOMBRAMOS Y DESCOMPRIMIMOS LOS ARCHIVOS
  # Obtener la lista de archivos en el directorio actual que terminan en .CEL.gz
  archivos <- list.files(path = "datos/", pattern = "\\.CEL\\.gz$", full.names = TRUE)
  # Renombrar los archivos usando gsub para extraer la parte inicial GSMxxxxxx y añadir ".CEL.gz"
  archivos_nuevos <- gsub("^(GSM\\d+).*\\.CEL\\.gz$", "\\1.CEL.gz", basename(archivos))
  # Ver los nombres originales y nuevos
  data.frame(Original = basename(archivos), Nuevo = archivos_nuevos)
  # Renombrar los archivos en el sistema
  file.rename(from = archivos, to = file.path(dirname(archivos), archivos_nuevos))
  # Listar archivos .CEL.gz en tu directorio
  gz_files <- list.files(path = "datos/", pattern = "\\.CEL\\.gz$", full.names = TRUE)
  # Descomprimir todos los archivos
  sapply(gz_files, gunzip, remove = TRUE)
  # nombramos aquellos datos seleccionados
  selected_files <- paste0(result$sample, ".CEL")
  # Mover los archivos seleccionados a la nueva carpeta
  # Definir las rutas completas de origen y destino
  destination_dir <- file.path(workingDir, "datos/selected/")
  # Mover los archivos seleccionados de la carpeta "datos" a la carpeta "datos_selected"
  for (file in selected_files) {
    # Generar las rutas completas de los archivos
    source_path <- file.path(dataDir, file)
    destination_path <- file.path(destination_dir, file)
    # Verificar si el archivo existe en la carpeta de origen
    if (file.exists(source_path)) {
      # Mover el archivo a la carpeta destino
      file.rename(source_path, destination_path)
      message(paste("Movido:", file))
    } else {
      message(paste("El archivo", file, "no existe en la carpeta de origen"))
    }
  }
}


library(oligo)
library(Biobase)
# carpeta donde están los archivos seleccionados
destination_dir <- file.path(workingDir, "datos/selected/")
# lista de los archivos .CEL seleccionados
CELfiles <- list.files(path = "datos/selected/", full.names = F)
# nombre de los sample
sampleNames <- result$sample
targets <- AnnotatedDataFrame(result)
rawData <- read.celfiles(file.path(destination_dir, CELfiles), phenoData = targets)
# guardamos el archivo en un .rda
save(rawData, file = "GSE38531_raw.Rda")
rawData

######### Análisis exploratorio y control de calidad

# Boxplot
boxplot(rawData,,which = "all", cex.axis = 0.6,las = 2,
        main = "Intensity distribution of RAW data", names = sampleNames)

# HIERARQUICAL CLUSTERING
clust.euclid.average <- hclust(dist(t(exprs(rawData))), method = "average")
plot(clust.euclid.average, labels = sampleNames, main = "Hierarchical clustering of RawData",
     cex = 0.7, hang = -1)

# PRINCIPAL COMPONENT ANALYSIS

library(ggplot2)

plotPCA_ggplot <- function(X, labels = NULL, colors = NULL, dataDesc = "", scale = FALSE,
                           formapunts = NULL, myCex = 0.8, ...) {
  # realizamos PCA
  pcX <- prcomp(t(X), scale. = scale)
  
  # extraemos las proporciones de varianza
  loads <- round(pcX$sdev^2 / sum(pcX$sdev^2) * 100, 1)
  pcData <- as.data.frame(pcX$x[, 1:2])  # toma solo los primeros 2 PCs
  pcData$labels <- labels
  pcData$colors <- colors
  
  # nombres de los ejes
  xlab <- paste("PC1", loads[1], "%")
  ylab <- paste("PC2", loads[2], "%")
  
  # Crear el gráfico con ggplot2
  p <- ggplot(pcData, aes(x = PC1, y = PC2, color = colors, shape = colors)) +
    geom_point(size = 3) +  # Puntos
    geom_text(aes(label = labels), vjust = -0.5, cex = myCex) +  # etiquetas
    labs(title = paste("Plot of first 2 PCs for expressions in", dataDesc),
         x = xlab, y = ylab) +
    scale_color_manual(values = c("red", "blue", "orange", "purple")) +  # colores
    theme_minimal() +  # Tema limpio
    theme(legend.position = "top",  # Posición de la leyenda
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
  
  print(p)  # muestra el gráfico
}

# usamos la función con los datos crudos
plotPCA_ggplot(exprs(rawData), labels = sampleNames, dataDesc = "Raw data", 
               colors = as.factor(result$agent), formapunts = c(rep(16, 4), rep(17, 4)), myCex = 0.8)



# SAVE TO A FILE
pdf(file.path(resultsDir, "QCPlots_Raw.pdf"))
boxplot(rawData,,which = "all", cex.axis = 0.6,las = 2,
        main = "Intensity distribution of RAW data", names = sampleNames)
plot(clust.euclid.average, labels = sampleNames, main = "Hierarchical clustering of samples of RawData",
     cex = 0.7, hang = -1)
plotPCA_ggplot(exprs(rawData), labels = sampleNames, dataDesc = "Raw data", 
               colors = as.factor(result$agent), formapunts = c(rep(16, 4), rep(17, 4)), myCex = 0.8)
dev.off()

######### Paquete arrayqualitymetrics

library("arrayQualityMetrics")
# Avoid re-running it each time the script is executed.
rerun <- F
if (rerun) {
  arrayQualityMetrics(rawData, reporttitle = "QC_RawData", force = TRUE)
}

######### Normalización

eset <- rma(rawData)
write.exprs(eset, file.path(resultsDir, "NormData.txt"))
boxplot(eset)
eset
# mediante boxplot:
boxplot(eset,,which = "all", cex.axis = 0.6,las = 2,
        main = "Intensity distribution of NORMALIZED data", names = sampleNames)

# Ahora, crear la gráfica con ggplot2
plotPCA_ggplot(exprs(eset), labels = sampleNames, dataDesc = "Normalized data", 
               colors = as.factor(result$agent), formapunts = c(rep(16, 4), rep(17, 4)), myCex = 0.8)


# Estudio efecto Batch

library(pvca)

# Asegúrate de que los factores son categóricos
pData(eset)$infection <- as.factor(pData(eset)$infection)
pData(eset)$time <- as.factor(pData(eset)$time)
pData(eset)$agent <- as.factor(pData(eset)$agent)

# Configurar factores y parámetros
pct_threshold <- 0.6
batch.factors <- c("infection", "time", "agent")

# Ejecutar PVCA
pvcaObj <- pvcaBatchAssess(eset, batch.factors, pct_threshold)

# Graficar resultados
bp <- barplot(
  pvcaObj$dat, 
  xlab = "Effects", 
  ylab = "Weighted average proportion variance", 
  main = "PVCA estimation",
  ylim= c(0,1),col = c("mediumorchid"), las=2
)
axis(1, at = bp, labels = pvcaObj$label, cex.axis = 0.75, las = 2)
text(bp, pvcaObj$dat, labels = round(pvcaObj$dat, 3), pos = 3, cex = 0.7)


######### Filtrado

library(genefilter)
library(mouse4302.db)
annotation(eset) <- "mouse4302"
eset_filtered <- nsFilter(eset, var.func = IQR, var.cutoff = 0.75, var.filter = TRUE, filterByQuantile = TRUE, require.entrez = TRUE)
# NUMBER OF GENES IN
print(eset_filtered$eset)
# extraemos la matriz de datos filtrados
filteredEset <- eset_filtered$eset
filteredData <- exprs(filteredEset)
colnames(filteredData) <- pData(eset_filtered$eset)$sample
print(eset_filtered$filter.log)

# guardamos los resultados normalizados y filtrados
write.csv(exprs(eset), file="./resultados/normalized.Data.csv")
write.csv(exprs(filteredEset), file="./resultados/normalized.Filtered.Data.csv")
save(eset, filteredEset, file="./resultados/normalized.Data.Rda")

######### Construcción de las matrices de diseño y de contrastes
library(limma)

# Creamos un factor combinado: infection + agent
group <- paste(pData(filteredEset)$infection, pData(filteredEset)$agent, sep = "_")
group <- make.names(group)
group <- factor(group)
# Creamos la matriz de diseño
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
rownames(design) <- sampleNames(filteredEset)
# Mostramos la matriz de diseño
print(design)

# Crear la matriz de contrastes
cont.matrix <- makeContrasts(
  Infected_vs_Uninfected_Untreated = S..aureus.USA300_untreated - uninfected_untreated,
  Infected_vs_Uninfected_Linezolid  = S..aureus.USA300_linezolid - uninfected_linezolid,
  Infected_vs_Uninfected_Vancomycin = S..aureus.USA300_vancomycin - uninfected_vancomycin,
  levels = design
)
# Mostrar la matriz de contrastes
print(cont.matrix)

# estimamos el modelo

fit<-lmFit(filteredEset, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)

# comparación entre listas de genes
# adjust.method = "fdr" ajusta los valores p usando el método FDR. 
# Si los valores ajustados son grandes, no se detectarán genes significativos.
res<-decideTests(fit.main, method="separate", adjust.method="none", p.value=0.1, lfc=1)
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))

# podemos hacer un diagrama de venn:
vennDiagram(res.selected[,1:3], include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"), cex=c(1,1,1))
title("Genes in common between the three comparisons\n Genes selected with FDR < 0.1 and logFC > 1")

######### Obtención de las listas de genes diferencialmente expresados para cada comparación

# obtenemos las listas de los DEGs sin anotar
topTab1 <- topTable(fit.main, number = nrow(fit.main), coef ="Infected_vs_Uninfected_Untreated", adjust = "fdr")
topTab2 <- topTable(fit.main, number = nrow(fit.main), coef ="Infected_vs_Uninfected_Linezolid", adjust = "fdr")
topTab3 <- topTable(fit.main, number = nrow(fit.main), coef ="Infected_vs_Uninfected_Vancomycin", adjust = "fdr")

######### Anotacion de los genes

# anotamos las listas mediante la función annotatedTopTable
annotatedTopTable <- function(topTab, anotPackage){
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
  return(annotatedTopTab)
  }

# generamos las listas anotadas
topAnnotated_Untreated <- annotatedTopTable(topTab1, anotPackage="mouse4302.db")
topAnnotated_Linezolid <- annotatedTopTable(topTab2, anotPackage="mouse4302.db")
topAnnotated_Vancomycin <- annotatedTopTable(topTab3, anotPackage="mouse4302.db")

# escribimos los resultados
library(xtable)
write.csv2(topAnnotated_Untreated, file = file.path(resultsDir, "Genes_UNTREATED.csv"))
print(xtable(topTab1, align = "lllllll"), type = "html", html.table.attributes = "",
      file = file.path(resultsDir, "Genes_UNTREATED.html"))
write.csv2(topAnnotated_Linezolid, file = file.path(resultsDir, "Genes_LINEZOLID.csv"))
print(xtable(topTab2, align = "lllllll"), type = "html", html.table.attributes = "",
      file = file.path(resultsDir, "Genes_LINEZOLID.html"))
write.csv2(topAnnotated_Vancomycin, file = file.path(resultsDir, "Genes_VANCOMYCIN.csv"))
print(xtable(topTab3, align = "lllllll"), type = "html", html.table.attributes = "",
      file = file.path(resultsDir, "Genes_VANCOMYCIN.html"))

# volcano plots
library(mouse4302.db)
geneSymbols <- select(mouse4302.db, rownames(fit.main), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL

volcanoplot(fit.main, coef=1, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))

volcanoplot(fit.main, coef=2, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[2], sep="\n"))
abline(v=c(-1,1))

volcanoplot(fit.main, coef=3, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[3], sep="\n"))
abline(v=c(-1,1))

# heatmap

library(pheatmap)
library(mouse4302.db)
library(AnnotationDbi)

# Obtener los genes diferencialmente expresados de las comparaciones
selected_genes <- unique(c(
  rownames(topTab1[topTab1$P.Value < 0.05 & abs(topTab1$logFC) > 1.5, ]),
  rownames(topTab2[topTab2$P.Value < 0.05 & abs(topTab2$logFC) > 1.5, ]),
  rownames(topTab3[topTab3$P.Value < 0.05 & abs(topTab3$logFC) > 1.5, ])
))


# Obtener los símbolos de los genes para los probes seleccionados
symbols <- select(mouse4302.db, keys = selected_genes, columns = "SYMBOL", keytype = "PROBEID")

# Filtrar para eliminar cualquier valor NA o duplicado en la columna SYMBOL
symbols <- symbols[!is.na(symbols$SYMBOL), ]
symbols <- unique(symbols$SYMBOL)
# Filtrar la matriz de expresión con los genes seleccionados
heatmap_data <- exprs(filteredEset)[selected_genes, ]
# las filas de la matriz de expresión coincidan con los símbolos
rownames(heatmap_data) <- symbols
# Escalar las filas (genes) para una mejor visualización
scaled_data <- t(scale(t(heatmap_data)))
# Crear el heatmap con los nombres de los genes anotados
pheatmap(
  scaled_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = pData(filteredEset)[, c("infection", "agent")],
  main = "Heatmap of Differentially Expressed Genes",
  color = colorRampPalette(c("blue", "white", "red"))(50)
)


######### Análisis de la significación biológica

library(AnnotationDbi)

listOfTables <- list(Infected_vs_Uninfected_Untreated = topTab1, 
                     Infected_vs_Uninfected_Linezolid  = topTab2, 
                     Infected_vs_Uninfected_Vancomycin = topTab3)
# iniciamos una lista vacía para guardar los resultados
listOfSelected <- list()
# seleccionamos los genes
for (i in 1:length(listOfTables)) {
  # seleccionamos la toptable
  topTab <- listOfTables[[i]]
  # seleccionamos los genes que superen el umbral estadístico
  whichGenes <- topTab[,"P.Value"] < 0.1
  selectedIDs <- rownames(topTab)[whichGenes]
  # convertimos los IDs de sondas en IDs de Entrez usando AnnotationDbi::select
  EntrezIDs <- AnnotationDbi::select(
    mouse4302.db,
    keys = selectedIDs,          # Input probe IDs
    columns = "ENTREZID",        # Column to retrieve
    keytype = "PROBEID"          # Type of input IDs
  )
  # extraemos solo los Entrez IDs (quitamos duplicados y NAs)
  EntrezIDs <- unique(na.omit(EntrezIDs$ENTREZID))
  # guardamos los Entrez IDs en una lista
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}
# mostramos el número de genes por condición
sapply(listOfSelected, length)


# Cargar bibliotecas necesarias
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(ggplot2)

# Crear un directorio para guardar los resultados
if (!dir.exists("resultados")) dir.create("resultados")

# Bucle para realizar análisis de enriquecimiento en GO, KEGG y Reactome
for (i in 1:length(listOfSelected)) {
  comparisonName <- names(listOfSelected)[i]
  genes <- listOfSelected[[i]]
  
  # Análisis GO
  ego <- enrichGO(
    gene          = genes,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  # Análisis KEGG
  ekegg <- enrichKEGG(
    gene          = genes,
    organism      = "mmu",
    keyType       = "ncbi-geneid",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  # Verifica si hay resultados antes de generar los gráficos
  if (nrow(ego) > 0) {
    go_plot <- dotplot(ego, showCategory = 10, title = paste("GO Enrichment (BP):", comparisonName))
    ggsave(filename = paste0("resultados/GO_", comparisonName, ".png"), plot = go_plot, width = 10, height = 8, dpi = 300)
    write.table(as.data.frame(ego), file = paste0("resultados/GO_", comparisonName, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    cat("No GO enrichment results for", comparisonName, "\n")
  }
  
  if (nrow(ekegg) > 0) {
    kegg_plot <- dotplot(ekegg, showCategory = 10, title = paste("KEGG Enrichment:", comparisonName))
    ggsave(filename = paste0("resultados/KEGG_", comparisonName, ".png"), plot = kegg_plot, width = 10, height = 8, dpi = 300)
    write.table(as.data.frame(ekegg), file = paste0("resultados/KEGG_", comparisonName, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    cat("No KEGG enrichment results for", comparisonName, "\n")
  }
  
  # Imprimir resumen
  cat("\n### Resultados para", comparisonName, "###\n")
  if (nrow(ego) > 0) print(head(ego))
  if (nrow(ekegg) > 0) print(head(ekegg))
}