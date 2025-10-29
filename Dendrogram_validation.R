#Working directory
setwd("/Users/martinbarreda/Desktop/Cecilia")

# libraries
library(factoextra)
library(pvclust)
library(dendextend)
library(cluster)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(ggdendro)
library(reshape2)
library(clValid)

# Data must be normalized before this statistics

### raw data----

data <- read.csv("sidas_rt_mz_progenesis_data_normalized.csv", row.names = 1,
                 check.names = FALSE)

labels <- as.character(data[1, ])
cat("\nFind groups:", unique(labels), "\n")

# Eliminar la fila de labels de los datos
data_numeric <- data[-1, ]

# Convert to numeric 
data_numeric <- apply(data_numeric, 2, as.numeric)
rownames(data_numeric) <- rownames(data)[-1]

cat("Data dimentions:", nrow(data_numeric), "metabolites x", 
    ncol(data_numeric), "samples\n\n")

# Transpose for sample analysis (common in metabolomics)
# now: rows = samples, columns = metabolites
data_t <- t(data_numeric)

# verify dimentions
cat("\nData dimentions:\n")
cat("Samples:", nrow(data_t), "\n")
cat("Metabolites:", ncol(data_t), "\n")

# ===== 2. COPHENETIC COEFFICIENT =====
cat("===== COPHENETIC COEFFICIENT =====\n")

metodos <- c("ward.D2")
resultados_cof <- data.frame(
  Metodo = metodos, 
  Coefenetico = NA,
  stringsAsFactors = FALSE
)

for(i in 1:length(metodos)){
  hc <- hclust(dist(data_t), method = metodos[i])
  coph <- cophenetic(hc)
  resultados_cof$Coefenetico[i] <- cor(dist(data_t), coph)
}

print(resultados_cof)

# Suggest optimal number based on silhouette
sil_scores <- sapply(2:10, function(k) {
  hc <- hclust(dist(data_t), method = metodos)
  clusters <- cutree(hc, k = k)
  mean(silhouette(clusters, dist(data_t))[, 3])
})
k_optimo <- which.max(sil_scores) + 1

cat("\n*** SUGGESTED OPTIMAL NUMBER:", k_optimo, "clusters ***\n")
cat("(Based on average silhouette)\n\n")

# ===== 5. VALIDATION INDICES =====
cat("===== VALIDATION INDICES FOR k =", k_optimo, "=====\n")

hc_final <- hclust(dist(data_t), method = metodos)
clusters <- cutree(hc_final, k = k_optimo)

# silhouette
sil <- silhouette(clusters, dist(data_t))
sil_promedio <- mean(sil[, 3])

cat("\n*** SILHOUETTE:", round(sil_promedio, 3), "***\n")
if(sil_promedio > 0.70) {
  cat("✓ EXCELLENT: Very strong clustering structure\n")
} else if(sil_promedio > 0.50) {
  cat("✓ GOOD: Clear clustering structure\n")
} else if(sil_promedio > 0.25) {
  cat("⚠ REGULAR: Weak clustering structure\n")
} else {
  cat("✗ POOR: There is no clear clustering structure\n")
}

# Dunn's Index
dunn_index <- dunn(distance = dist(data_t), clusters = clusters)
cat("\n*** Dunn's Index:", round(dunn_index, 3), "***\n")
cat("(Higher values indicate better separation)\n\n")

### label data----

data <- read.csv("sidas_identificaciones_data_normalized.csv",
                 row.names = 1,
                 check.names = FALSE)

labels <- as.character(data[1, ])
cat("\nFind groups:", unique(labels), "\n")

# Eliminar la fila de labels de los datos
data_numeric <- data[-1, ]

# Convert to numeric 
data_numeric <- apply(data_numeric, 2, as.numeric)
rownames(data_numeric) <- rownames(data)[-1]

cat("Data dimentions:", nrow(data_numeric), "metabolites x", 
    ncol(data_numeric), "samples\n\n")

# Transpose for sample analysis (common in metabolomics)
# now: rows = samples, columns = metabolites
data_t <- t(data_numeric)

# verify dimentions
cat("\nData dimentions:\n")
cat("Samples:", nrow(data_t), "\n")
cat("Metabolites:", ncol(data_t), "\n")

# ===== 2. COPHENETIC COEFFICIENT =====
cat("===== COPHENETIC COEFFICIENT =====\n")

metodos <- c("ward.D2")
resultados_cof <- data.frame(
  Metodo = metodos, 
  Coefenetico = NA,
  stringsAsFactors = FALSE
)

for(i in 1:length(metodos)){
  hc <- hclust(dist(data_t), method = metodos[i])
  coph <- cophenetic(hc)
  resultados_cof$Coefenetico[i] <- cor(dist(data_t), coph)
}

print(resultados_cof)

# Suggest optimal number based on silhouette
sil_scores <- sapply(2:10, function(k) {
  hc <- hclust(dist(data_t), method = metodos)
  clusters <- cutree(hc, k = k)
  mean(silhouette(clusters, dist(data_t))[, 3])
})
k_optimo <- which.max(sil_scores) + 1

cat("\n*** SUGGESTED OPTIMAL NUMBER:", k_optimo, "clusters ***\n")
cat("(Based on average silhouette)\n\n")

# ===== 5. VALIDATION INDICES =====
cat("===== VALIDATION INDICES FOR k =", k_optimo, "=====\n")

hc_final <- hclust(dist(data_t), method = metodos)
clusters <- cutree(hc_final, k = k_optimo)

# silhouette
sil <- silhouette(clusters, dist(data_t))
sil_promedio <- mean(sil[, 3])

cat("\n*** SILHOUETTE:", round(sil_promedio, 3), "***\n")
if(sil_promedio > 0.70) {
  cat("✓ EXCELLENT: Very strong clustering structure\n")
} else if(sil_promedio > 0.50) {
  cat("✓ GOOD: Clear clustering structure\n")
} else if(sil_promedio > 0.25) {
  cat("⚠ REGULAR: Weak clustering structure\n")
} else {
  cat("✗ POOR: There is no clear clustering structure\n")
}

# Dunn's Index
dunn_index <- dunn(distance = dist(data_t), clusters = clusters)
cat("\n*** Dunn's Index:", round(dunn_index, 3), "***\n")
cat("(Higher values indicate better separation)\n\n")

### Sidas markers----

data <- read.csv("sidas_markers_data_normalized.csv", 
                 row.names = 1, 
                 check.names = FALSE)

# labels 
labels <- data$Label
cat("\nfind groupss:\n")
print(table(labels))

# keep just numeric data
data_numeric <- data[, -1]  

# Convert all to numeric 
data_numeric <- as.data.frame(lapply(data_numeric, as.numeric))
rownames(data_numeric) <- rownames(data)

# verify dimentions
cat("\nData dimentions:\n")
cat("Samples:", nrow(data_t), "\n")
cat("Metabolites:", ncol(data_t), "\n")

# ===== 2. COPHENETIC COEFFICIENT =====
cat("===== COPHENETIC COEFFICIENT =====\n")

metodos <- c("ward.D2")
resultados_cof <- data.frame(
  Metodo = metodos, 
  Coefenetico = NA,
  stringsAsFactors = FALSE
)

for(i in 1:length(metodos)){
  hc <- hclust(dist(data_numeric), method = metodos[i])
  coph <- cophenetic(hc)
  resultados_cof$Coefenetico[i] <- cor(dist(data_numeric), coph)
}

print(resultados_cof)

# Suggest optimal number based on silhouette
sil_scores <- sapply(2:10, function(k) {
  hc <- hclust(dist(data_numeric), method = metodos)
  clusters <- cutree(hc, k = k)
  mean(silhouette(clusters, dist(data_numeric))[, 3])
})
k_optimo <- which.max(sil_scores) + 1

cat("\n*** SUGGESTED OPTIMAL NUMBER:", k_optimo, "clusters ***\n")
cat("(Based on average silhouette)\n\n")

# ===== 5. VALIDATION INDICES =====
cat("===== VALIDATION INDICES FOR k =", k_optimo, "=====\n")

hc_final <- hclust(dist(data_numeric), method = metodos)
clusters <- cutree(hc_final, k = k_optimo)

# silhouette
sil <- silhouette(clusters, dist(data_numeric))
sil_promedio <- mean(sil[, 3])

cat("\n*** SILHOUETTE:", round(sil_promedio, 3), "***\n")
if(sil_promedio > 0.70) {
  cat("✓ EXCELLENT: Very strong clustering structure\n")
} else if(sil_promedio > 0.50) {
  cat("✓ GOOD: Clear clustering structure\n")
} else if(sil_promedio > 0.25) {
  cat("⚠ REGULAR: Weak clustering structure\n")
} else {
  cat("✗ POOR: There is no clear clustering structure\n")
}

# Dunn's Index
dunn_index <- dunn(distance = dist(data_numeric), clusters = clusters)
cat("\n*** Dunn's Index:", round(dunn_index, 3), "***\n")
cat("(Higher values indicate better separation)\n\n")

### Sida phenolics----

data <- read.csv("sidas_phenolics_data_normalized.csv", 
                 row.names = 1, 
                 check.names = FALSE)

# labels 
labels <- data$Label
cat("\nfind groupss:\n")
print(table(labels))

# keep just numeric data
data_numeric <- data[, -1]  

# Convert all to numeric 
data_numeric <- as.data.frame(lapply(data_numeric, as.numeric))
rownames(data_numeric) <- rownames(data)

# verify dimentions
cat("\nData dimentions:\n")
cat("Samples:", nrow(data_t), "\n")
cat("Metabolites:", ncol(data_t), "\n")

# ===== 2. COPHENETIC COEFFICIENT =====
cat("===== COPHENETIC COEFFICIENT =====\n")

metodos <- c("ward.D2")
resultados_cof <- data.frame(
  Metodo = metodos, 
  Coefenetico = NA,
  stringsAsFactors = FALSE
)

for(i in 1:length(metodos)){
  hc <- hclust(dist(data_numeric), method = metodos[i])
  coph <- cophenetic(hc)
  resultados_cof$Coefenetico[i] <- cor(dist(data_numeric), coph)
}

print(resultados_cof)

# Suggest optimal number based on silhouette
sil_scores <- sapply(2:10, function(k) {
  hc <- hclust(dist(data_numeric), method = metodos)
  clusters <- cutree(hc, k = k)
  mean(silhouette(clusters, dist(data_numeric))[, 3])
})
k_optimo <- which.max(sil_scores) + 1

cat("\n*** SUGGESTED OPTIMAL NUMBER:", k_optimo, "clusters ***\n")
cat("(Based on average silhouette)\n\n")

# ===== 5. VALIDATION INDICES =====
cat("===== VALIDATION INDICES FOR k =", k_optimo, "=====\n")

hc_final <- hclust(dist(data_numeric), method = metodos)
clusters <- cutree(hc_final, k = k_optimo)

# silhouette
sil <- silhouette(clusters, dist(data_numeric))
sil_promedio <- mean(sil[, 3])

cat("\n*** SILHOUETTE:", round(sil_promedio, 3), "***\n")
if(sil_promedio > 0.70) {
  cat("✓ EXCELLENT: Very strong clustering structure\n")
} else if(sil_promedio > 0.50) {
  cat("✓ GOOD: Clear clustering structure\n")
} else if(sil_promedio > 0.25) {
  cat("⚠ REGULAR: Weak clustering structure\n")
} else {
  cat("✗ POOR: There is no clear clustering structure\n")
}

# Dunn's Index
dunn_index <- dunn(distance = dist(data_numeric), clusters = clusters)
cat("\n*** Dunn's Index:", round(dunn_index, 3), "***\n")
cat("(Higher values indicate better separation)\n\n")

##### References --------

pkgs<- c("factoextra","pvclust","dendextend", "cluster","pheatmap","clValid")

# obtain all references
bib_list <- lapply(pkgs, citation)
bibtex_all <- unlist(lapply(bib_list, toBibtex))

# Save as .bib
writeLines(bibtex_all, "referencias.bib")

