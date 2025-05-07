#libs
library(readxl)
library(dplyr)
library(edgeR)
library(MetaCycle)
library(tidyr)
library(ggplot2)

#file load
counts <- read_excel("GSE196728_Peru_raw.xlsx")
colnames(counts)[1] <- "Gene"
rownames(counts) <- counts$Gene
counts <- counts[, -1]

#metadata
sample_names <- colnames(counts)
metadata <- data.frame(
  Sample = sample_names,
  ID = sub("_.*", "", sample_names),
  Altitude = sub(".*_(SL|PU|RI).*", "\\1", sample_names),
  Time = as.numeric(sub(".*_(SL|PU|RI)(\\d{1,2})\\.raw", "\\2", sample_names)),
  Donor = sub("_.*", "", sample_names),
  stringsAsFactors = FALSE
)

#shifts 2:00 to 26 to reflect end of circadian cycle
metadata$Time <- ifelse(metadata$Time == 2, 26, metadata$Time)
counts <- counts[, metadata$Sample]

#edgeR
dge <- DGEList(counts = counts)
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)

#MC input
time_vector <- metadata$Time
logCPM_with_time <- rbind(time_vector, logCPM)
rownames(logCPM_with_time)[1] <- ""
write.table(logCPM_with_time, file = "expr_matrix_line1.txt", sep = "\t", quote = FALSE, col.names = FALSE)

#MC
meta2d(
  infile = "expr_matrix_line1.txt",
  outdir = "MetaCycle_out",
  filestyle = "txt",
  timepoints = "line1",
  cycMethod = c("LS", "ARS"),
  minper = 20,
  maxper = 28,
  outputFile = TRUE,
  adjustPhase = "predictedPer"
)

#top rhithmic genes
meta_result <- read.table("MetaCycle_out/meta2d_expr_matrix_line1.txt", header = TRUE, sep = "\t", check.names = FALSE)
top3_genes <- as.character(meta_result$CycID[order(meta_result$LS_BH.Q)][1:3])
expr_subset <- logCPM[top3_genes, ]
expr_df <- as.data.frame(t(expr_subset))

#metadata
expr_df$Time <- metadata$Time
expr_df$Altitude <- recode(metadata$Altitude, "SL" = "Sea Level (200 m)", "PU" = "Puno (3,800 m)", "RI" = "Rinconada (5,100 m)")
expr_df$Sample <- metadata$Sample

#plot
long_df <- expr_df %>%
  pivot_longer(cols = all_of(top3_genes), names_to = "Gene", values_to = "Expression") %>%
  group_by(Altitude, Gene, Time) %>%
  summarise(MeanExpr = mean(Expression), .groups = "drop")

for (alt in unique(long_df$Altitude)) {
  df_alt <- filter(long_df, Altitude == alt)
  
  p <- ggplot(df_alt, aes(x = Time, y = MeanExpr, color = Gene)) +
    geom_point(size = 2) +
    geom_smooth(method = "loess", span = 0.75, se = FALSE, linetype = "dashed") +
    scale_x_continuous(
      breaks = c(6, 10, 14, 18, 22, 26),
      labels = c("06:00", "10:00", "14:00", "18:00", "22:00", "02:00")
    ) +
    labs(title = paste("Top 3 Oscillating Genes –", alt),
         x = "Time (hours)",
         y = "Mean log2 CPM Expression",
         color = "Gene") +
    theme_minimal(base_size = 14)
  
  ggsave(paste0("GSE196728_top3_", gsub(" ", "_", alt), ".png"),
         plot = p, width = 8, height = 5.5, dpi = 300)
}
