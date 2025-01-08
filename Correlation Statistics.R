ASV <- read.csv("~/ASVtable_rhizosbulk_forR_noCorM.csv", header=TRUE)

metadata <- read.csv("~/Lily_Rhizo_Metadata_Correlation.csv")%>%
  select(Sample, ph)

merged_data <- merge(ASV, metadata, by = "Sample")


merged_data[, -1] <- lapply(merged_data[, -1], as.numeric)

str(merged_data)

# Perform Spearman's rank correlation for each ASV against pH
cor_results <- lapply(names(merged_data)[2:ncol(merged_data)], function(asv_col) {
  cor.test(merged_data[[asv_col]], merged_data$ph, method = "spearman")
})

rho_values <- unlist(lapply(cor_results, function(x) x$estimate))

weighted_avg_rho <- weighted.mean(abs(rho_values), weights = abs(rho_values), na.rm = TRUE)

median_rho <- median(rho_values, na.rm = TRUE)

print(weighted_avg_rho)
print(rho_values)

