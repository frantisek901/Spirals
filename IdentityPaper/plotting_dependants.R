#To be run after script for step3.
ESBG_Set <- ts[ts("ESBG")]
diversity_set <- ts[ts("diversity")]
extremness_set <- ts[ts("extremness")]
dependant_measures <- data.frame(ESBG_Set, diversity_set, extremness_set)

# plotting
png(filename="step3_Diversity_vs_ESBG_allRuns.png")
plot(dependant_measures[,1], dependant_measures[,2],
     main="Diversity vs Polarization for every run",
     xlab="ESBG polarization",
     ylab="Diversity")
dev.off()

png(filename="step3_Extremness_vs_ESBG_allRuns.png")
plot(dependant_measures[,1], dependant_measures[,3],
     main="Extremness vs Polarization for every run",
     xlab="ESBG polarization",
     ylab="Extremness")
dev.off()

png(filename="step3_Extremness_vs_Diversity_allRuns.png")
plot(dependant_measures[,2], dependant_measures[,3],
     main="Extremness vs Diversity for every run",
     xlab="Diversity",
     ylab="Extremness")
dev.off()
