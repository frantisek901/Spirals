#To be run after script for step3.

#1. plotting and saving dependant measures of every simulation:

ESBG_Set <- ts[ts("ESBG")]
diversity_set <- ts[ts("diversity")]
extremness_set <- ts[ts("extremness")]
dependant_measures <- data.frame(ESBG_Set, diversity_set, extremness_set)


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


# 2. plotting aggregates for common parameter combinations

ESBG_mean_set <- tm[7]
ESBG_sd_set <- tm[8]

extremness_mean_set <- tm[5]
extremness_sd_set <- tm[6]


diversity_mean_set <- tm[3]
diversity_sd_set <- tm[4]

dependant_measures_aggregate <- data.frame(ESBG_mean_set,
                                           ESBG_sd_set,
                                           diversity_mean_set,
                                           diversity_sd_set,
                                           extremness_mean_set,
                                           extremness_sd_set)

#plotting means


png(filename="step3_Diversity_Mean_vs_ESBG_Mean.png")
plot(dependant_measures_aggregate[,1], dependant_measures_aggregate[,3],
     main="Diversity vs Polarization",
     xlab="ESBG polarization_Mean",
     ylab="Diversity_Mean")
dev.off()

png(filename="step3_Extremness_Mean_vs_ESBG_Mean.png")
plot(dependant_measures_aggregate[,1], dependant_measures_aggregate[,5],
     main="Extremness_Mean vs Polarization_Mean for every run",
     xlab="ESBG_Mean polarization",
     ylab="Extremness_Mean")
dev.off()

png(filename="step3_Extremness_Mean_vs_Diversity_Mean.png")
plot(dependant_measures_aggregate[,3], dependant_measures_aggregate[,5],
     main="Extremness_Mean vs Diversity_Mean for every run",
     xlab="Diversity_Mean",
     ylab="Extremness_Mean")
dev.off()

