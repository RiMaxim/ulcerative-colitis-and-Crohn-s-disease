# PERMANOVA: объясняет ли генотип общий состав микробиома (роды, CLR; расстояние Эйчисона)
# с поправкой на диагноз, возраст и пол. Запуск: Rscript genotype_permanova.R
suppressPackageStartupMessages(library(vegan))
d <- read.csv("results/genotype_omics_input.csv", check.names = FALSE, row.names = 1)
meta_cols <- c("diagnosis", "age", "sex_male", "FUT2_nonsecretor", "NOD2_carrier", "ATG16L1_T300A", "Shannon")
X <- as.matrix(d[, setdiff(names(d), meta_cols)])
set.seed(1)
out <- do.call(rbind, lapply(c("FUT2_nonsecretor", "NOD2_carrier", "ATG16L1_T300A"), function(g) {
  ok <- complete.cases(d[, c("diagnosis", "age", "sex_male", g)])
  m <- d[ok, ]; D <- dist(X[ok, ])
  a <- adonis2(as.formula(paste("D ~ diagnosis + age + sex_male +", g)), data = m, by = "margin", permutations = 9999)
  data.frame(genotype = g, n = sum(ok), R2 = a[g, "R2"], F = a[g, "F"], p = a[g, "Pr(>F)"],
             R2_diagnosis = a["diagnosis", "R2"], p_diagnosis = a["diagnosis", "Pr(>F)"])
}))
print(out, digits = 3)
write.csv(out, "results/genotype_permanova.csv", row.names = FALSE)
