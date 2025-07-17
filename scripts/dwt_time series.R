
pkg <- c("readxl","zoo","wavelets","cluster","stats","ggplot2","lubridate")
inst <- pkg[!pkg %in% installed.packages()[,"Package"]]
if(length(inst)) install.packages(inst)
lapply(pkg, library, character.only = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 1) Lettura dei file e normalizzazione
cartella <- "Shanghai_T2DM"
files    <- list.files(cartella, pattern="\\.xls[x]?$", full.names=TRUE)

lista_pazienti <- list()
for(f in files) {
  id <- tools::file_path_sans_ext(basename(f))
  df <- read_excel(f)
  if(nrow(df) < 2) next
  
  # rinomina e converte CGM → glicemia
  names(df)[grepl("^CGM", names(df))] <- "glicemia"
  df$glicemia <- as.numeric(df$glicemia)
  
  # mantieni Date come character per il parsing flessibile
  df$Date <- as.character(df$Date)
  
  lista_pazienti[[id]] <- df
}

# ──────────────────────────────────────────────────────────────────────────────
# 2) Costruzione di zoo regolari a 15 minuti
ts_list <- lapply(lista_pazienti, function(df) {
  # 2.0 filtra righe incomplete
  df <- df[!is.na(df$Date) & !is.na(df$glicemia), ]
  if(nrow(df) < 2) return(NULL)
  
  # 2.1 parsing ISO con as.POSIXct (formato "YYYY-MM-DD HH:MM:SS")
  times <- as.POSIXct(df$Date, format = "%Y-%m-%d %H:%M:%S", tz = "")
  
  # fallback su dmy_hm per i pochi NA rimasti (es. "DD/MM/YYYY HH:MM")
  idx_na <- which(is.na(times))
  if(length(idx_na) > 0) {
    times[idx_na] <- suppressWarnings(
      dmy_hm(df$Date[idx_na], tz = "")
    )
  }
  
  # rimuovo definitivamente i NA
  valid <- !is.na(times)
  times <- times[valid]
  vals  <- df$glicemia[valid]
  if(length(times) < 2) return(NULL)
  
  # 2.2 ordina
  ord   <- order(times)
  times <- times[ord]
  vals  <- vals[ord]
  
  # 2.3 griglia regolare a 15 minuti
  grid  <- seq(min(times), max(times), by = "15 min")
  
  # 2.4 crea e interpola lo zoo
  z_reg <- na.approx(zoo(vals, times), xout = grid, na.rm = FALSE)
  return(z_reg)
})

# scarta i null (pazienti con dati insufficienti)
ts_list <- ts_list[!sapply(ts_list, is.null)]

# ──────────────────────────────────────────────────────────────────────────────
# 3) DWT di ciascuna serie CGM
filter_type <- "la8"
J_max       <- 4
dwt_list    <- lapply(ts_list, function(z) {
  x     <- coredata(z)
  thisJ <- min(J_max, floor(log2(length(x))))
  dwt(x, filter = filter_type, n.levels = thisJ, boundary = "reflection")
})

# ──────────────────────────────────────────────────────────────────────────────
# 4) Feature-matrix (media e varianza di D1…DJ e A_J)
J <- max(sapply(dwt_list, function(d) length(d@W)))
feature_matrix2 <- t(sapply(dwt_list, function(d) {
  # dettaglio
  stats_det <- unlist(lapply(d@W, function(w) {
    c(mean(w, na.rm=TRUE), var(w, na.rm=TRUE))
  }))
  if(length(stats_det) < 2*J) {
    stats_det <- c(stats_det, rep(NA, 2*J - length(stats_det)))
  }
  # approssimazione
  aJ        <- d@V[[length(d@V)]]
  stats_app <- c(mean(aJ, na.rm=TRUE), var(aJ, na.rm=TRUE))
  c(stats_det, stats_app)
}))
colnames(feature_matrix2) <- c(
  as.vector(t(outer(1:J, c("M_D","V_D"), paste0))),
  "M_A","V_A"
)
rownames(feature_matrix2) <- names(dwt_list)

# ──────────────────────────────────────────────────────────────────────────────
# 5) Standardizzazione delle feature
X <- scale(feature_matrix2)   # media 0, sd 1 per ciascuna colonna

silhouette_avg <- function(dist_matrix) {
  sapply(2:10, function(k) {
    pam_res <- pam(dist_matrix, k = k)
    mean(silhouette(pam_res$clustering, dist_matrix)[, 3])
  })
}

# Distances
dist_euc <- dist(X, method="euclidean")
dist_man <- dist(X, method="manhattan")
corr_mat <- cor(t(X))
dist_cor <- as.dist(1 - corr_mat)

# Silhouette per ciascuna metrica
sil_euc_k <- silhouette_avg(dist_euc)
sil_man_k <- silhouette_avg(dist_man)
sil_cor_k <- silhouette_avg(dist_cor)

# Best k
best_k_euc <- which.max(sil_euc_k) + 1
best_k_man <- which.max(sil_man_k) + 1
best_k_cor <- which.max(sil_cor_k) + 1

cat("Best k (Euclidea):", best_k_euc, "→ Sil =", max(sil_euc_k), "\n")
cat("Best k (Manhattan):", best_k_man, "→ Sil =", max(sil_man_k), "\n")
cat("Best k (Correlazione):", best_k_cor, "→ Sil =", max(sil_cor_k), "\n")

# ──────────────────────────────────────────────────────────────────────────────
# 7) Clustering PAM e Gerarchico su 1–Correlazione
dist_to_use <- dist_cor
k <- best_k_cor

# PAM
pam_res <- pam(dist_to_use, k = k)
clusters_pam <- pam_res$clustering
cat("\nCluster assegnati (PAM):\n")
print(table(clusters_pam))

# Gerarchico
hc_res <- hclust(dist_to_use, method = "ward.D2")
plot(hc_res, main = paste("Dendrogramma (distanza selezionata)"))
clusters_hc <- cutree(hc_res, k = k)
cat("\nCluster assegnati (Gerarchico):\n")
print(table(clusters_hc))

# Silhouette finale
plot(silhouette(clusters_pam, dist_to_use), main="Silhouette PAM (distanza selezionata)")


library(readxl)

# Carica il file con i dati clinici
summary <- read_excel("Shanghai_T2DM_Summary.xlsx")

# Aggiungi i cluster al dataframe
summary$Cluster <- clusters_hc[as.character(summary$`Patient Number`)]

# Controlla che non ci siano NA
sum(is.na(summary$Cluster))

# Converto le variabili in fattori
summary$Hypoglycemia <- as.factor(summary$`Hypoglycemia (yes/no)`)
summary$Drinker <- as.factor(summary$`Alcohol Drinking History (drinker/non-drinker)`)
summary$Cluster <- as.factor(summary$Cluster)

tab_hypo <- table(summary$Cluster, summary$Hypoglycemia)
print(tab_hypo)
chisq.test(tab_hypo)

tab_drink <- table(summary$Cluster, summary$Drinker)
print(tab_drink)
chisq.test(tab_drink)

summary$GA <- summary$`Glycated Albumin (%)`

# Controlla la presenza di NA
summary$Cluster <- as.factor(summary$Cluster)
summary$GA <- as.numeric(summary$`Glycated Albumin (%)`)
summary_clean <- summary[!is.na(summary$GA), ]

by(summary_clean$GA, summary_clean$Cluster, shapiro.test)

kruskal.test(GA ~ Cluster, data = summary_clean)

install.packages("FSA")
library(FSA)

# Test di Dunn con correzione Bonferroni
dunnTest(GA ~ Cluster, data = summary_clean, method = "bonferroni")

ggplot(summary_clean, aes(x = Cluster, y = GA)) +
  geom_boxplot(aes(fill = Cluster)) +
  labs(title = "Glycated Albumin per Cluster") +
  theme_minimal()

names(summary)

summary$HbA1c_real <- as.numeric(summary$`HbA1c (mmol/mol)`)
summary_hba1c <- summary[!is.na(summary$HbA1c_real), ]

by(summary_hba1c$HbA1c_real, summary_hba1c$Cluster, shapiro.test)

kruskal.test(HbA1c_real ~ Cluster, data = summary_hba1c)

dunnTest(HbA1c_real ~ Cluster, data = summary_hba1c, method = "bonferroni")
