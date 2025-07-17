install.packages(c("readxl","dplyr","purrr"))
library(readxl)
library(dplyr)
library(purrr)

cartella_cgm <- "Shanghai_T2DM"
file_summary <- "Shanghai_T2DM_Summary.xlsx"

# 1) Lettura del summary e forzatura di HbA1c come numeric
summary_df <- read_excel(file_summary) %>%
  rename(
    PatientID = `Patient Number`,
    HbA1c     = `HbA1c (mmol/mol)`
  ) %>%
  mutate(HbA1c = as.numeric(HbA1c)) %>%
  select(PatientID, HbA1c)

# 2) Lettura + parsing dei file CGM
files    <- list.files(path = cartella_cgm, pattern = "\\.xls[x]?$", full.names = TRUE)
lista_pazienti <- list()
errori   <- c()
vuoti    <- c()

for (file in files) {
  id_paz <- tools::file_path_sans_ext(basename(file))
  
  tryCatch({
    df_raw <- read_excel(file)
    
    if (nrow(df_raw)==0 || ncol(df_raw)==0) {
      vuoti <- c(vuoti, id_paz)
      next
    }
    
    # 2a) rinomina colonna CGM → "glicemia"
    cgm_col <- grep("^CGM", names(df_raw), value = TRUE)
    if (length(cgm_col)>0) {
      names(df_raw)[names(df_raw)==cgm_col[1]] <- "glicemia"
    } else {
      stop("Colonna CGM non trovata")
    }
    
    date_col <- names(df_raw)[1]
    df_parsed <- df_raw %>%
      transmute(
        PatientID = id_paz,
        Timestamp = as.POSIXct(.data[[date_col]],
                               format = "%d/%m/%Y %H:%M", tz = ""),
        glicemia  = as.numeric(sub(",", ".", glicemia))
      ) %>%
      arrange(Timestamp)
    
    lista_pazienti[[id_paz]] <- df_parsed
    
  }, error = function(e) {
    errori <<- c(errori, id_paz)
  })
}

cat("File letti: ", length(lista_pazienti), "\n")
cat("File vuoti: ", length(vuoti), "->", if(length(vuoti)) paste(vuoti, collapse=", "), "\n")
cat("Errori:     ", length(errori), "->", if(length(errori)) paste(errori, collapse=", "), "\n")

# 3) Estrazione feature di variabilità glicemica
features_df <- map_df(lista_pazienti, function(df) {
  g <- df$glicemia
  tibble(
    PatientID   = df$PatientID[1],
    mean_gl     = mean(g, na.rm = TRUE),
    var_gl      = var(g,  na.rm = TRUE),
    sd_gl       = sd(g,   na.rm = TRUE),
    TIR_70_180  = mean(g >= 70 & g <= 180, na.rm = TRUE),
    TAR_over180 = mean(g > 180,            na.rm = TRUE)
  )
})

# 4) Fusione features + HbA1c
data_tot <- features_df %>%
  left_join(summary_df, by = "PatientID")

library(dplyr)

n_na <- sum(is.na(data_tot$HbA1c))
cat("Totale pazienti senza HbA1c:", n_na, "\n")

# Rimuovo i NA e ricavo un nuovo data set
data_clean <- data_tot %>%
  filter(!is.na(HbA1c))

cat("Pazienti usati ora:", nrow(data_clean), "\n")

# Ricalcolo le correlazioni su data_clean
corr_results_clean <- data_clean %>%
  summarise(
    cor_mean  = cor(mean_gl,    HbA1c, use="complete.obs"),
    cor_sd    = cor(sd_gl,      HbA1c, use="complete.obs"),
    cor_TIR   = cor(TIR_70_180, HbA1c, use="complete.obs"),
    cor_TAR   = cor(TAR_over180,HbA1c, use="complete.obs")
  )

print(corr_results_clean)
# 5) Test di significatività delle correlazioni (cor.test)

cat("\n=== Test di significatività delle correlazioni con HbA1c ===\n")

test_mean <- cor.test(data_clean$mean_gl,    data_clean$HbA1c)
test_sd   <- cor.test(data_clean$sd_gl,      data_clean$HbA1c)
test_TIR  <- cor.test(data_clean$TIR_70_180, data_clean$HbA1c)
test_TAR  <- cor.test(data_clean$TAR_over180,data_clean$HbA1c)

# Stampa i risultati principali
print(list(
  cor_mean  = list(estimate = test_mean$estimate,  p.value = test_mean$p.value,  conf.int = test_mean$conf.int),
  cor_sd    = list(estimate = test_sd$estimate,    p.value = test_sd$p.value,    conf.int = test_sd$conf.int),
  cor_TIR   = list(estimate = test_TIR$estimate,   p.value = test_TIR$p.value,   conf.int = test_TIR$conf.int),
  cor_TAR   = list(estimate = test_TAR$estimate,   p.value = test_TAR$p.value,   conf.int = test_TAR$conf.int)
))

library(ggplot2)
# Scatterplot: glicemia media vs HbA1c
ggplot(data_clean, aes(x = mean_gl, y = HbA1c)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "HbA1c vs Glicemia Media", x = "Glicemia media (mg/dL)", y = "HbA1c (mmol/mol)")


# Scatterplot: deviazione standard vs HbA1c
ggplot(data_clean, aes(x = sd_gl, y = HbA1c)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen") +
  labs(title = "HbA1c vs Deviazione Standard Glicemia", x = "Deviazione Standard", y = "HbA1c (mmol/mol)")

# Scatterplot: TIR vs HbA1c
ggplot(data_clean, aes(x = TIR_70_180, y = HbA1c)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "purple") +
  labs(title = "HbA1c vs Tempo nel Range 70-180", x = "TIR (proporzione)", y = "HbA1c (mmol/mol)")

# Scatterplot: TAR vs HbA1c
ggplot(data_clean, aes(x = TAR_over180, y = HbA1c)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "HbA1c vs Tempo sopra 180", x = "TAR (proporzione)", y = "HbA1c (mmol/mol)")

# 1.1) Carica i pacchetti necessari
install.packages(c("lmtest","car"))
library(lmtest)
library(car)

# 1) Stimo il modello
mod_multi <- lm(HbA1c ~ mean_gl + sd_gl + TIR_70_180 + TAR_over180,
                data = data_clean)
# 2) Controllo output
summary(mod_multi)

library(car)
vif(mod_multi)
#-------------------------------------------------------------------------------
mod_noTAR <- lm(HbA1c ~ mean_gl + sd_gl + TIR_70_180, data = data_clean)
car::vif(mod_noTAR)

summary(mod_noTAR)
#-------------------------------------------------------------------------------
mod_noMean <- lm(HbA1c ~ sd_gl + TIR_70_180, data = data_clean)
car::vif(mod_noMean)

summary(mod_noMean)
#-------------------------------------------------------------------------------
mod_noTIR <- lm(HbA1c ~ mean_gl + sd_gl, data = data_clean)
car::vif(mod_noTIR)
summary(mod_noTIR)
#-------------------------------------------------------------------------------

par(mfrow = c(2,2))
plot(mod_noTAR)

res  <- residuals(mod_noTAR)      
sh   <- shapiro.test(res)
sh$p.value

res  <- residuals(mod_noMean)
sh   <- shapiro.test(res)
sh$p.value
