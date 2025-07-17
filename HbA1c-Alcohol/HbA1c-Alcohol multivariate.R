# === LIBRERIE ===
library(readxl)     # Lettura file Excel
library(dplyr)      # Manipolazione dati
library(lmtest)     # Test statistici sui modelli
library(broom)      # Estrazione ordinata risultati dei modelli
library(ggplot2)    # Visualizzazione grafici
library(car)        # Diagnostica avanzata dei modelli

# === 1. CARICAMENTO E PULIZIA DEL DATASET ===
summary <- read_excel("C:/Users/rossa/library/magistrale/I ANNO I SEMESTRE/Statistica e Analisi dei Dati/datasets/Shanghai_T2DM_Summary.xlsx")

# Riconoscimento dei valori mancanti in varie forme
valori_na <- c("/", "none", "None", "-", "NA", "na", "N/A", "", " ", "—")
summary[] <- lapply(summary, function(col) {
  col <- ifelse(trimws(tolower(as.character(col))) %in% tolower(valori_na), NA, col)
  return(col)
})

# === 2. CONVERSIONE VARIABILI ===
summary$HbA1c   <- as.numeric(gsub(",", ".", summary$`HbA1c (mmol/mol)`))
summary$Age     <- as.numeric(gsub(",", ".", summary$`Age (years)`))
summary$BMI     <- as.numeric(gsub(",", ".", summary$`BMI (kg/m2)`))
summary$Smoking <- as.numeric(gsub(",", ".", summary$`Smoking History (pack year)`))

# Variabili categoriali
summary$Gender  <- factor(summary$`Gender (Female=1, Male=2)`, levels = c(1, 2), labels = c("F", "M"))
summary$Alcohol <- factor(summary$`Alcohol Drinking History (drinker/non-drinker)`)
summary$Terapia <- factor(summary$`Hypoglycemic Agents`)

# === 3. RIMOZIONE DEI MISSING E SELEZIONE VARIABILI ===
vars_modello <- c("HbA1c", "Age", "BMI", "Smoking", "Gender", "Alcohol", "Terapia")
sapply(summary[vars_modello], function(x) sum(is.na(x)))  # controllo NA

dati_model <- summary %>%
  select(all_of(vars_modello)) %>%
  na.omit()

# === 4. ESPOLORAZIONE INIZIALE ===
# Correlazioni tra variabili numeriche
cor(dati_model %>% select(Age, BMI, Smoking), use = "complete.obs")

# Istogramma e boxplot di HbA1c
hist(dati_model$HbA1c, main = "Distribuzione HbA1c", col = "skyblue")
boxplot(dati_model$HbA1c, main = "Boxplot HbA1c")

# === 5. MODELLO INIZIALE: REGRESSIONE LINEARE COMPLETA ===
modello <- lm(HbA1c ~ Age + BMI + Smoking + Gender + Alcohol + Terapia, data = dati_model)

# Diagnostica del modello
par(mfrow = c(2, 2))
plot(modello)

# Test diagnostici
bptest(modello)           # Omocedasticità (Breusch-Pagan)
shapiro.test(resid(modello))  # Normalità dei residui
dwtest(modello)           # Autocorrelazione (Durbin-Watson)
plot(modello, which = 4)  # Residui vs leverage (Cook’s distance)

# === 6. GESTIONE OSSERVAZIONI INFLUENTI ===
influential <- which(cooks.distance(modello) > 4 / length(modello$residuals))
dati_clean <- dati_model[-influential, ]

# Nuovo modello con dati puliti
modello_sens <- lm(HbA1c ~ Age + BMI + Smoking + Gender + Alcohol + Terapia, data = dati_clean)
summary(modello_sens)

# === 7. VISUALIZZAZIONE COEFFICIENTI ===
tidy(modello) %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
  coord_flip() +
  labs(title = "Effetti stimati sulla HbA1c", x = "Variabile", y = "Stima")

# === 8. RICATEGORIZZAZIONE TERAPIA ===
dati_clean$Terapia_cat <- case_when(
  grepl("^metformin$", dati_clean$Terapia, ignore.case = TRUE) ~ "Solo Metformin",
  grepl("insulin", dati_clean$Terapia, ignore.case = TRUE) &
    grepl("metformin|acarbose|glimepiride|sitagliptin|voglibose|gliclazide|glulisine", dati_clean$Terapia, ignore.case = TRUE) ~ "Insulina + Orali",
  grepl("insulin", dati_clean$Terapia, ignore.case = TRUE) ~ "Solo Insulina",
  TRUE ~ "Orali"
)
dati_clean$Terapia_cat <- factor(dati_clean$Terapia_cat)
table(dati_clean$Terapia_cat)

# === 9. MODELLO SEMPLIFICATO CON TERAPIA_CAT ===
mod_simpl <- lm(HbA1c ~ Age + BMI + Smoking + Gender + Alcohol + Terapia_cat, data = dati_clean)
summary(mod_simpl)

# Diagnostica del nuovo modello semplificato
par(mfrow = c(2, 2))
plot(mod_simpl)

# Boxplot per tipo di terapia
ggplot(dati_clean, aes(x = Terapia_cat, y = HbA1c, fill = Terapia_cat)) +
  geom_boxplot() +
  labs(title = "HbA1c per categoria terapeutica", y = "HbA1c (mmol/mol)", x = "Categoria terapia") +
  theme_minimal()

# === 10. ANALISI DELL’EFFETTO DELL’ALCOL ===
mod_alcohol <- lm(HbA1c ~ Alcohol, data = dati_clean)
summary(mod_alcohol)

mod_completo <- lm(HbA1c ~ Age + BMI + Smoking + Gender + Alcohol + Terapia_cat, data = dati_clean)
summary(mod_completo)

# Confronto grafico tra modelli
bind_rows(
  tidy(mod_alcohol, conf.int = TRUE) %>% filter(term == "Alcoholnon-drinker") %>% mutate(Modello = "Solo Alcohol"),
  tidy(mod_completo, conf.int = TRUE) %>% filter(term == "Alcoholnon-drinker") %>% mutate(Modello = "Modello Completo")
) %>%
  ggplot(aes(x = Modello, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(color = "tomato", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Effetto del consumo di alcol su HbA1c", y = "Stima dell'effetto (mmol/mol)") +
  theme_minimal()

# === 11. TRASFORMAZIONE LOGARITMICA DI HbA1c ===
dati_clean$log_HbA1c <- log(dati_clean$HbA1c)

mod_log <- lm(log_HbA1c ~ Age + BMI + Smoking + Gender + Alcohol + Terapia_cat, data = dati_clean)

# Diagnostica del modello logaritmico
par(mfrow = c(2, 2))
plot(mod_log)

# Test di normalità residui del modello log
shapiro.test(resid(mod_log))

# Output del modello log
summary(mod_log)
