library(readxl)
library(dplyr)
library(ggplot2)

# Caricamento
summary <- read_excel("C:/Users/rossa/library/magistrale/I ANNO I SEMESTRE/Statistica e Analisi dei Dati/datasets/Shanghai_T2DM_Summary.xlsx")

# Definizione di valori testuali che devono essere trattati come NA (es. "/", "none", ecc.)
valori_na <- c("/", "none", "None", "-", "NA", "na", "N/A", "", " ", "—")

# Applica la pulizia: ogni colonna viene trasformata, e i "falsi NA" vengono convertiti in veri NA
summary[] <- lapply(summary, function(col) {
  col <- ifelse(trimws(tolower(as.character(col))) %in% tolower(valori_na), NA, col)
  return(col)
})

# Converte la colonna "HbA1c" da testo a numero (virgole convertite in punti)
summary$HbA1c <- as.numeric(gsub(",", ".", summary$`HbA1c (mmol/mol)`))

# Converte la variabile 'Alcohol Drinking History' in fattore con due livelli
summary$Alcohol <- factor(summary$`Alcohol Drinking History (drinker/non-drinker)`,
                          levels = c("non-drinker", "drinker"))

# Crea un boxplot dei valori di HbA1c per ciascun gruppo (drinker vs non-drinker)
# Il punto nero rappresenta la mediana
ggplot(summary, aes(x = Alcohol, y = HbA1c, fill = Alcohol)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "point", shape = 20, size = 4, color = "black") +
  labs(title = "HbA1c nei pazienti con e senza consumo di alcol", 
       y = "HbA1c (mmol/mol)", x = "Consumo di alcol") +
  theme_minimal()

# Test di normalità (Shapiro-Wilk) per i valori di HbA1c nel gruppo "drinker"
shapiro.test(summary$HbA1c[summary$Alcohol == "drinker"])

# Test di normalità (Shapiro-Wilk) per i valori di HbA1c nel gruppo "non-drinker"
shapiro.test(summary$HbA1c[summary$Alcohol == "non-drinker"])

# Esecuzione del test di Wilcoxon (non parametrico) per verificare differenze tra i due gruppi
wilcox.test(HbA1c ~ Alcohol, data = summary)

# Calcolo delle statistiche descrittive per ciascun gruppo: mediana, media e deviazione standard
summary %>%
  group_by(Alcohol) %>%
  summarise(
    Mediana = median(HbA1c, na.rm = TRUE),
    Media = mean(HbA1c, na.rm = TRUE),
    SD = sd(HbA1c, na.rm = TRUE)
  )

