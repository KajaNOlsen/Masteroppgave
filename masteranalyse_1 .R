
library(haven) 


# Laster inn datasett og gir det navn masterdata
masterdata <- read_dta("N:/durable/Project_A/project_A4/7 - Tilrettelagte data/Analysefil.dta")

library(WeightedCluster)
library(TraMineR)
library(TraMineRextras)
library(cluster)
library(foreign)
library(ggplot2)
library(nnet)
library(tidyr) 
library(dplyr)

#Statuser som finnes i datasettet
seqstatl(masterdata[, 22:32])

#Definerer alfabet - de faktiske verdiene i datasettet som tekst
analysis.alphabet <- c("0", "1", "2", "3", "4", "5", "6")

analysis.labels <- c("(0%)",
                     "(100%)",
                     "(80-99%)",
                     "(60-79%)",
                     "(40-59%)",
                     "(20-39%)",
                     "(1-19%)")

analysis.scode <- c("IN", "FU", "HO", "MH", "MI", "LA", "SL")

# Definerer variabler / sekvensobjekter
mor_kolonner <- c("Innt_status_mor_2012", "Innt_status_mor_2013",
                  "Innt_status_mor_2014", "Innt_status_mor_2015",
                  "Innt_status_mor_2016")

far_kolonner <- c("Innt_status_far_2012", "Innt_status_far_2013",
                  "Innt_status_far_2014", "Innt_status_far_2015",
                  "Innt_status_far_2016")

# Sekvens mor og far
seq_mor <- seqdef(masterdata,
                  var      = mor_kolonner,
                  alphabet = analysis.alphabet,
                  states   = analysis.scode,
                  labels   = analysis.labels,
                  xtstep   = 1,
                  x1       = 2012)

seq_far <- seqdef(masterdata,
                  var      = far_kolonner,
                  alphabet = analysis.alphabet,
                  states   = analysis.scode,
                  labels   = analysis.labels,
                  xtstep   = 1,
                  x1       = 2012)



agg_mor <- wcAggregateCases(masterdata[, mor_kolonner])
agg_far <- wcAggregateCases(masterdata[, far_kolonner])

print(agg_mor)  # Ser hvor mange unike sekvenser som finnes
print(agg_far)

unique_mor_data <- masterdata[agg_mor$aggIndex, mor_kolonner]
unique_far_data <- masterdata[agg_far$aggIndex, far_kolonner]


uniqueSeq_mor <- seqdef(unique_mor_data,
                        weights  = agg_mor$aggWeights,
                        alphabet = analysis.alphabet,
                        states   = analysis.scode,
                        labels   = analysis.labels,
                        xtstep   = 1,
                        x1       = 2012)

uniqueSeq_far <- seqdef(unique_far_data,
                        weights  = agg_far$aggWeights,
                        alphabet = analysis.alphabet,
                        states   = analysis.scode,
                        labels   = analysis.labels,
                        xtstep   = 1,
                        x1       = 2012)



cost_mor <- seqcost(uniqueSeq_mor, method = "INDELSLOG", weighted = TRUE)
cost_far <- seqcost(uniqueSeq_far, method = "INDELSLOG", weighted = TRUE)

# Beregner avstandsmatriser med Optimal Matching (OM)
dist_mor <- seqdist(uniqueSeq_mor, method = "OM",
                    indel = cost_mor$indel, sm = cost_mor$sm)

dist_far <- seqdist(uniqueSeq_far, method = "OM",
                    indel = cost_far$indel, sm = cost_far$sm)

dist_mor_full <- seqdist(seq_mor, method = "OM",
                         indel = cost_mor$indel, sm = cost_mor$sm)

dist_far_full <- seqdist(seq_far, method = "OM",
                         indel = cost_far$indel, sm = cost_far$sm)


# Ward-klynger 
ward_mor <- hclust(as.dist(dist_mor), method = "ward.D2",
                   members = agg_mor$aggWeights)

ward_far <- hclust(as.dist(dist_far), method = "ward.D2",
                   members = agg_far$aggWeights)


# PAM-klynger 
pam_mor <- wcKMedRange(dist_mor, kvals = 2:8,
                       weights       = agg_mor$aggWeights,
                       initialclust  = ward_mor)

pam_far <- wcKMedRange(dist_far, kvals = 2:8,
                       weights       = agg_far$aggWeights,
                       initialclust  = ward_far)


klynge_mor <- pam_mor$clustering[[which(pam_mor$kvals == 3)]]
klynge_far <- pam_far$clustering[[which(pam_far$kvals == 2)]]



# Plot kvalitetsmal for valv av antall klynger
plot(pam_mor, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore",
     main = "Mor: klyngekvalitet")
plot(pam_far, stat = c("ASWw", "HG", "PBC", "HC"), norm = "zscore",
     main = "Far: klyngekvalitet")


# Se topp 5 løsninger rangert etter kvalitetsmål
summary(pam_mor, max.rank = 5)
summary(pam_far, max.rank = 5)

# Setter farge for plottene - fargene tilsvarer statusene i rekkefølgen som ble satt over
cpal(uniqueSeq_mor) <- c("#969696", "#1a9641", "#a6d96a",
                         "#ffffbf", "#fdae61", "#d7191c", "#7b3294")
cpal(uniqueSeq_far) <- c("#969696", "#1a9641", "#a6d96a",
                         "#ffffbf", "#fdae61", "#d7191c", "#7b3294")

# Tilstandsfordeling over tid per klynge (seqdplot)
par(mar = c(3, 3, 3, 3))

# Sekvensplott mor - alle tre samlet
seqdplot(uniqueSeq_mor,
         group   = pam_mor$clustering$cluster3,
         sortv   = sortv(uniqueSeq_mor, start = "beg"),
         border  = NA, cols = 2,
         cex.axis = 0.5, cex.main = 0.8,
         main = "Mor: inntektsutvikling per klynge",
         ylab = "", xlab = "")

#Sekvensplott mor med riktig navn
par(mar = c(3, 3, 3, 6))

klynger_mor <- pam_mor$clustering$cluster3

#Plot varig-klyngen
seqdplot(uniqueSeq_mor[klynger_mor == 336, ],
         sortv   = sortv(uniqueSeq_mor[klynger_mor == 336, ], start = "beg"),
         border  = NA, cols = 2,
         cex.axis = 0.5, cex.main = 0.8,
         main = "Varig-klyngen",
         ylab = "", xlab = "")

#Plot moderat-klyngen
seqdplot(uniqueSeq_mor[klynger_mor == 527, ],
         sortv   = sortv(uniqueSeq_mor[klynger_mor == 527, ], start = "beg"),
         border  = NA, cols = 2,
         cex.axis = 0.5, cex.main = 0.8,
         main = "Moderat-klyngen",
         ylab = "", xlab = "")

#Plot referanseklyngen
seqdplot(uniqueSeq_mor[klynger_mor == 73, ],
         sortv   = sortv(uniqueSeq_mor[klynger_mor == 73, ], start = "beg"),
         border  = NA, cols = 2,
         cex.axis = 0.5, cex.main = 0.8,
         main = "Referanseklyngen",
         ylab = "", xlab = "")


#Plott for far - begge samlet
seqdplot(uniqueSeq_far,
         group   = pam_far$clustering$cluster2,
         sortv   = sortv(uniqueSeq_far, start = "beg"),
         border  = NA, cols = 2,
         cex.axis = 0.5, cex.main = 0.8,
         main = "Far: inntektsutvikling per klynge",
         ylab = "", xlab = "")

#Fars klynger med riktige navn
klynger_far <- pam_far$clustering$cluster2

#Fars moderat-klynge
seqdplot(uniqueSeq_far[klynger_far == 315, ],
         sortv   = sortv(uniqueSeq_far[klynger_far == 315, ], start = "beg"),
         border  = NA, cols = 2,
         cex.axis = 0.5, cex.main = 0.8,
         main = "Moderat-klyngen",
         ylab = "", xlab = "")

#Referanseklyngen
seqdplot(uniqueSeq_far[klynger_far == 58, ],
         sortv   = sortv(uniqueSeq_far[klynger_far == 58, ], start = "beg"),
         border  = NA, cols = 2,
         cex.axis = 0.5, cex.main = 0.8,
         main = "Referanseklyngen",
         ylab = "", xlab = "")

# Representativ sekvens per klynge (seqmsplot)
seqmsplot(uniqueSeq_mor,
          group   = pam_mor$clustering$cluster3,
          border  = NA, cols = 2,
          cex.axis = 0.5, cex.main = 0.8,
          main = "Mor: representativ sekvens per klynge",
          ylab = "", xlab = "")

seqmsplot(uniqueSeq_far,
          group   = pam_far$clustering$cluster2,
          border  = NA, cols = 2,
          cex.axis = 0.5, cex.main = 0.8,
          main = "Far: representativ sekvens per klynge",
          ylab = "", xlab = "")


# Gjennomsnittlig tid i hver tilstand per klynge
seqmtplot(uniqueSeq_mor,
          group   = pam_mor$clustering$cluster3,
          border  = NA, cols = 2,
          cex.axis = 0.5, cex.main = 0.8,
          main = "Mor: gjennomsnittlig tid per tilstand",
          ylab = "", xlab = "")

# Indeksplott - alle individuelle sekvenser
seqiplot(uniqueSeq_mor,
         group   = pam_mor$clustering$cluster3,
         tlim    = 0, border = NA, cols = 2,
         cex.axis = 0.5, cex.main = 0.8,
         main = "Mor: alle sekvenser per klynge",
         ylab = "", xlab = "")

# Klynger kobles tilbake til hele datasettet
klynge_mor <- pam_mor$clustering$cluster3
klynge_far <- pam_far$clustering$cluster2

masterdata$klynge_mor <- klynge_mor[agg_mor$disaggIndex]
masterdata$klynge_far <- klynge_far[agg_far$disaggIndex]


# Sjekk fordelingen av klynger etter aggregering
table(masterdata$klynge_mor)

table(masterdata$klynge_far)


# Konverterer til faktor med riktig referanse
masterdata$klynge_mor <- relevel(as.factor(masterdata$klynge_mor), ref = "73")
masterdata$klynge_far <- relevel(as.factor(masterdata$klynge_far), ref = "58")


# Forbereder variabler til regresjon
tapply(rowMeans(masterdata[, mor_kolonner], na.rm = TRUE),
       masterdata$klynge_mor, mean)

# Referansekategorien settes til klyngen med høyest inntekt
masterdata$klynge_mor <- relevel(masterdata$klynge_mor, ref = "73")
masterdata$klynge_far <- relevel(masterdata$klynge_far, ref = "58")

#############################################################################
#Gir verdiene navn
 
# Psykososial belastning
masterdata$forskjell_strain_5 <- factor(
  masterdata$forskjell_strain_5,
  levels = 1:5,
  labels = c("(Q1)", "(Q2)", "(Q3)",
             "(Q4)", "(Q5)"))

masterdata$forskjell_strain_5 <- relevel(
  masterdata$forskjell_strain_5, ref = "(Q2)")

# Mekanisk belastning
masterdata$forskjell_mekanisk_5 <- factor(
  masterdata$forskjell_mekanisk_5,
  levels = 1:5,
  labels = c("(Q1)", "(Q2)", "(Q3)",
             "(Q4)", "(Q5)"))

masterdata$forskjell_mekanisk_5 <- relevel(
  masterdata$forskjell_mekanisk_5, ref = "(Q3)")


# Kryssbelastning
masterdata$forskjell_kryssbelastning_5 <- factor(
  masterdata$forskjell_kryssbelastning_5,
  levels = 1:5,
  labels = c("(Q1)", "(Q2)", "(Q3)",
             "(Q4)", "(Q5)"))

masterdata$forskjell_kryssbelastning_5 <- relevel(
  masterdata$forskjell_kryssbelastning_5, ref = "(Q3)")

# Inntektsvariabel
masterdata$forskjell_inntekt_5 <- factor(
  masterdata$forskjell_inntekt_5,
  levels = 1:5,
  labels = c("(Q1)", "(Q2)", "(Q3)",
             "(Q4)", "(Q5)")
)

masterdata$forskjell_inntekt_5 <- relevel(
  masterdata$forskjell_inntekt_5, ref = "(Q4)")


##################################################################################

library(broom)

# SELVE REGRESJONENE - multinom logistisk regresjon

# Hjelpefunksjon for p-verdier som brukes for alle modellene
pval <- function(mod) {
  z <- summary(mod)$coefficients / summary(mod)$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  return(round(p, 3))}

# MODELL M1 uten kontrollvariabler

mod1 <- multinom(klynge_mor ~ klynge_far,
                 data = masterdata)
summary(mod1)
cat("P-verdier modell 1:\n"); print(pval(mod1))

# Krysstabell med enkel beskrivelse av kj??nnsforskjeller

krysstabell_m1 <- table(masterdata$klynge_far, masterdata$klynge_mor)
prop.table(krysstabell_m1, margin = 1)
chisq.test(krysstabell_m1)

# MODELL M2

mod2 <- multinom(klynge_mor ~ forskjell_inntekt_5,
                 data = masterdata)
summary(mod2)
cat("P-verdier modell 2:\n"); print(pval(mod2))

# Krysstabell
krysstabell_m2 <- table(masterdata$forskjell_inntekt_5, masterdata$klynge_mor)
prop.table(krysstabell_m2, margin = 1)
chisq.test(krysstabell_m2)

# Modell M3

mod3 <- multinom(klynge_mor ~ forskjell_strain_5,
                 data = masterdata)
summary(mod3)
cat("P-verdier modell 3:\n"); print(pval(mod3))


# MODELL M4

mod4 <- multinom(klynge_mor ~ forskjell_mekanisk_5,
                 data = masterdata)
summary(mod4)
cat("P-verdier modell 4:\n"); print(pval(mod4))

# MODELL M5

mod5 <- multinom(klynge_mor ~ forskjell_kryssbelastning_5,
                 data = masterdata)
summary(mod5)
cat("P-verdier modell 5:\n"); print(pval(mod5))

### MODELL M6

mod6 <- multinom(klynge_far ~ forskjell_kryssbelastning_5,
                 data = masterdata, trace = FALSE)
summary(mod6)
cat("P-verdier modell 6:\n"); print(pval(mod6))

######################################################################################
# FULLSTENDIG MODELLSTRUKTUR MED ALLE KONTROLLVARIABLER

library(nnet)
library(dplyr)

masterdata <- masterdata %>%
  mutate(
    klynge_mor = relevel(as.factor(klynge_mor), ref = "73"),
    klynge_far = relevel(as.factor(klynge_far), ref = "58"),
    forskjell_kryssbelastning_5 = relevel(as.factor(forskjell_kryssbelastning_5), ref = "(Q3)"),
    forskjell_inntekt_5 = relevel(as.factor(forskjell_inntekt_5), ref = "(Q4)"),
    innvandrer_mor = as.factor(ifelse(invkat_mor == "A", 0, 1)),
    innvandrer_far = as.factor(ifelse(invkat_far == "A", 0, 1)),
    sykefravar_mor = as.factor(sykefravar_mor),
    sykefravar_far = as.factor(sykefravar_far),
    yrkinntekt_mor_2009_t = yrkinntekt_mor_2009 / 1000,
    yrkinntekt_far_2009_t = yrkinntekt_far_2009 / 1000,
    forskjell_strain_5 = relevel(as.factor(forskjell_strain_5), ref = "(Q2)"),
    forskjell_mekanisk_5 = relevel(as.factor(forskjell_mekanisk_5), ref = "(Q3)"))

# Under lager jeg to ulike sett med kontrollvariabler
# Kontrollvariabler for modell M1 og M2
kontroller_full <- "+ utd_nivaa_mor + utd_nivaa_far + 
                    innvandrer_mor + innvandrer_far + 
                    sykefravar_mor + sykefravar_far + 
                    yrkinntekt_mor_2009_t + yrkinntekt_far_2009_t + 
                    mekanisk_yrke_indeks_mor + mekanisk_yrke_indeks_far + 
                    strain_yrke_indeks_mor + strain_yrke_indeks_far"

#Kontrollvariabler for modell M3 til M8
kontroller_uten_belastning <- "+ utd_nivaa_mor + utd_nivaa_far + 
                                innvandrer_mor + innvandrer_far + 
                                sykefravar_mor + sykefravar_far + 
                                yrkinntekt_mor_2009_t + yrkinntekt_far_2009_t"

################################################################################
#MODELLENE MED KONTROLLVARIABLER

#Modell M1
mod1k <- multinom(as.formula(paste("klynge_mor ~ klynge_far", kontroller_full)),
                  data = masterdata, trace = FALSE)

#Modell M2 
mod2k <- multinom(as.formula(paste("klynge_mor ~ forskjell_inntekt_5", kontroller_full)),
                  data = masterdata, trace = FALSE)

#Modell M3 
mod3k <- multinom(as.formula(paste("klynge_mor ~ forskjell_strain_5", kontroller_uten_belastning)),
                  data = masterdata, trace = FALSE)

#Modell M4
mod4k <- multinom(as.formula(paste("klynge_mor ~ forskjell_mekanisk_5", kontroller_uten_belastning)),
                  data = masterdata, trace = FALSE)

#Modell M5
mod5k <- multinom(as.formula(paste("klynge_mor ~ forskjell_kryssbelastning_5", kontroller_uten_belastning)),
                  data = masterdata, trace = FALSE)

#Modell M6
mod6k <- multinom(as.formula(paste("klynge_far ~ forskjell_kryssbelastning_5", kontroller_uten_belastning)),
                  data = masterdata, trace = FALSE)


#Modell M7
mod7k <- multinom(as.formula(paste("klynge_mor ~ forskjell_kryssbelastning_5 * forskjell_inntekt_5 ", kontroller_uten_belastning)),
                  data = masterdata, trace = FALSE)


#Modell M8
mod8k <- multinom(as.formula(paste("klynge_far ~ forskjell_kryssbelastning_5 * forskjell_inntekt_5 ", kontroller_uten_belastning)),
                  data = masterdata, trace = FALSE)


#Modell M1 (kjønnsasymmetri)
summary(mod1k)
cat("P-verdier modell 1k:\n"); print(pval(mod1k))

#Modell M2 (relative inntektsforskjeller)
summary(mod2k)
cat("P-verdier modell 2k:\n"); print(pval(mod2k))

#Modell M3 (relativ psykososial belastning)
summary(mod3k)
cat("P-verdier modell 3k:\n"); print(pval(mod3k))

#Modell M4 (relativ mekanisk belastning)
summary(mod4k)
cat("P-verdier modell 4k:\n"); print(pval(mod4k))

#Modell M5 (relativ kryssbelastning)
summary(mod5k)
cat("P-verdier modell 5k:\n"); print(pval(mod5k))

#Modell M6 (kryssbelastning far)
summary(mod6k)
cat("P-verdier modell 6k:\n"); print(pval(mod6k))

#########################

#Samle alle modellene i felles tabeller

library(stargazer)

# Tabell 1: modell M1 og modell M2
stargazer(
  mod1k, mod2k,
  type = "text",  
  title = "Modell 1 og 2",
  column.labels = c("M1 varig", "M1 moderat", "M2 varig", "M2 moderat"),
  digits = 3,
  star.cutoffs = c(0.05, 0.01, 0.001),
  keep.stat = c("n", "aic"),
  no.space = TRUE
)

# Tabell 2: modell M3, M4 og M5
stargazer(
  mod3k, mod4k, mod5k,
  type = "text",  
  title = "Modell 3, 4 og 5",
  column.labels = c("M3 varig", "M3 moderat","M4 varig", "M4 moderat", "M5 varig", "M5 moderat"),
  digits = 3,
  star.cutoffs = c(0.05, 0.01, 0.001),
  keep.stat = c("n", "aic"),
  no.space = TRUE
  )

# Tabell 3: modell M6
stargazer(
  mod6k,
  type = "text",  
  title = "Modell 6",
  column.labels = c("M6"),
  digits = 3,
  star.cutoffs = c(0.05, 0.01, 0.001),
  keep.stat = c("n", "aic"),
  no.space = TRUE
)

# Tabell 4: modell M7
stargazer(
  mod7k,
  type = "text",  
  title = "Modell 7",
  column.labels = c("M7 varig", "M7 moderat"),
  digits = 3,
  star.cutoffs = c(0.05, 0.01, 0.001),
  keep.stat = c("n", "aic"),
  no.space = TRUE
)

# Tabell 5: modell M8
stargazer(
  mod8k,
  type = "text",  
  title = "Modell 8",
  column.labels = c("M8"),
  digits = 3,
  star.cutoffs = c(0.05, 0.01, 0.001),
  keep.stat = c("n", "aic"),
  no.space = TRUE
)


##############################################################################
#MODELLER TIL GJENNOMSNITTLIGE MARGINALE EFFEKTER (AME)

library(marginaleffects)

#Modell M1
gme_mod1k <- avg_slopes(mod1k)
print(gme_mod1k)

#Modell M2
gme_mod2k <- avg_slopes(mod2k)
print(gme_mod2k)

# Modell M3
gme_mod3k <- avg_slopes(mod3k)
print(gme_mod3k)

# Modell M4
gme_mod4k <- avg_slopes(mod4k)
print(gme_mod4k)

# Modell M5
gme_mod5k <- avg_slopes(mod5k)
print(gme_mod5k)

# Modell M6
gme_mod6k <- avg_slopes(mod6k)
print(gme_mod6k)

# Modell M7
gme_mod7k <- avg_slopes(mod7k)
print(gme_mod7k)

# Modell M8
gme_mod8k <- avg_slopes(mod8k)
print(gme_mod8k)

# Tabell som kan limes inn i oppgaven

library(dplyr)
library(gt)

#################################################################################
# Lager tabell til del én av oppgaven

# Fikser tabell modell M3, M4 og M5
tabell_m345 <- bind_rows(
  gme_mod3k |> mutate(Modell = "M3 Psykososial"),
  gme_mod4k |> mutate(Modell = "M4 Mekanisk"),
  gme_mod5k |> mutate(Modell = "M5 Kryssbelastning")) |>
  filter(term %in% c("forskjell_strain_5", 
                     "forskjell_mekanisk_5", 
                     "forskjell_kryssbelastning_5")) |>
  select(Modell, group, contrast, estimate, std.error, p.value) |>
  mutate(
    group = case_when(
      group == "336" ~ "Varig",
      group == "527" ~ "Moderat",
      group == "73"  ~ "Referanse",
      TRUE ~ group ),
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ ""))

tabell_m345 |>
  gt() |>
  fmt_number(columns = c(estimate, std.error, p.value), decimals = 3) |>
  cols_label(
    Modell = "Modell",
    group = "Klynge",
    contrast = "Kontrast",
    estimate = "AME",
    std.error = "Std. feil",
    p.value = "p-verdi",
    stars = "") |>
  gtsave("tabell_m345.docx")

# Printer kontrollvariabler for å se på de
gme_mod3k |>
  filter(term %in% c("utd_nivaa_mor", "utd_nivaa_far",
                     "innvandrer_mor", "innvandrer_far",
                     "sykefravar_mor", "sykefravar_far",
                     "yrkinntekt_mor_2009_t", "yrkinntekt_far_2009_t"))

gme_mod4k |>
  filter(term %in% c("utd_nivaa_mor", "utd_nivaa_far",
                     "innvandrer_mor", "innvandrer_far",
                     "sykefravar_mor", "sykefravar_far",
                     "yrkinntekt_mor_2009_t", "yrkinntekt_far_2009_t"))

gme_mod5k |>
  filter(term %in% c("utd_nivaa_mor", "utd_nivaa_far",
                     "innvandrer_mor", "innvandrer_far",
                     "sykefravar_mor", "sykefravar_far",
                     "yrkinntekt_mor_2009_t", "yrkinntekt_far_2009_t"))

###################################################################################
#Lager tabell for del to av masteren

#Tabell modell M6
tabell_m6 <- gme_mod6k |>
  dplyr::filter(term == "forskjell_kryssbelastning_5") |>
  dplyr::select(term, group, contrast, estimate, std.error, p.value) |>
  dplyr::mutate(
    group = dplyr::case_when(
      group == "315" ~ "Moderat",
      TRUE ~ group),
    stars = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ ""))

tabell_m6 |>
  gt() |>
  fmt_number(columns = c(estimate, std.error, p.value), decimals = 3) |>
  cols_label(
    term = "Variabel",
    group = "Klynge",
    contrast = "Kontrast",
    estimate = "AME",
    std.error = "Std. feil",
    p.value = "p-verdi",
    stars = "") |>
  gtsave("tabell_m6.docx")


###############################################################
#APPENDIKSER

library(flextable)
library(broom)

# Gruppenavn
gruppenavn_mor <- c("336" = "Varig", "527" = "Moderat", "73" = "Referanse")
gruppenavn_far <- c("315" = "Moderat", "58" = "Referanse")

# Funksjon for AME-tabeller
lag_ame_tabell <- function(gme_objekter, modellnavn, tittel, filnavn, gruppenavn) {
  bind_rows(
    mapply(function(obj, navn) {
      as.data.frame(obj) |> mutate(Modell = navn)
    }, gme_objekter, modellnavn, SIMPLIFY = FALSE)) |>
    mutate(
      group = recode(group, !!!gruppenavn),
      p.value_fmt = case_when(
        p.value < 0.001 ~ paste0(round(estimate, 3), "***"),
        p.value < 0.01  ~ paste0(round(estimate, 3), "**"),
        p.value < 0.05  ~ paste0(round(estimate, 3), "*"),
        TRUE            ~ as.character(round(estimate, 3))
      )) |>
    select(Modell, group, term, contrast, p.value_fmt, std.error, p.value) |>
    flextable() |>
    colformat_double(j = c("std.error", "p.value"), digits = 3) |>
    set_header_labels(
      Modell      = "Modell",
      group       = "Klynge",
      term        = "Variabel",
      contrast    = "Kontrast",
      p.value_fmt = "AME",
      std.error   = "Std. feil",
      p.value     = "p-verdi" ) |>
    add_header_lines(values = tittel) |>
    add_footer_lines(values = "* p<0.05, ** p<0.01, *** p<0.001") |>
    autofit() |>
    save_as_docx(path = filnavn)}

#Tabell A1 (M1 og M2)
lag_ame_tabell(
  gme_objekter = list(gme_mod1k, gme_mod2k),
  modellnavn   = c("M1", "M2"),
  tittel       = "Tabell A1: M1 og M2 - Mors klynge",
  filnavn      = "appendiks_A1.docx",
  gruppenavn   = gruppenavn_mor)

#Tabell A2 (M3, M4, M5)
lag_ame_tabell(
  gme_objekter = list(gme_mod3k, gme_mod4k, gme_mod5k),
  modellnavn   = c("M3", "M4", "M5"),
  tittel       = "Tabell A2: M3, M4 og M5 - Mors klynge",
  filnavn      = "appendiks_A2.docx",
  gruppenavn   = gruppenavn_mor)

# Tabell A3 (M6) 
lag_ame_tabell(
  gme_objekter = list(gme_mod6k),
  modellnavn   = c("M6"),
  tittel       = "Tabell A3: M6 - Fars klynge",
  filnavn      = "appendiks_A3.docx",
  gruppenavn   = gruppenavn_far)

# Logg-odds-tabeller
lag_logo_tabell <- function(modell, tittel, filnavn, gruppenavn) {
  tidy(modell) |>
    mutate(
      group = recode(y.level, !!!gruppenavn),
      p.value_fmt = case_when(
        p.value < 0.001 ~ paste0(round(estimate, 3), "***"),
        p.value < 0.01  ~ paste0(round(estimate, 3), "**"),
        p.value < 0.05  ~ paste0(round(estimate, 3), "*"),
        TRUE            ~ as.character(round(estimate, 3))
      )) |>
    select(group, term, p.value_fmt, std.error, p.value) |>
    flextable() |>
    colformat_double(j = c("std.error", "p.value"), digits = 3) |>
    set_header_labels(
      group       = "Klynge",
      term        = "Variabel",
      p.value_fmt = "Logg-odds",
      std.error   = "Std. feil",
      p.value     = "p-verdi") |>
    add_header_lines(values = tittel) |>
    add_footer_lines(values = "* p<0.05, ** p<0.01, *** p<0.001") |>
    autofit() |>
    save_as_docx(path = filnavn)}

# Tabell A4 (M7) 
lag_logo_tabell(
  modell   = mod7k,
  tittel   = "Tabell A4: M7 - Trippelbyrde, mors klynge (logg-odds)",
  filnavn  = "appendiks_A4.docx",
  gruppenavn = gruppenavn_mor)

# Tabell A5 (M8)
lag_logo_tabell(
  modell   = mod8k,
  tittel   = "Tabell A5: M8 - Trippelbyrde, fars klynge (logg-odds)",
  filnavn  = "appendiks_A5.docx",
  gruppenavn = gruppenavn_far)
