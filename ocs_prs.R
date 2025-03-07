BHRC_Controles_OCD -> BHRC_Controles_OCD_2var
BHRC_Controles_OCD_2var$ident <- NULL
BHRC_Controles_OCD_2var$subjectid <- NULL
BHRC_Controles_OCD_2var$IID <- NULL
BHRC_Controles_OCD_2var$AnyDis <- NULL
BHRC_Controles_OCD_2var$OCD<- NULL
fid <- rep(0, 1425)
BHRC_Controles_OCD_2var$FID <- fid
write.table(BHRC_Controles_OCD_2var, file = 'PRS-CS//BHRC_Controles_OCD.txt', quote = F, row.names = F, col.names = F)

ocd_prscs <- read.table("OCD_BHRC_Probands_prscs.profile", header = T)
BHRC_eigenvec_ocd <- read.table("arquivo_pca_OCD.eigenvec", h=F)
colnames(BHRC_eigenvec_ocd) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
library(plyr)
prscs_alvos_OCD <- join_all(list(BHRC_Controles_OCD, ocd_prscs, BHRC_Probands_TABELONA, BHRC_eigenvec_ocd),by = "IID", type="inner")
prscs_alvos_OCD$ident <- NULL
prscs_alvos_OCD$subjectid <- NULL
prscs_alvos_OCD$AnyDis <- NULL
prscs_alvos_OCD$CNT2 <- NULL
prscs_alvos_OCD$PHENO <- NULL

# Grafico 
library(dplyr)
library(plyr)
library(ggplot2)
library(rcompanion)

fit <- lm (SCORESUM ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=prscs_alvos_OCD)
prscs_alvos_OCD$CS <- residuals(fit)

mean <- ddply(tabelonaOCD_scores_allpgs, "OCD", summarise, grp.mean=mean(CS))
table(tabelonaOCD_scores_allpgs$OCD)


ggplot(tabelonaOCD_scores_allpgs, aes(x=CS, group=OCD, fill=OCD)) +
  geom_density(alpha=0.5) +
  geom_vline(data=mean, aes(xintercept=grp.mean, colour=OCD), linetype="dashed", linewidth=0.8) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_classic() +
  labs(x="Polygenic Score", y="Density")

Null <- glm(formula = OCD ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex, family = binomial(link="logit"), data = prscs_alvos_OCD)
Full <- glm(formula = OCD ~ SCORESUM + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex, family = binomial(link="logit"), data = prscs_alvos_OCD)
nagelkerke(Full, null=Null)$Pseudo.R.squared.for.model.vs.null[3]
summary (Full)


## OCD PRSice
library(plyr)
BHRC_sex_pcs_OCD <- join_all(list(BHRC_Sex, BHRC_eigenvec_ocd),by = "IID", type="inner")
View(BHRC_sex_pcs_OCD)
fid <- rep(0, 1213)
BHRC_sex_pcs_OCD$FID <- fid
BHRC_sex_pcs_OCD <- BHRC_sex_pcs_OCD [c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Sex")]
# BHRC_Controles_OCD_2var usada pra tabela com os fenotipos e a que eu fiz em cima pras covariaveis
write.table (BHRC_sex_pcs_OCD, file='C:/Users/diego/OneDrive/Área de Trabalho/Unifesp/OCD///BHRC_Sex_pcs_OCD.txt', quote = F, row.names = F)
write.table (BHRC_Controles_OCD_2var, file='C:/Users/diego/OneDrive/Área de Trabalho/Unifesp/OCD///BHRC_fenotipo_OCD.txt', quote = F, row.names = F)
# Pos PRS
setwd("C:/Users/diego/OneDrive/Área de Trabalho/Unifesp/OCD")
prsice_scores_ocd <- read.table("covsexpcs_bhrc_OCD.best", h=T)
prsice_all_scores_ocd <- read.table("covsexpcs_bhrc_OCD.all_score", h=T)
prsice_scores_ocd$FID <- NULL
scores_ocd <- join_all(list(prscs_alvos_OCD, prsice_scores_ocd), by= "IID", type="inner")

# ggrafico


fit <- lm (PRS ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=scores_ocd)
scores_ocd$Sice <- residuals(fit)

mean <- ddply(tabelonaOCD_scores_allpgs, "OCD", summarise, grp.mean=mean(Sice))
table(scores_ocd$OCD)

library(ggplot2)
ggplot(tabelonaOCD_scores_allpgs, aes(x=Sice, group=OCD, fill=OCD)) +
 geom_histogram()+
 #geom_density(alpha=0.5) +
  geom_vline(data=mean, aes(xintercept=grp.mean, colour=OCD), linetype="dashed", linewidth=0.8) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_classic() +
  labs(x="Polygenic Score", y="Density")


tabelonaOCD_scores_allpgs$OCD <- ifelse (tabelonaOCD_scores_allpgs$OCD == 1, "Case", "Control")
tabelonaOCD_scores_allpgs$OCD <- as.factor(tabelonaOCD_scores_allpgs$OCD)

Null <- glm(formula = OCD ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex, family = binomial(link="logit"), data = tabelonaOCD_scores_allpgs)
Full <- glm(formula = OCD ~ ZScore_PRSice + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex, family = binomial(link="logit"), data = tabelonaOCD_scores_allpgs)
nagelkerke(Full, null=Null)$Pseudo.R.squared.for.model.vs.null[3]
summary (Full)


############## BOX PLOTS OCD ###################
scores_ocd$Zcs <- scale(scores_ocd$CS)
scores_ocd$Zsice <- scale(scores_ocd$Sice)
scores_ocd -> tabelonaOCD_scores_allpgs
scores_ocd <- tabelonaOCD_scores_allpgs[, c("IID", "Zcs", "Zsice", "OCD")]
scores_ocd -> scores_cs_ocd 
scores_ocd -> scores_sice_ocd 
scores_cs_ocd$Zsice <- NULL
scores_sice_ocd$Zcs <- NULL
PRSCS <- rep("PRS-CS", 1213) 
scores_cs_ocd$Ferramenta <- PRSCS
PRSice <- rep("PRSice", 1213) 
scores_sice_ocd$Ferramenta <- PRSice

scores_cs_ocd$Zsice <- rep(NA, 1213)
scores_cs_ocd$Zsice <- scale(scores_cs_ocd$Zsice)
scores_sice_ocd$Zcs <- rep(NA, 1213)
scores_sice_ocd$Zcs <- scale(scores_sice_ocd$Zcs)
scores_ocd <- rbind(scores_cs_ocd, scores_sice_ocd)
# Uma variavel chamada Score
scores_ocd[is.na(scores_ocd)] <- 0
scores_ocd$Score <- scores_ocd$Zcs + scores_ocd$Zsice

scores_ocd$OCD <- ifelse (scores_ocd$OCD == 1, "Case", "Control")

(scores_ocd %>% # Data.frame
    ggplot(aes(y = Score, x = OCD, color = OCD)) + # Coluna de interesse 
    geom_boxplot()+ # Formato de boxplot
    geom_jitter(alpha = 0.4, position = position_jitter(width = 0.2, seed = 0)) +
    theme_bw()+ # Tema com fundo branco e contorno preto
    theme(plot.title = element_text(size=20, face="bold.italic"), # Aumentar o tamanho do título e deixar em negrito/itálico
          plot.subtitle = element_text(size=12, color="gray34", face="bold.italic"), # Aumentar o tamanho do subtítulo e deixar em negrito/itálico
          strip.text = element_text(size=12, color="white"), # Alterar a cor do texto da faceta
          strip.background = element_rect (fill="gray20"))+ # Alterar a cor de fundo da faceta 
    labs(title = "Polygenic Score", # Título
         subtitle = "Comparing PGS for OCD", # Subtítulo
         x = "", # Título do eixo x
         y = "Z-Score", # Título do eixo y
         color = "OCD") + # Título da legenda
    scale_color_brewer(palette="Set1")+ # Cor dos gráficos
    facet_wrap(~Ferramenta)) # Separar em facetas

## Um pt mais rigoroso

prsice_all_scores_ocd <- join_all(list(prscs_alvos_OCD, prsice_all_scores_ocd), by= "IID", type="inner") # Outro pt

fit <- lm (Pt_0.001 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=prsice_all_scores_ocd)
prsice_all_scores_ocd$Sice_Pt0.001 <- residuals(fit)

mean <- ddply(prsice_all_scores_ocd, "OCD", summarise, grp.mean=mean(Sice_Pt0.001))
table(scores_ocd$OCD)

library(ggplot2)
ggplot(prsice_all_scores_ocd_SemFID, aes(x=Sice_Pt0.001, group=OCD, fill=OCD)) +
 geom_histogram(alpha=0.5) +
  # geom_density(alpha=0.5) +
  geom_vline(data=mean, aes(xintercept=grp.mean, colour=OCD), linetype="dashed", linewidth=0.8) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_classic() +
  labs(x="Polygenic Score", y="Density")


tabelonaOCD_scores_allpgs$OCD <- ifelse (tabelonaOCD_scores_allpgs$OCD == 1, "Case", "Control")
tabelonaOCD_scores_allpgs$OCD <- as.factor(tabelonaOCD_scores_allpgs$OCD)

Null <- glm(formula = OCD ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex, family = binomial(link="logit"), data = tabelonaOCD_scores_allpgs)
Full <- glm(formula = OCD ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Sex, family = binomial(link="logit"), data = tabelonaOCD_scores_allpgs)
nagelkerke(Full, null=Null)$Pseudo.R.squared.for.model.vs.null[3]
summary (Full)


#### Pt mais restrito 
mean <- ddply(prsice_all_scores_ocd, "OCD", summarise, grp.mean=mean(Sice_Pt1e05))
table(prsice_all_scores_ocd$OCD)

library(ggplot2)
ggplot(prsice_all_scores_ocd_SemFID, aes(x=Sice_Pt1e05, group=OCD, fill=OCD)) +
  geom_histogram(alpha=0.5) +
  # geom_density(alpha=0.5) +
  geom_vline(data=mean, aes(xintercept=grp.mean, colour=OCD), linetype="dashed", linewidth=0.8) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_classic() +
  labs(x="Polygenic Score", y="Density")
