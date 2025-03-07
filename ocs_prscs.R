BHRC_Controles_OCD -> BHRC_Controles_OCD_2var
BHRC_Controles_OCD_2var$ident <- NULL
BHRC_Controles_OCD_2var$subjectid <- NULL
BHRC_Controles_OCD_2var$IID <- NULL
BHRC_Controles_OCD_2var$AnyDis <- NULL
BHRC_Controles_OCD_2var$OCD<- NULL
fid <- rep(0, 1425)
BHRC_Controles_OCD_2var$FID <- fid
write.table(BHRC_Controles_OCD_2var, file = 'C:/Users/diego/OneDrive/Área de Trabalho/Unifesp/PRS-CS//BHRC_Controles_OCD.txt', quote = F, row.names = F, col.names = F)

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

mean <- ddply(prscs_alvos_OCD, "OCD", summarise, grp.mean=mean(CS))
table(prscs_alvos_OCD$OCD)


ggplot(prscs_alvos_OCD, aes(x=CS, group=OCD, fill=OCD)) +
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







R2_prscs_ocd <- (glm(OCD ~ CS + Sex, data=prscs_alvos_OCD, family=binomial("logit")))
summary(R2_prscs_ocd)
NagelkerkeR2(R2_prscs_ocd)