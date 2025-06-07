library(tidyverse)
library(sjmisc)
library(hablar)
library(TAM)
library(psych)
library(lordif)

library(gridExtra)
library(cowplot)
library(ggpubr)
library(extrafont)
library(paletteer)
library(MOTE)
library(flextable)

#
# grab data & responses----
samb<-read_csv("SAmb_R.csv", col_names = TRUE)
resp<-samb[1:19]

psych::describe(samb)
psych::describe(resp)

table(samb$gender) # 0=male, 1=female
table(samb$ethnic) # 1=EuroAm, 2=AfAm, 3=Hispn/Latinx, 4=AsianAm, 5=NativeAm, 6=Other/Decline


#
# reliability & dimensionality----
alpha_resp<-psych::alpha(resp)
alpha_resp$total

omega(resp, fm="ml", poly=TRUE)

# item-scale correlations, from low to high
sort(alpha_resp$item.stats$raw.r) %>% round(2)


# parallel analysis
fpar<-fa.parallel(resp, fm="ml", fa="fa", cor="tet", sim=FALSE,
            error.bars=FALSE, se.bars=FALSE, n.iter=50, correct=0)
fpar$fa.values

# MAP
vss(resp, fm="ml", cor="tet", n=4, correct=0)

# essential unidimensionality: ratio of first and second eigenvalues
fpar$fa.values
(fpar$fa.values[1]/fpar$fa.values[2]) %>% round(3)

# EFA with bifactor rotation
bi_efa<-fa(resp, fm="ml", cor="tet", nfactors = 4, rotate = "bifactor", correct=0)
bi_efa

bi_efa$Structure


#
# IRT Models in TAM----

# Rasch model
tr<-TAM::tam.mml(resp, control = list(Msteps=10))

# 2PL IRT model
tp<-TAM::tam.mml.2pl(resp, irtmodel = "2PL", control = list(Msteps=10))

summary(tr)
summary(tp)


# information criteria fit measures
tr$ic
tp$ic
anova(tr, tp) # 2PL fits better based on AIC, BIC, GHP

# 2PL reliability
tp$EAP.rel


# local dependence
tres<-TAM::tam.modelfit(tp)
tres$Q3_summary
tres$stat.itempair[1:10,c(1,2,6)] # specific item pairs, sorted by absolute size


# item fit
i.fit<-msq.itemfit(tp) # infit, outfit
r.fit<-IRT.itemfit(tp) # RMSD

summary(i.fit)
summary(r.fit)


#
# DIF via lordif----
table(samb$ethnic)
table(samb$gender)

# recodes ethnicity (Decline/Other category is set to missing)
sambdif<-samb %>%
  select(1:19, SAmb_tot, gender, ethnic) %>%
  mutate(ethnic_fac = rec(ethnic, rec = "1=0; 6=NA; else=1")) %>%
  convert(fct(gender, ethnic_fac))

table(sambdif$gender) # 0=women, 1=men
table(sambdif$ethnic_fac) #0=EuroAm, 1=PoC

# lordif is picky about "data frame" vs matrix, tibble, etc...
sgender<-sambdif$gender
sethnic<-sambdif$ethnic_fac
respdif<-data.frame(sambdif[1:19])

# gender DIF
gendif.r2<-lordif(respdif, sgender, model = "GPCM",
               criterion = "R2", 
               pseudo.R2 = "McFadden", # this is the criterion
               R2.change =.01, # this is the % R2 to use as a flagging threshold
               maxIter = 1000, minCell = 25)
summary(gendif.r2)

ethdif.r2<-lordif(respdif, sethnic, model = "GPCM",
                  criterion = "R2", 
                  pseudo.R2 = "McFadden", # this is the criterion
                  R2.change =.01, # this is the % R2 to use as a flagging threshold
                  maxIter = 1000, minCell = 25)
summary(ethdif.r2)

# viewing the 2 items flagged at 1% DIF
ethdif.r2$stats[c(13,16),c(1,7:9)]

# some plots for these 2 items
plot.lordif(ethdif.r2, labels = c("White","PoC"))


#
# plots----

# you might not have these palettes and fonts on your computer

# colors
paletteer_d("yarrr::info")

# Oswald is freely available on Google Fonts
# view fonts on your computer
extrafont::fonts()

# scree plot
scree<-data.frame("Actual" = fpar$fa.values,
                  "Resampled" = fpar$fa.simr,
                  "Item" = seq.int(1, 19, 1))

splot<-pivot_longer(scree, cols=1:2, names_to="Method") %>%
  filter(Item<8)

eig<-ggline(splot, x="Item", y="value", group = "Method",
               color="Method", size=1.1,
               palette = paletteer_d("yarrr::info", direction=1),
               xlab="Factor", ylab="Eigenvalue")+
  theme_minimal_hgrid(font_size=12, font_family = "Oswald")
eig<-ggpar(eig, ylim = c(-.5, 7.5), yticks.by=2, 
      legend = "top", legend.title = "")
eig

# item fit plots
fits<-as.data.frame(i.fit$itemfit[,c(1,3,6)])
fitplot<-full_join(fits, r.fit$RMSD)
fitplot$item<-seq.int(1, 19, 1)
colnames(fitplot)<-c("Item","Outfit","Infit","RMSD")

# infit outift
inout<-gather(fitplot, "Outfit","Infit",
              key="Metric", value="Fit")

item_inout<-ggbarplot(inout, x="Item", y="Fit",
                       fill = "Metric", group.by = "Metric",
                       position = position_dodge(0.77),
                       color="white", 
                       palette = paletteer_d("yarrr::info", direction=1),
                       xlab = "SAS Item")+
  theme_minimal_hgrid(font_size=12, font_family = "Oswald")
item_inout<-ggpar(item_inout, ylim = c(.8, 1.1), yticks.by=.1, 
      legend = "top", legend.title = "")
item_inout

item_rmsd<-ggbarplot(fitplot, x="Item", y="RMSD",
                      color="white", 
                      fill = "#6B8993FF",
                      xlab = "SAS Item", ylab = "RMSD")+
  theme_minimal_hgrid(font_size=12, font_family = "Oswald")
item_rmsd<-ggpar(item_rmsd, ylim = c(.0, .02), yticks.by=.005, 
      legend = "none", legend.title = "")
item_rmsd


# difficulty values from TAM
tamdat<-tp$item_irt
tamdat$item<-seq.int(1, 19, 1)

diffdot<-ggdotchart(tamdat, x="item", y="beta", size=3,
                  color = "#E7695DFF",
                  dot.size = 3.5, rotate = FALSE,
                  xlab = "SAS Item", ylab = "Difficulty")+
    theme_minimal_grid(font_family = "Oswald")
diffdot<-ggpar(diffdot, x.text.angle = 0, legend = "top", legend.title = "")
diffdot

# discrimination values from TAM
discdot<-ggdotchart(tamdat, x="item", y="alpha", size=3,
                      color = "#6B8993FF",
                      dot.size = 3.5, rotate = FALSE,
                      xlab = "SAS Item", ylab = "Discrimination")+
    theme_minimal_grid(font_family = "Oswald")
discdot<-ggpar(discdot, x.text.angle = 0, legend = "top", legend.title = "")
discdot


# test information function
tinfo<-IRT.informationCurves(tp,theta=seq(-5,5, len=100))
ttheta<-round(seq(-5,5, len=100),2)
test.info<-data.frame(tinfo$theta, tinfo$test_info_curve)
colnames(test.info)<-c("Theta","scamb")

tinfo_plot<-ggplot(test.info)+
  geom_line(size=1.25,mapping=aes(x=Theta, y = scamb), 
            color="#6B8993FF")+
  labs(y="Test Information")+
  scale_x_continuous(name = "Theta", breaks=c(-5,-4,-3,-2,-1,0,1,2,3,4,5))+
  theme_minimal_hgrid(font_size=12, font_family = "Oswald")+
  theme(legend.position = "top", legend.title =element_blank())
tinfo_plot


#
# table of item features----
item_words<-c("Often I feel like I hate even my favorite activities.",
              "My thoughts and feelings always seem to be contradictory.",
              "My feelings about my own worth as a person are constantly changing back and forth.",
              "Very often when I feel like doing something, at the same time I don’t feel like doing it.",
              "When I am trying to make a decision, it almost feels like I am physically switching from side to side.",
              "It’s impossible to know how you feel because the people around you are constantly changing.",
              "I always seem to be the most unsure of myself at the same time that I am most confident of myself.",
              "I always seem to have difficulty deciding what I would like to do.",
              "Most people seem to know what they’re feeling more easily than I do.",
              "Love and hate tend to go together.",
              "Love never seems to last very long.",
              "The closer I get to people, the more I am annoyed by their faults.",
              "Everyone has a lot of hidden resentment toward his loved one.",
              "I have noticed that feelings of tenderness often turn into feelings of anger.",
              "My experiences with love have always been muddled with great frustrations.",
              "I usually find that feelings of hate will interfere when I have grown to love someone.",
              "A sense of shame has often interfered with my accepting words of praise from others.",
              "I usually experience doubt when I have accomplished something that I have worked on for a long time.",
              "I doubt if I can ever be sure exactly what my true interests are.")

descr_resp<-psych::describe(resp)

table_data<-tibble("Item" = tamdat$item,
                   "Text" = item_words,
                   "Mean (Percent Endorsed)" = apa(descr_resp$mean, 2, F),
                   "SD" = apa(descr_resp$sd, 2, F),
                   "Item-Scale Correlation" = apa(alpha_resp$item.stats$raw.r, 2, F),
                   "IRT Difficulty" = apa(tamdat$beta, 2, T),
                   "IRT Discrimination" = apa(tamdat$alpha, 2, T),
                   "Gender DIF (R2)" = apa(gendif.r2$stats$pseudo13.McFadden,3, F),
                   "Race-Ethnicity DIF (R2)" = apa(ethdif.r2$stats$pseudo13.McFadden,3, F))

#
# export plots and tables----

# table 1
flextable(table_data) %>% 
  save_as_docx(path = "Samb_Table1.docx")

# scree plot
tiff("Fig1_Scree.tif", res=300, width=4, height=3, units="in")
eig
dev.off()

# item features
tiff("Fig2_DiffDisc.tif", res=300, width=5.5, height=6, units="in")
grid.arrange(diffdot, discdot, nrow=2, ncol=1)
dev.off()

# test information
tiff("Fig3_TestInfo.tif", res=300, width=4, height=3, units="in")
tinfo_plot
dev.off()

# item fit
tiff("Fig4_InfitOutfit.tif", res=300, width=5.5, height=3, units="in")
item_inout
dev.off()

tiff("Fig5_RMSDfit.tif", res=300, width=5.5, height=3, units="in")
item_rmsd
dev.off()

# END #