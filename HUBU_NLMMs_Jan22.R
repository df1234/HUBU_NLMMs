##########################################################
# Project: HUBU - CORTICAL THINNING
# Date:    20/01/22
# Author:  Delia Fuhrmann, delia.fuhrmann@kcl.ac.uk
##########################################################

# Load libraries
packages <- c("tidyverse", "viridis", "qpcR",'reshape',
              'ggplot2','minpack.lm','saemix','lavaan',
              'psych','tidyverse','plyr','minpack.lm',
              'saemix','corrplot', 'ggseg', 'RVAideMemoire',
              'BayesFactor','patchwork','doBy')
invisible(lapply(packages, library, character.only = TRUE))
rm(packages)

# Read in data
HUBU  = read.csv("~/ownCloud/Rogier_Delia/Project HUBU Brains/Data/HUBU_freesurfer_data_for_Rogier_210620_euler20211014.csv")

# Set WD
setwd("~/ownCloud/Rogier_Delia/Project HUBU Brains/Scripts_JAN22")

# Filter out ppts with exclusion criterion
HUBU <- HUBU[HUBU$Include_freesurfer != 0,]
length(unique(HUBU$subject_id)) # 90

# Format data
cols.num <- c(6,8,10:296)
HUBU[cols.num] <- sapply(HUBU[cols.num],as.numeric)
cols.fac <- c(1:5,9)
HUBU[cols.fac] <- sapply(HUBU[cols.fac],as.character)
HUBU[cols.fac] <- sapply(HUBU[cols.fac],as.factor)
HUBU$sex[HUBU$sex == 0] = 'male'
HUBU$sex[HUBU$sex == 1] = 'female'
HUBU$sex = as.factor(HUBU$sex)

# Save data
saveRDS(HUBU, 'HUBU.rds')

##########################################################
# MCT schematic illustration

range03_5 <- function(x){(x-min(x))/(max(x)-min(x))}+2.2 #rescale to thickness normal values

x1 = -1*SSlogis(seq(from=0,to=30,length.out =1000), Asym = 100, xmid = 15, scal =2)

year<-seq(from=0,to=30,length.out =1000) 
id<-as.factor(rep(1,each=1000))
age<-year+rnorm(length(year),0,.1)
thickness<-range03_5(c(x1))

sigmoids_thinning<-data.frame(id,age,thickness)

ggplot(sigmoids_thinning,aes(age,thickness))+
  geom_smooth(col='grey', size = 0.9) +
  geom_point(aes(x=15, y=2.7), size = 2.1) +
  annotate(
    "text", label = "MCT",
    x = 17, y = 2.7, size = 4)+
  theme_classic()+
  theme(plot.title=element_text(size=18),
        axis.text=element_blank(),
        axis.ticks = element_blank()) +
  xlab("Age (late childhood to early adulthood)") +
  ylab("Cortical thickness")

ggsave(file='MCT.png',width=5,height=3)

##########################################################
# Summarize thickness data

HUBUThickness_temp = HUBU[,grep(c("*_thickness"), names(HUBU), value=T)]
HUBUThickness = HUBUThickness_temp[c(1:35)]

for(i in 1:35) {
  HUBUThickness[i] = (rowMeans(HUBUThickness_temp[,c(i,(i+35))], na.rm=T))
}

for ( col in 1:ncol(HUBUThickness)){
  colnames(HUBUThickness)[col] <-  sub("_thickness", "", colnames(HUBUThickness)[col])
  colnames(HUBUThickness)[col] <-  sub("rh_", "", colnames(HUBUThickness)[col])
}

HUBUThickness$subject_id = HUBU$subject_id
HUBUThickness$age_mri = HUBU$age_mri
HUBUThickness$sex = HUBU$sex
HUBUThickness$mr_round = as.factor(HUBU$mr_round)
levels(HUBUThickness$mr_round) = c(1:12)
HUBUThickness$mr_round = as.numeric(as.character(HUBUThickness$mr_round ))
HUBUThickness$euler = HUBU$mean_euler

# MRI data summary
sumdata1 =
  summaryBy(age_mri ~ subject_id, data = HUBUThickness,
            FUN = function(x) { c(min = min(x)) } )
sumdata1 = sumdata1[order(sumdata1$age_mri.min),] 
sumdata1$Participant = c(1:90)

sumdata2 = HUBU[,c(2,6)]
sumdata = merge(sumdata2, sumdata1, by = "subject_id")
sumdata$Age = sumdata$age_mri

s1 =
  ggplot(data = sumdata, aes(x = Age, y = Participant, group = Participant),
       show.legend = FALSE) +
  geom_line(aes(group=Participant, colour = Participant)) +
  geom_point(aes(group=Participant, colour = Participant)) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(7:22)) 

sumdata =
  summaryBy(mr_round + age_mri ~ subject_id, data = HUBUThickness,
            FUN = function(x) { c(m = mean(x), sum = sum(x),
                                  max = max(x), min = min(x)) } )

s2 =
  ggplot(sumdata, aes(x=mr_round.max)) + 
  geom_histogram(binwidth=1, fill="white", color="black") +
  scale_x_continuous(breaks = c(1:12)) +
  labs(x = "Number of MRI Scans", y ="Number of participants")

patchworks1 = s1/s2
patchworks1 + plot_annotation(tag_levels = 'A') + plot_layout(ncol=1,heights=c(3,1))
ggsave("HUBU_AgeMrRounds.png", width = 8, height = 10)

##########################################################
# Thinning across regions

HUBUThicknessDataPlot = melt(HUBUThickness[c(1:34,36:37)], id.vars = c("subject_id","age_mri"))

ggplot(data = HUBUThicknessDataPlot, aes(x=age_mri, y=value)) +
  geom_point(position= "jitter", aes(group = subject_id), size=1, alpha=.1) +
  geom_line (position= "jitter", aes(group = subject_id), alpha=.1) +
  stat_smooth(aes(group = 1, colour = variable), method="loess", se=F, size=2, formula = 'y~x') +
  facet_wrap(~ variable, scales = "free_y", ncol=4) +
  theme(legend.position="none", text = element_text(size=18)) +
  labs(x = "Age (years)", y = "Apparent cortical thickness (mm)")+
  scale_colour_viridis(discrete=T, option="D")

ggsave("HUBU_LoessPlot.png", width = 12, height = 15)

##########################################################
# Get starting values from NLM

M.4pl <- function(x, lower.asymp, upper.asymp, inflec, hill){
  f <- lower.asymp + ((upper.asymp - lower.asymp)/
                        (1 + (x / inflec)^-hill))
  return(f)
} # 4-parameter logistic function

x = HUBUThickness$age_mri
y = HUBUThickness$MeanThickness

nlslmfit = nlsLM(y ~ M.4pl(x, lower.asymp, upper.asymp, inflec, hill),
                 data = data.frame(x=x, y=y),
                 start = c(lower.asymp=min(y)+1E-10, upper.asymp=max(y)-1E-10, inflec=mean(x), hill=1),
                 control = nls.control(maxiter=1000, warnOnly=TRUE) )

summary(nlslmfit)

##########################################################
# Estimating inflection point and asymptotic for mean thickness

# Order data for saemix
HUBU_sub <- HUBUThickness[,c(36,37,35,38,40)]
HUBU_ordered <- HUBU_sub[order(HUBU_sub[,1], HUBU_sub[,2], HUBU_sub[,3], HUBU_sub[,4]),]

# Read data in
HUBU.data <- saemixData(name.data       = HUBU_ordered,
                        name.group      = "subject_id",
                        name.predictors = "age_mri",
                        name.response   = "MeanThickness",
                        name.covariates = c("sex","euler"))

# Define 4 parameter logistic function
logistic4.model <- function(psi, id, xidep){
  age_mri     <- xidep[, 1]   
  lower.asymp <- psi[id, 1]
  upper.asymp <- psi[id, 2]
  inflec      <- psi[id, 3]
  hill        <- psi[id, 4]
  resp <- lower.asymp + ((upper.asymp - lower.asymp)/(1 + (age_mri / inflec)^-hill))
  return(resp)
}

# Set up the saemix model
HUBU.model <- saemixModel(model = logistic4.model,
                          description = "4-param. logistic", 
                          psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                        ncol = 4, 
                                        byrow = TRUE, 
                                        dimnames = list(NULL, 
                                                        c("lower.asymp","upper.asymp","inflec","hill"))),
                          omega.init = diag(rep(0.5, 4)),
                          covariate.model = matrix(c(0, 1, 1, 0, 1, 0, 0, 0), ncol = 4, byrow = TRUE))

opt <- list(seed = 94352514, save = FALSE, save.graphs = FALSE)

HUBU.fit <- saemix(HUBU.model, HUBU.data, opt)

# pdf("obsvspred.pdf", width = 6.4, height = 3.8) 
# saemix.plot.select(HUBU.fit, observations.vs.predictions=T)
# dev.off() 
# 
# pdf("indfit.pdf", width = 6.4, height = 5.4) 
# saemix.plot.select(HUBU.fit, individual.fit = T)
# dev.off() 

# Get individual estimates
est <- psi(HUBU.fit)

cor.test(est$inflec, est$upper.asymp)
cor.test(est$inflec, est$lower.asymp)
cor.test(est$inflec, est$hill)

min(est$inflec)
max(est$inflec)

# Plot average thickness
est$subject_id = unique(HUBU_ordered$subject_id)
new = merge(est, HUBU_ordered, by = 'subject_id')
new$Thickness <- new$lower.asymp + ((new$upper.asymp - new$lower.asymp)/(1 + (new$age_mri / new$inflec)^-new$hill))

new$DeltaThickness = new$upper.asymp - new$lower.asymp

inflec_y <- 2.59 + ((2.95 - 2.59)/(1 + (14.72 / 14.72)^6.04))

f <- function(.x) 2.59 + ((2.95 - 2.59)/(1 + ((.x) / 14.72)^6.04))

p1 = ggplot(data = new, aes(x=age_mri, y=Thickness)) + 
  geom_point(position= "jitter", aes(group = subject_id), size=1, alpha=.1) + 
  geom_line (position= "jitter", aes(group = subject_id), alpha=.1) +
  stat_function(fun = f ) +
  theme(legend.position="none",text = element_text(size=16)) +
  labs(x = "Age (years)", y = "Predicted mean apparent thickness (mm)") +
  scale_colour_viridis(discrete=T,option="C") +
  annotate(geom = "point", x = 14.59, y = inflec_y, colour = "#FF7F50", size = 3) + 
  annotate(geom = "text",  x = 16.4, y = 2.75, label="MCT",
           color= "#FF7F50", size = 4)+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

p2 = ggplot(new, aes(x=inflec)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF7F50") +
  labs(x = "MCT (years)", y = "Density") +
  theme(legend.position="none",text = element_text(size=16))+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

my_cols <- c("#534666","#CD7672")

p3 = ggplot(data = new, aes(x=age_mri, y=Thickness)) +
  geom_point(position= "jitter", aes(group = subject_id, colour = sex), size=1, alpha=.4) +
  geom_line (position= "jitter", aes(group = subject_id, colour = sex), alpha=.4) +
  stat_smooth(aes(colour = sex), method="loess", se=T, size=1, formula = 'y~x') +
  theme(text = element_text(size=16)) +
  labs(x = "Age (years)", y = "Predicted mean apparent thickness (mm)",
       colour ="Gender")+
  scale_color_manual(values = my_cols) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

patchwork1 = p1 /(p2 + p3)
patchwork1 + plot_annotation(tag_levels = 'A')
ggsave("HUBU_MeanThickness.png", width = 8, height = 7)

##########################################################
# Run non-linear mixed effects models for each region

# results    <- list()
# inflection <- data.frame(start =factor(1:90))
# lowerAsymp <- data.frame(start =factor(1:90))
# upperAsymp <- data.frame(start =factor(1:90))
# hill       <- data.frame(start =factor(1:90))
# 
# # Define 4 parameter logistic function
# logistic4.model <- function(psi, id, xidep){
#   age_mri     <- xidep[, 1]
#   lower.asymp <- psi[id, 1]
#   upper.asymp <- psi[id, 2]
#   inflec      <- psi[id, 3]
#   hill        <- psi[id, 4]
#   resp <- lower.asymp + ((upper.asymp - lower.asymp)/(1 + (age_mri / inflec)^-hill))
#   return(resp)
# }
# 
# for (r in 1:35){
# # Order data for saemix
# HUBU_sub <- HUBUThickness[,c(36,37,r,38,40)]
# HUBU_ordered <- HUBU_sub[order(HUBU_sub[,1], HUBU_sub[,2], HUBU_sub[,3], HUBU_sub[,4]),]
# 
# # Read data in
# HUBU.data <- saemixData(name.data       = HUBU_ordered,
#                         name.group      = "subject_id",
#                         name.predictors = "age_mri",
#                         name.response   = colnames(HUBUThickness[r]),
#                         name.covariates = c("sex","euler"))
# 
# HUBU.model <- saemixModel(model = logistic4.model,
#                           description = "4-param. logistic",
#                           psi0 = matrix(c(2.6, 2.9, 15.3, -6.8),
#                                         ncol = 4,
#                                         byrow = TRUE,
#                                         dimnames = list(NULL,             c("lower.asymp","upper.asymp","inflec","hill"))),
#                           omega.init = diag(rep(0.5,4)),
#                           covariate.model = matrix(c(0,1,1,0,1,0,0,0), ncol = 4, byrow = TRUE))
# 
# opt <- list(seed = 94352514, save = FALSE, save.graphs = FALSE)
# 
# HUBU.fit <- saemix(HUBU.model, HUBU.data, opt)
# 
# # Get individual estimates
# est <- psi(HUBU.fit)
# 
# temp = list(data = HUBU.data, model = HUBU.fit, estimates = est)
# 
# results[[colnames(HUBUThickness[r])]] <- temp
# 
# inflection[,r] = est$inflec
# lowerAsymp[,r] = est$lower.asymp
# upperAsymp[,r] = est$upper.asymp
# hill[,r]       = est$hill
# 
# est  = NULL
# temp = NULL
# }
# 
# save(results,    file="NLMM_cov_results.RData")
# save(inflection, file="NLMM_cov_inflection.RData")
# save(lowerAsymp, file="NLMM_cov_lowerAsymp.RData")
# save(upperAsymp, file="NLMM_cov_upperAsymp.RData")
# save(hill,       file="NLMM_cov_hill.RData")

load("NLMM_cov_results.RData")
load("NLMM_cov_inflection.RData")
load("NLMM_cov_lowerAsymp.RData")
load("NLMM_cov_upperAsymp.RData")
load("NLMM_cov_hill.RData")

# Check CV (Coefficient of variation, should be < 20% https://www.scielo.br/scielo.php?script=sci_arttext&pid=S0103-84782013000600003&lng=en&tlng=en)

cv = c()
for (r in 1:35){
  x = summary(results[[r]][["model"]])
  y = t(data.frame(x[["fixed.effects"]][["CV(%)"]]))
  cv= rbind(cv,y)
}

write.csv(cv,"cv.csv")


estimates = c()
for (r in 1:35){
  x = summary(results[[r]][["model"]])
  y = t(data.frame(x[["fixed.effects"]][["Estimate"]]))
  estimates= rbind(estimates,y)
}

write.csv(estimates,"estimates.csv")

# EFA
fa.parallel(inflection[,c(1:4,6:11,13:14,16:19,21:31,34)])
fit_efa_i <- factanal(inflection[,c(1:4,6:11,13:14,16:19,21:31,34)], factors= 1, rotation="varimax")
print(fit_efa_i)

loadings = data.frame(unclass(fit_efa_i$loadings))
loadings$region = colnames(HUBUThickness[,c(1:4,6:11,13:14,16:19,21:31,34)])
write.csv(loadings,"loadings.csv")

# Plot regional effects
colnames(inflection) = colnames(HUBUThickness[,1:35])
inflection$subject_id = unique(HUBUThickness$subject_id)

InflectionDataPlot = melt(inflection[,c(1:4,6:11,13:14,16:19,21:31,34,36)], id.vars = c('subject_id'))

raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

sumld<- ddply(InflectionDataPlot, ~variable, summarise, mean = mean(value), median = median(value), lower = lb(value), upper = ub(value))

sumld <- sumld[order(sumld$mean),] 

inflecNew = merge(sumld, InflectionDataPlot, by='variable')

attach(inflecNew)
inflecNew = inflecNew[order(mean),]
detach(inflecNew)

inflecNew$split = c(rep(1,1260),rep(2,1260))
inflecNew1 = inflecNew[inflecNew$split == "1",]
inflecNew2 = inflecNew[inflecNew$split == "2",]

inflecNew1$split = factor(inflecNew1$split)
inflecNew2$split = factor(inflecNew2$split)

mean_fun <- function(x){
  return(data.frame(y = mean(x), label = round(mean(x), digits=2)))
}

source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

cloud1 =
  inflecNew1 %>%
  mutate(variable = fct_reorder(variable, value, .fun='mean')) %>%   
  ggplot(aes(y = value, x = reorder(variable, value), fill = variable)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = value, color = variable), position = position_jitter(width = .15),  size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  theme_bw() +
  raincloud_theme + 
  coord_flip()+
  scale_color_viridis(discrete=TRUE, option = 'inferno', begin = 0, end = 0.5) +
  scale_fill_viridis(discrete=TRUE, option = 'inferno', begin = 0, end = 0.5) +
  stat_summary(fun.data = mean_fun, geom = "text", position=position_nudge(x = -0.2, y = 0))+
  labs(y='MCT (years)', x='Brain region') 
cloud1

cloud2 =
  inflecNew2 %>%
  mutate(variable = fct_reorder(variable, value, .fun='mean')) %>%   
  ggplot(aes(y = value, x = reorder(variable, value), fill = variable)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = value, color = variable), position = position_jitter(width = .15),  size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  theme_bw() +
  raincloud_theme + 
  coord_flip()+
  scale_color_viridis(discrete=TRUE, option = 'inferno', begin = 0.5, end = 1) +
  scale_fill_viridis(discrete=TRUE, option = 'inferno', begin = 0.5, end = 1) +
  stat_summary(fun.data = mean_fun, geom = "text", position=position_nudge(x = -0.2, y = 0))+
  labs(y='MCT (years)', x='Brain region') 
cloud2

factor_regions = loadings$Factor1

region= c(
  "bankssts"                  , "caudal anterior cingulate",
  "caudal middle frontal"     , "cuneus"                   ,
  "fusiform"                 ,  
  "inferior parietal"         , "inferior temporal"        ,  
  "isthmus cingulate"         , "lateral occipital"        ,  
  "lateral orbitofrontal"     ,  
  "medial orbitofrontal"      , "middle temporal"           
  , "paracentral"              ,
  "pars opercularis"          , "pars orbitalis"           ,
  "pars triangularis"         , 
  "postcentral"               , "posterior cingulate"      ,  
  "precentral"                , "precuneus"                ,  
  "rostral anterior cingulate", "rostral middle frontal"   ,  
  "superior frontal"          , "superior parietal"        ,  
  "superior temporal"         , "supramarginal"            ,   
  "frontal pole"              ,
  "insula")  

brain_data_inflec = data.frame(region, factor_regions,
                               stringsAsFactors = FALSE)
rownames(brain_data_inflec) = NULL

brain =
  ggseg(.data=brain_data_inflec, mapping=aes(fill=factor_regions), 
        colour="grey", hemisphere = "left") + 
  labs(title="", fill="Loading") +
  scale_fill_viridis_c(option="mako", begin = 0.4) 
brain 

patchwork <- (cloud1 + cloud2) /  brain
patchwork +  plot_layout(height = c(6, 1))
ggsave("HUBU_InflectionRegions.pdf", width= 10, height = 15)

##########################################################
# CFA

model_OneFactor <-
  '
  g_inflection =~ 
                 bankssts
+ caudalanteriorcingulate
+     caudalmiddlefrontal
+                  cuneus
+                fusiform
+        inferiorparietal
+        inferiortemporal
+        isthmuscingulate
+        lateraloccipital
+    lateralorbitofrontal
+     medialorbitofrontal
+          middletemporal
+             paracentral
+         parsopercularis
+           parsorbitalis
+        parstriangularis
+             postcentral
+      posteriorcingulate
+              precentral
+               precuneus
+rostralanteriorcingulate
+    rostralmiddlefrontal
+         superiorfrontal
+        superiorparietal
+        superiortemporal
+           supramarginal
+             frontalpole
+                  insula                        
  '
fit_OneFactor_i <- cfa(model_OneFactor, data=inflection, missing="fiml")

model_OneFactor_intfixed_i <-
  '
  g_inflection =~ 
                 bankssts
+ caudalanteriorcingulate
+     caudalmiddlefrontal
+                  cuneus
+                fusiform
+        inferiorparietal
+        inferiortemporal
+        isthmuscingulate
+        lateraloccipital
+    lateralorbitofrontal
+     medialorbitofrontal
+          middletemporal
+             paracentral
+         parsopercularis
+           parsorbitalis
+        parstriangularis
+             postcentral
+      posteriorcingulate
+              precentral
+               precuneus
+rostralanteriorcingulate
+    rostralmiddlefrontal
+         superiorfrontal
+        superiorparietal
+        superiortemporal
+           supramarginal
+             frontalpole
+                  insula  

                bankssts ~ 14.72
 caudalanteriorcingulate ~ 14.72 
     caudalmiddlefrontal ~ 14.72 
                  cuneus ~ 14.72 
                fusiform ~ 14.72 
        inferiorparietal ~ 14.72 
        inferiortemporal ~ 14.72 
        isthmuscingulate ~ 14.72 
        lateraloccipital ~ 14.72 
    lateralorbitofrontal ~ 14.72 
     medialorbitofrontal ~ 14.72 
          middletemporal ~ 14.72 
             paracentral ~ 14.72 
         parsopercularis ~ 14.72 
           parsorbitalis ~ 14.72 
        parstriangularis ~ 14.72 
             postcentral ~ 14.72 
      posteriorcingulate ~ 14.72 
              precentral ~ 14.72 
               precuneus ~ 14.72 
rostralanteriorcingulate ~ 14.72 
    rostralmiddlefrontal ~ 14.72 
         superiorfrontal ~ 14.72 
        superiorparietal ~ 14.72 
        superiortemporal ~ 14.72 
           supramarginal ~ 14.72 
             frontalpole ~ 14.72 
                  insula ~ 14.72 
  '

fit_OneFactor_i <- cfa(model_OneFactor, data=inflection, missing="fiml")
fit_OneFactor_intfixed_i <- cfa(model_OneFactor_intfixed_i, data=inflection, missing="fiml")
anova(fit_OneFactor_intfixed_i, fit_OneFactor_i)

summary(fit_OneFactor_i, fit.measures=TRUE)

##########################################################
# Model selection

# Order data for saemix
HUBU_sub <- HUBUThickness[,c(36,37,35,38,40)]
HUBU_ordered <- HUBU_sub[order(HUBU_sub[,1], HUBU_sub[,2], HUBU_sub[,3], HUBU_sub[,4]),]

# Read data in
HUBU.data <- saemixData(name.data       = HUBU_ordered,
                        name.group      = "subject_id",
                        name.predictors = "age_mri",
                        name.response   = "MeanThickness",
                        name.covariates = c("sex","euler"))

# Define 4 parameter logistic function
logistic4.model <- function(psi, id, xidep){
  age_mri     <- xidep[, 1]   
  lower.asymp <- psi[id, 1]
  upper.asymp <- psi[id, 2]
  inflec      <- psi[id, 3]
  hill        <- psi[id, 4]
  resp <- lower.asymp + ((upper.asymp - lower.asymp)/(1 + (age_mri / inflec)^-hill))
  return(resp)
}

# Set up the saemix model
HUBU.model0000 <- saemixModel(model = logistic4.model,
                          description = "4-param. logistic", 
                          psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                        ncol = 4, 
                                        byrow = TRUE, 
                                        dimnames = list(NULL, 
                                                        c("lower.asymp","upper.asymp","inflec","hill"))),
                          omega.init = diag(rep(0.5, 4)),
                          covariate.model = matrix(c(0, 0, 0, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE))

opt <- list(seed = 94352514, save = FALSE, save.graphs = FALSE)

HUBU.fit0000 <- saemix(HUBU.model0000, HUBU.data, opt)



HUBU.model10000000 <- saemixModel(model = logistic4.model,
                              description = "4-param. logistic", 
                              psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                            ncol = 4, 
                                            byrow = TRUE, 
                                            dimnames = list(NULL, 
                                                            c("lower.asymp","upper.asymp","inflec","hill"))),
                              omega.init = diag(rep(0.5, 4)),
                              covariate.model = matrix(c(1, 0, 0, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE))

HUBU.fit1000000 <- saemix(HUBU.model10000000, HUBU.data, opt)


HUBU.model01000000 <- saemixModel(model = logistic4.model,
                                  description = "4-param. logistic", 
                                  psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                                ncol = 4, 
                                                byrow = TRUE, 
                                                dimnames = list(NULL, 
                                                                c("lower.asymp","upper.asymp","inflec","hill"))),
                                  omega.init = diag(rep(0.5, 4)),
                                  covariate.model = matrix(c(0, 1, 0, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE))

HUBU.fit0100000 <- saemix(HUBU.model01000000, HUBU.data, opt)


HUBU.model01100000 <- saemixModel(model = logistic4.model,
                                  description = "4-param. logistic", 
                                  psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                                ncol = 4, 
                                                byrow = TRUE, 
                                                dimnames = list(NULL, 
                                                                c("lower.asymp","upper.asymp","inflec","hill"))),
                                  omega.init = diag(rep(0.5, 4)),
                                  covariate.model = matrix(c(0, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE))

HUBU.fit0110000 <- saemix(HUBU.model01100000, HUBU.data, opt)


HUBU.model01110000 <- saemixModel(model = logistic4.model,
                                  description = "4-param. logistic", 
                                  psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                                ncol = 4, 
                                                byrow = TRUE, 
                                                dimnames = list(NULL, 
                                                                c("lower.asymp","upper.asymp","inflec","hill"))),
                                  omega.init = diag(rep(0.5, 4)),
                                  covariate.model = matrix(c(0, 1, 1, 1, 0, 0, 0, 0), ncol = 4, byrow = TRUE))

HUBU.fit0111000 <- saemix(HUBU.model01110000, HUBU.data, opt)


HUBU.model01101000 <- saemixModel(model = logistic4.model,
                                  description = "4-param. logistic", 
                                  psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                                ncol = 4, 
                                                byrow = TRUE, 
                                                dimnames = list(NULL, 
                                                                c("lower.asymp","upper.asymp","inflec","hill"))),
                                  omega.init = diag(rep(0.5, 4)),
                                  covariate.model = matrix(c(0, 1, 1, 0, 1, 0, 0, 0), ncol = 4, byrow = TRUE))

HUBU.fit01101000 <- saemix(HUBU.model01101000, HUBU.data, opt)


HUBU.model01101100 <- saemixModel(model = logistic4.model,
                                  description = "4-param. logistic", 
                                  psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                                ncol = 4, 
                                                byrow = TRUE, 
                                                dimnames = list(NULL, 
                                                                c("lower.asymp","upper.asymp","inflec","hill"))),
                                  omega.init = diag(rep(0.5, 4)),
                                  covariate.model = matrix(c(0, 1, 1, 0, 1, 1, 0, 0), ncol = 4, byrow = TRUE))

HUBU.fit01101100 <- saemix(HUBU.model01101100, HUBU.data, opt)


HUBU.model01101010 <- saemixModel(model = logistic4.model,
                                  description = "4-param. logistic", 
                                  psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                                ncol = 4, 
                                                byrow = TRUE, 
                                                dimnames = list(NULL, 
                                                                c("lower.asymp","upper.asymp","inflec","hill"))),
                                  omega.init = diag(rep(0.5, 4)),
                                  covariate.model = matrix(c(0, 1, 1, 0, 1, 0, 1, 0), ncol = 4, byrow = TRUE))

HUBU.fit01101010 <- saemix(HUBU.model01101010, HUBU.data, opt)


HUBU.model01101001 <- saemixModel(model = logistic4.model,
                                  description = "4-param. logistic", 
                                  psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                                ncol = 4, 
                                                byrow = TRUE, 
                                                dimnames = list(NULL, 
                                                                c("lower.asymp","upper.asymp","inflec","hill"))),
                                  omega.init = diag(rep(0.5, 4)),
                                  covariate.model = matrix(c(0, 1, 1, 0, 1, 0, 0, 1), ncol = 4, byrow = TRUE))

HUBU.fit01101001 <- saemix(HUBU.model01101001, HUBU.data, opt)


# Final model
HUBU.model01101000 <- saemixModel(model = logistic4.model,
                                  description = "4-param. logistic", 
                                  psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                                ncol = 4, 
                                                byrow = TRUE, 
                                                dimnames = list(NULL, 
                                                                c("lower.asymp","upper.asymp","inflec","hill"))),
                                  omega.init = diag(rep(0.5, 4)),
                                  covariate.model = matrix(c(0, 1, 1, 0, 1, 0, 0, 0), ncol = 4, byrow = TRUE))

opt <- list(seed = 94352514, save = FALSE, save.graphs = FALSE)

HUBU.fit01101000 <- saemix(HUBU.model01101000, HUBU.data, opt)
