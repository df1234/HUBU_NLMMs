##########################################################
# Project: HUBU - CORTICAL THINNING
# Date:    05/05/22
# Author:  Delia Fuhrmann, delia.fuhrmann@kcl.ac.uk
##########################################################

# Set WD
setwd("~/ownCloud/Rogier_Delia/Project HUBU Brains/Scripts_May22")

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

ggsave(file='MCT.tiff',width=5,height=3)

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

# Euler summary
sumdataE =
  summaryBy(euler ~ subject_id, data = HUBUThickness)

mean(sumdataE$euler.mean)
require(plotrix)
std.error(sumdataE$euler.mean)
max(sumdataE$euler.mean)
min(sumdataE$euler.mean)

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
patchworks1 + plot_annotation(tag_levels = 'A') + plot_layout(ncol=1,heights=c(3,1)) & 
  theme(plot.tag = element_text(size = 13, face="bold"))
ggsave("HUBU_AgeMrRounds.tiff", width = 8, height = 10)

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

ggsave("HUBU_LoessPlot.tiff", width = 12, height = 15)

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
                          covariate.model = matrix(c(1, 1, 1, 0, 1, 0, 0, 0), ncol = 4, byrow = TRUE))

opt <- list(seed = 94352514, save = FALSE, save.graphs = FALSE)

HUBU.fit <- saemix(HUBU.model, HUBU.data, opt)

pdf("obsvspred.pdf", width = 6.4, height = 3.8) 
saemix.plot.select(HUBU.fit, observations.vs.predictions=T)
dev.off() 

pdf("indfit.pdf", width = 6.4, height = 5.4) 
saemix.plot.select(HUBU.fit, individual.fit = T)
dev.off() 

# Get individual estimates
est <- psi(HUBU.fit)

cor.test(est$inflec, est$upper.asymp)
cor.test(est$inflec, est$lower.asymp)
cor.test(est$inflec, est$hill)

min(est$inflec)
max(est$inflec)

##########################################################
# Plot average thickness
est$subject_id = unique(HUBU_ordered$subject_id)
new = merge(est, HUBU_ordered, by = 'subject_id')
new$Thickness <- new$lower.asymp + ((new$upper.asymp - new$lower.asymp)/(1 + (new$age_mri / new$inflec)^-new$hill))

new$DeltaThickness = new$upper.asymp - new$lower.asymp

inflec_y <- 2.62 + ((2.95 - 2.62)/(1 + (14.36 / 14.36)^6.34))

f <- function(.x) 2.62 + ((2.95 - 2.62)/(1 + ((.x) / 14.36)^6.34))

p1 = ggplot(data = new, aes(x=age_mri, y=Thickness)) + 
  geom_point(position= "jitter", aes(group = subject_id), size=1, alpha=.1) + 
  geom_line (position= "jitter", aes(group = subject_id), alpha=.1) +
  stat_function(fun = f ) +
  theme(legend.position="none",text = element_text(size=16)) +
  labs(x = "Age (years)", y = "Predicted mean apparent thickness (mm)") +
  scale_colour_viridis(discrete=T,option="C") +
  annotate(geom = "point", x = 14.36, y = inflec_y, colour = "#FF7F50", size = 3) + 
  annotate(geom = "text",  x = 15.4, y = 2.785, label="MCT",
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
       colour ="Sex")+
  scale_color_manual(values = my_cols) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

patchwork1 = p1 /(p2 + p3)
patchwork1 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 13, face="bold"))
ggsave("HUBU_MeanThickness.tiff", width = 8, height = 7)

##########################################################
# Comparison to Linear Model

lin.model <- function(psi, id, xidep){
  age_mri     <- xidep[, 1]   
  int         <- psi[id, 1]
  slope       <- psi[id, 2]
  resp <- int + slope*age_mri
  return(resp)
}

HUBU.modellin <- saemixModel(model = lin.model,
                             description = "linear", 
                             psi0 = matrix(c(3.06131, -0.01979), 
                                           ncol = 2, 
                                           byrow = TRUE, 
                                           dimnames = list(NULL,
                                                           c("int","slope"))),
                             covariate.model = matrix(c(1, 1, 1, 1), ncol = 2, byrow = TRUE))

opt <- list(seed = 94352514, save = FALSE, save.graphs = FALSE)

HUBU.fit.lin <- saemix(HUBU.modellin, HUBU.data, opt)

AIC(HUBU.fit) - AIC(HUBU.fit.lin)

##########################################################
# Comparison to linear model for pericalcarine

# Order data for saemix
HUBU_subPC <- HUBUThickness[,c(36,37,20,38,40)]

HUBU_orderedPC <-HUBU_subPC[order(HUBU_subPC[,1], HUBU_subPC[,2], HUBU_subPC[,3], HUBU_subPC[,4]),]

HUBU.dataPC <- saemixData(name.data       = HUBU_orderedPC,
                        name.group      = "subject_id",
                        name.predictors = "age_mri",
                        name.response   = "pericalcarine",
                        name.covariates = c("sex","euler"))

# Set up the saemix model
HUBU.modelPC <- saemixModel(model = logistic4.model,
                          description = "4-param. logistic", 
                          psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                        ncol = 4, 
                                        byrow = TRUE, 
                                        dimnames = list(NULL, 
                                                        c("lower.asymp","upper.asymp","inflec","hill"))),
                          omega.init = diag(rep(0.5, 4)),
                          covariate.model = matrix(c(1, 1, 1, 0, 1, 0, 0, 0), ncol = 4, byrow = TRUE))


HUBU.fitPC <- saemix(HUBU.modelPC, HUBU.dataPC, opt)

HUBU.modellinPC <- saemixModel(model = lin.model,
                             description = "linear", 
                             psi0 = matrix(c(3.06131, -0.01979), 
                                           ncol = 2, 
                                           byrow = TRUE, 
                                           dimnames = list(NULL,
                                                           c("int","slope"))),
                             covariate.model = matrix(c(1, 1, 1, 1), ncol = 2, byrow = TRUE))

HUBU.fit.linPC <- saemix(HUBU.modellinPC, HUBU.dataPC, opt)

AIC(HUBU.fitPC) - AIC(HUBU.fit.linPC)

##########################################################
# Run non-linear mixed effects models for each region

# results    <- list()
# inflection <- data.frame(start =factor(1:90))
# lowerAsymp <- data.frame(start =factor(1:90))
# upperAsymp <- data.frame(start =factor(1:90))
# hill       <- data.frame(start =factor(1:90))
# 
# for (r in 1:35){
# # Order data for saemix
# HUBU_sub.r <- HUBUThickness[,c(36,37,r,38,40)]
# HUBU_ordered.r <- HUBU_sub.r[order(HUBU_sub[,1], HUBU_sub.r[,2], HUBU_sub.r[,3], HUBU_sub.r[,4]),]
# 
# # Read data in
# HUBU.data.r <- saemixData(name.data       = HUBU_ordered.r,
#                         name.group      = "subject_id",
#                         name.predictors = "age_mri",
#                         name.response   = colnames(HUBUThickness[r]),
#                         name.covariates = c("sex","euler"))
# 
# HUBU.model.r <- saemixModel(model = logistic4.model,
#                           description = "4-param. logistic",
#                           psi0 = matrix(c(2.6, 2.9, 15.3, -6.8),
#                                         ncol = 4,
#                                         byrow = TRUE,
#                                         dimnames = list(NULL, c("lower.asymp","upper.asymp","inflec","hill"))),
#                           omega.init = diag(rep(0.5,4)),
#                           covariate.model = matrix(c(1,1,1,0,1,0,0,0), ncol = 4, byrow = TRUE))
# 
# HUBU.fit.r <- saemix(HUBU.model.r, HUBU.data.r, opt)
# 
# # Get individual estimates
# est <- psi(HUBU.fit.r)
# 
# temp = list(data = HUBU.data.r, model = HUBU.fit.r, estimates = est)
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
# save(results,    file="NLMM_cov_results_May22.RData")
# save(inflection, file="NLMM_cov_inflection_May22.RData")
# save(lowerAsymp, file="NLMM_cov_lowerAsymp_May22.RData")
# save(upperAsymp, file="NLMM_cov_upperAsymp_May22.RData")
# save(hill,       file="NLMM_cov_hill_May22.RData")

load("NLMM_cov_results_May22.RData")
load("NLMM_cov_inflection_May22.RData")
load("NLMM_cov_lowerAsymp_May22.RData")
load("NLMM_cov_upperAsymp_May22.RData")
load("NLMM_cov_hill_May22.RData")

colnames(inflection) <- names(results[])

# Check CV (Coefficient of variation, should be < 20% https://www.scielo.br/scielo.php?script=sci_arttext&pid=S0103-84782013000600003&lng=en&tlng=en)

cv = c()
for (r in 1:35){
  x = summary(results[[r]][["model"]])
  y = t(data.frame(x[["fixed.effects"]][["CV(%)"]]))
  cv= rbind(cv,y)
}
rownames(cv) <- colnames(inflection[1:35])
write.csv(cv,"cv_May22.csv")


estimates = c()
for (r in 1:35){
  x = summary(results[[r]][["model"]])
  y = t(data.frame(x[["fixed.effects"]][["Estimate"]]))
  estimates= rbind(estimates,y)
}
rownames(estimates) <- colnames(inflection[1:35])
write.csv(estimates,"estimates_May22.csv")

# Check correlations across regions
cors = data.frame()
for (c in 1:35){
  x1 = inflection[c]
  y1 = upperAsymp[c]
  cors[c,1] = cor(x1,y1)[1]

  x2 = inflection[c]
  y2 = lowerAsymp[c]
  cors[c,2] = cor(x2,y2)[1]

  x3 = inflection[c]
  y3 = hill[c]
  cors[c,3] = cor(x3,y3)[1]
}
rownames(cors) <- colnames(inflection[1:35])
write.csv(cors,"cors_May22.csv")

##########################################################
# EFA
fa.parallel(inflection[,c(1:4,6:11,13:19,21:31,34)])
fit_efa_i <- factanal(inflection[,c(1:4,6:11,13:19,21:31,34)], factors = 2, rotation="varimax", )
print(fit_efa_i, cutoff= 0.6)

loadings = data.frame(unclass(fit_efa_i$loadings))
loadings$region = colnames(HUBUThickness[,c(1:4,6:11,13:19,21:31,34)])
write.csv(loadings,"loadings_May22.csv")

# EFA excluding rACC and FP
fa.parallel(inflection[,c(1:4,6:11,13:19,21:30,34)])
fit_efa_i2 <- factanal(inflection[,c(1:4,6:11,13:19,21:30,34)], factors= 2, rotation="varimax")
print(fit_efa_i2, cutoff= 0.6)

loadings2 = data.frame(unclass(fit_efa_i2$loadings))
loadings2$region = colnames(HUBUThickness[,c(1:4,6:11,13:19,21:30,34)])
write.csv(loadings2,"loadings2_May22.csv")

##########################################################
# Plot regional effects for MCT

inflection$subject_id = unique(HUBUThickness$subject_id)

InflectionDataPlot = melt(inflection[,c(1:4,6:11,13:19,21:31,34,36)], id.vars = c('subject_id'))

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

sumld <- ddply(InflectionDataPlot, ~variable, summarise, mean = mean(value), median = median(value), lower = lb(value), upper = ub(value))

# get estimates from saemix
sumldmeanest = merge(sumld, estimates, by.x = "variable", by.y =0, all.x = T)
sumldmeanest$mean = sumldmeanest$V6
sumldmeanest <- sumldmeanest[order(sumldmeanest$V6),] 

inflecNew = merge(sumldmeanest, InflectionDataPlot, by='variable')
inflecNew <- inflecNew[order(inflecNew$V6),] 
inflecNew$V6[duplicated(inflecNew$V6)] <- NA

inflecNew$split = c(rep(1,1350),rep(2,1260))
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
  ggplot(aes(y = value, x = reorder(variable, mean), fill = reorder(variable, mean), label=round(V6, digits = 2))) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = value, color = reorder(variable, mean)), position = position_jitter(width = .15),  size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  theme_bw() +
  raincloud_theme + 
  coord_flip()+
  scale_color_viridis(discrete=TRUE, option = 'inferno', begin = 0.5, end = 0.9, direction = -1) +
  scale_fill_viridis(discrete=TRUE, option = 'inferno', begin = 0.5, end = 0.9, direction = -1) +
  geom_text(position=position_nudge(x = -0.3, y = 1), check_overlap = TRUE)+
  labs(y='MCT (years)', x='Brain region') 
cloud1

cloud2 =
  inflecNew2 %>%
  ggplot(aes(y = value, x = reorder(variable, mean), fill = reorder(variable, mean), label=round(V6, digits = 2))) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = value, color = reorder(variable, mean)), position = position_jitter(width = .15),  size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  theme_bw() +
  raincloud_theme + 
  coord_flip()+
  scale_color_viridis(discrete=TRUE, option = 'inferno', begin = 0, end = 0.5, direction = -1) +
  scale_fill_viridis(discrete=TRUE, option = 'inferno', begin = 0, end = 0.5, direction = -1) +
  geom_text(position=position_nudge(x = -0.3, y = 1), check_overlap = TRUE)+
  labs(y='MCT (years)', x='Brain region') 
cloud2

patchwork <- (cloud1 + cloud2) 
patchwork
ggsave("HUBU_InflectionRegions_May22.tiff", width= 10, height = 15)


##########################################################
# Plot factor loadings across regions

factor_regions1 = loadings$Factor1
factor_regions1[factor_regions1 < 0.6] = NA
factor_regions2 = loadings$Factor2
factor_regions2[factor_regions2 < 0.6] = NA

region= c(
  "bankssts"                  , "caudal anterior cingulate",
  "caudal middle frontal"     , "cuneus"                   ,
  "fusiform"                  ,  
  "inferior parietal"         , "inferior temporal"        ,  
  "isthmus cingulate"         , "lateral occipital"        ,  
  "lateral orbitofrontal"     ,  
  "medial orbitofrontal"      , "middle temporal"          , 
  "parahippocampal"           , "paracentral"              ,
  "pars opercularis"          , "pars orbitalis"           ,
  "pars triangularis"         , 
  "postcentral"               , "posterior cingulate"      ,  
  "precentral"                , "precuneus"                ,  
  "rostral anterior cingulate", "rostral middle frontal"   ,  
  "superior frontal"          , "superior parietal"        ,  
  "superior temporal"         , "supramarginal"            ,   
  "frontal pole"              ,
  "insula")  

brain_data_inflec1 = data.frame(region, factor_regions1,
                               stringsAsFactors = FALSE)
rownames(brain_data_inflec1) = NULL

brain1 =
  ggseg(.data=brain_data_inflec1, mapping=aes(fill=factor_regions1), 
        colour="grey", hemisphere = "left") + 
  labs(title="", fill="Loading F1") +
  scale_fill_viridis_c(option="mako", direction = -1, begin = 0.75) 
brain1 

brain_data_inflec2 = data.frame(region, factor_regions2,
                                stringsAsFactors = FALSE)
rownames(brain_data_inflec2) = NULL

brain2 =
  ggseg(.data=brain_data_inflec2, mapping=aes(fill=factor_regions2), 
        colour="grey", hemisphere = "left") + 
  labs(title="", fill="Loading F2") +
  scale_fill_viridis_c(option="mako", direction = -1, end = 0.5) 
brain2 

patchwork <- brain1 + brain2 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20, face="bold"))
patchwork
ggsave("HUBU_LoadingsRegions_May22.pdf", width= 10, height = 3.5)

##########################################################
# Plot parameters across brain

estimatesdf = data.frame(estimates)
inflec_average     = estimatesdf$X6
hill_average       = estimatesdf$X8
lowerAsymp_average = estimatesdf$X1
upperAsymp_average = estimatesdf$X4
regions_name       = colnames(inflection[-36])

param_average = data.frame(regions_name, inflec_average, hill_average, lowerAsymp_average, upperAsymp_average)

inflec_regions = param_average$inflec_average[c(1:4,6:11,13:19,21:31,34)]

brain_data_inflec_av = data.frame(region, inflec_regions,
                               stringsAsFactors = FALSE)
rownames(brain_data_inflec_av) = NULL

brain_av_i =
  ggseg(.data=brain_data_inflec_av, mapping=aes(fill=inflec_regions), 
        colour="grey", hemisphere = "left") + 
  labs(title="", fill="MCT") +
  scale_fill_viridis_c(option="rocket", direction = -1) 
brain_av_i

hill_regions = param_average$hill_average[c(1:4,6:11,13:19,21:31,34)]

brain_data_hill_av = data.frame(region, hill_regions,
                                  stringsAsFactors = FALSE)
rownames(brain_data_hill_av) = NULL

brain_av_hill =
  ggseg(.data=brain_data_hill_av, mapping=aes(fill=hill_regions), 
        colour="grey", hemisphere = "left") + 
  labs(title="", fill="Hill") +
  scale_fill_viridis_c(option="rocket", direction = 1) 
brain_av_hill

lowerAsymp_regions = param_average$lowerAsymp_average[c(1:4,6:11,13:19,21:31,34)]

brain_data_lowerAsymp_av = data.frame(region, lowerAsymp_regions,
                                stringsAsFactors = FALSE)
rownames(brain_data_lowerAsymp_av) = NULL

brain_av_lowerAsymp =
  ggseg(.data=brain_data_lowerAsymp_av, mapping=aes(fill=lowerAsymp_regions), 
        colour="grey", hemisphere = "left") + 
  labs(title="", fill="Lower Asymptote") +
  scale_fill_viridis_c(option="rocket", direction = -1) 
brain_av_lowerAsymp


upperAsymp_regions = param_average$upperAsymp_average[c(1:4,6:11,13:19,21:31,34)]

brain_data_upperAsymp_av = data.frame(region, upperAsymp_regions,
                                stringsAsFactors = FALSE)
rownames(brain_data_upperAsymp_av) = NULL

brain_av_upperAsymp =
  ggseg(.data=brain_data_upperAsymp_av, mapping=aes(fill=upperAsymp_regions), 
        colour="grey", hemisphere = "left") + 
  labs(title="", fill="Upper Asymptote") +
  scale_fill_viridis_c(option="rocket", direction = -1) 
brain_av_upperAsymp

patchwork_av <- (brain_av_i + brain_av_hill) /  (brain_av_upperAsymp + brain_av_lowerAsymp)+ plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20, face="bold"))
patchwork_av
ggsave("HUBU_ParamRegions_May22.pdf", width= 10, height = 7)


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
+                  parahippocampal 
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
+                  parahippocampal 

                bankssts ~ 14.36
 caudalanteriorcingulate ~ 14.36 
     caudalmiddlefrontal ~ 14.36 
                  cuneus ~ 14.36 
                fusiform ~ 14.36 
        inferiorparietal ~ 14.36 
        inferiortemporal ~ 14.36 
        isthmuscingulate ~ 14.36 
        lateraloccipital ~ 14.36 
    lateralorbitofrontal ~ 14.36 
     medialorbitofrontal ~ 14.36 
          middletemporal ~ 14.36 
             paracentral ~ 14.36 
         parsopercularis ~ 14.36 
           parsorbitalis ~ 14.36 
        parstriangularis ~ 14.36 
             postcentral ~ 14.36 
      posteriorcingulate ~ 14.36 
              precentral ~ 14.36 
               precuneus ~ 14.36 
rostralanteriorcingulate ~ 14.36 
    rostralmiddlefrontal ~ 14.36 
         superiorfrontal ~ 14.36 
        superiorparietal ~ 14.36 
        superiortemporal ~ 14.36 
           supramarginal ~ 14.36 
             frontalpole ~ 14.36 
                  insula ~ 14.36
         parahippocampal ~ 14.36
  '

fit_OneFactor_i <- cfa(model_OneFactor, data=inflection, missing="fiml", estimator = "mlr")
fit_OneFactor_intfixed_i <- cfa(model_OneFactor_intfixed_i, data=inflection, missing="fiml", estimator = "mlr")

summary(fit_OneFactor_i, fit =T)

anova(fit_OneFactor_intfixed_i, fit_OneFactor_i)

##########################################################
# Gradient

require(brainGraph)
data("dk")

dk$area <- c(
  "bank of the superior temporal sulcus"
  ,"caudal anterior cingulate"           
  ,"caudal middle frontal gyrus"         
  ,"cuneus"                              
  ,"entorhinal"                          
  ,"fusiform"                            
  ,"inferior parietal lobule"            
  ,"inferior temporal gyrus"             
  ,"isthmus cingulate cortex"            
  , "lateral occipital gyrus"             
  , "lateral orbitofrontal"               
  , "lingual"                             
  , "medial orbitofrontal"                
  , "middle temporal gyrus"               
  , "parahippocampal"                     
  , "paracentral"                         
  , "pars opercularis"                    
  , "pars orbitalis"                      
  , "pars triangularis"                   
  , "pericalcarine"                       
  , "postcentral"                         
  , "posterior cingulate cortex"          
  , "precentral"                          
  , "precuneus"                           
  , "rostral anterior cingulate cortex"   
  , "rostral middle frontal gyrus"        
  , "superior frontal gyrus"              
  , "superior parietal lobule"            
  , "superior temporal gyrus"             
  , "supramarginal gyrus"                 
  , "frontal pole"                        
  , "temporal pole"                       
  , "transverse temporal"                 
  , "insula"                              
  , "bank of the superior temporal sulcus"
  , "caudal anterior cingulate"           
  , "caudal middle frontal gyrus"         
  , "cuneus"                              
  , "entorhinal"                          
  , "fusiform"                            
  , "inferior parietal lobule"            
  , "inferior temporal gyrus"             
  , "isthmus cingulate cortex"            
  , "lateral occipital gyrus"             
  , "lateral orbitofrontal"               
  , "lingual"                             
  , "medial orbitofrontal"                
  , "middle temporal gyrus"               
  , "parahippocampal"                     
  , "paracentral"                         
  , "pars opercularis"                    
  , "pars orbitalis"                      
  , "pars triangularis"                   
  , "pericalcarine"                       
  , "postcentral"                         
  , "posterior cingulate cortex"          
  , "precentral"                          
  , "precuneus"                           
  , "rostral anterior cingulate cortex"   
  , "rostral middle frontal gyrus"        
  , "superior frontal gyrus"              
  , "superior parietal lobule"            
  , "superior temporal gyrus"             
  , "supramarginal gyrus"                 
  , "frontal pole"                        
  , "temporal pole"                       
  , "transverse temporal"                 
  , "insula"
)        

mni_y_av <- aggregate(dk$y.mni, list(dk$area), mean) #average the y coordinate for both hemispheres 
colnames(mni_y_av) <- c("area", "mni.y")

mni_y_av[-2] <- data.frame(lapply(mni_y_av[-2], gsub, pattern = " gyrus", replacement = "", fixed = TRUE))
mni_y_av[-2] <- data.frame(lapply(mni_y_av[-2], gsub, pattern = " lobule", replacement = "", fixed = TRUE))
mni_y_av[-2] <- data.frame(lapply(mni_y_av[-2], gsub, pattern = "bank of the superior temporal sulcus", replacement = "bankssts", fixed = TRUE))
mni_y_av[-2] <- data.frame(lapply(mni_y_av[-2], gsub, pattern = " cortex", replacement = "", fixed = TRUE))
mni_y_av[-2] <- data.frame(lapply(mni_y_av[-2], gsub, pattern = " ", replacement = "", fixed = TRUE))

postant = merge(sumld, mni_y_av, by.x = 'variable', by.y = 'area', all.x = T)

cor.test(postant$mean, postant$mni.y)

regressionBF(
  formula = mean ~ mni.y,
  data = postant
)


mni_z_av <- aggregate(dk$z.mni, list(dk$area), mean) #average the z coordinate for both hemispheres 
colnames(mni_z_av) <- c("area", "mni.z")

mni_z_av[-2] <- data.frame(lapply(mni_z_av[-2], gsub, pattern = " gzrus", replacement = "", fixed = TRUE))
mni_z_av[-2] <- data.frame(lapply(mni_z_av[-2], gsub, pattern = " lobule", replacement = "", fixed = TRUE))
mni_z_av[-2] <- data.frame(lapply(mni_z_av[-2], gsub, pattern = "bank of the superior temporal sulcus", replacement = "bankssts", fixed = TRUE))
mni_z_av[-2] <- data.frame(lapply(mni_z_av[-2], gsub, pattern = " cortex", replacement = "", fixed = TRUE))
mni_z_av[-2] <- data.frame(lapply(mni_z_av[-2], gsub, pattern = " ", replacement = "", fixed = TRUE))

dorsvent = merge(sumld, mni_z_av, by.x = 'variable', by.y = 'area', all.x = T)

cor.test(dorsvent$mean, dorsvent$mni.z)

regressionBF(
  formula = mean ~ mni.z,
  data = na.omit(dorsvent)
)

##########################################################
# Replicate main model with min T=3

HUBU3 = merge(HUBUThickness, sumdata, by = "subject_id", all =T)

HUBU3 = HUBU3[which(HUBU3$mr_round.max > 2), ]

length(unique(HUBU3$subject_id)) # 84

# Order data for saemix
HUBU3_sub <- HUBU3[,c(1,37,36,38,40)]
HUBU3_ordered <- HUBU3_sub[order(HUBU3_sub[,1], HUBU3_sub[,2], HUBU3_sub[,3], HUBU3_sub[,4]),]

# Read data in
HUBU3.data <- saemixData(name.data       = HUBU3_ordered,
                         name.group      = "subject_id",
                         name.predictors = "age_mri",
                         name.response   = "MeanThickness",
                         name.covariates = c("sex","euler"))


# Set up the saemix model
HUBU3.model <- saemixModel(model = logistic4.model,
                           description = "4-param. logistic", 
                           psi0 = matrix(c(2.6, 2.9, 15.4, -7.2), 
                                         ncol = 4, 
                                         byrow = TRUE, 
                                         dimnames = list(NULL, 
                                                         c("lower.asymp","upper.asymp","inflec","hill"))),
                           omega.init = diag(rep(0.5, 4)),
                           covariate.model = matrix(c(1, 1, 1, 0, 1, 0, 0, 0), ncol = 4, byrow = TRUE))

HUBU3.fit <- saemix(HUBU3.model, HUBU3.data, opt)

