#______________________________________________________________________________#
#--------- MARSS HYDROCLIMATE EFFECT COVARIATE Time Series Analysis ----------
#______________________________________________________________________________#

# This script tests for the effects of hydroclimate on river NO3 dynamics using multivariate autoregressive state space models. This work relies heavily on the the MARSS package written by E. E. Holmes, M. D. Scheuerell, and E. J. Ward. Please find links to the package user guide and online book below.

# MARSS User guide: https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf

# MARSS Online book: https://atsa-es.github.io/atsa-labs/chap-mss.html

#______________________________________________________________________________#
# SECTION 1: SETTING UP R AND LOADING NECESSARY PACKAGES ------------------
#______________________________________________________________________________#

# install and/or load the necessary R packages
if("pacman" %in% installed.packages() == FALSE){install.packages("pacman")}

pacman::p_load(MARSS, here, tidyverse, reshape2, zoo, broom, janitor) 

## set directory locations
nitrate <- here("data")

#______________________________________________________________________________#
# SECTION 2: READ IN AND PREP RESPONSE VARIABLE ---------------------------
#______________________________________________________________________________#

#Read data
data <- read_csv(file.path(nitrate, "model_data/MARSS_NO3response_data.csv"))

#Select columns
yy <-data %>% dplyr::select(time,
                            NO3yield_km2_composite,#Change this
                            HEELCode)

#Transposing and zscore matrix
dat <- reshape2::acast(yy, HEELCode ~ time, value.var = "NO3yield_km2_composite")#Change this

#Change this 
dat.log <- log(dat) #NO3mgL/NO3yield
#dat.log <- dat #d15N/d18O (bc neg values)

#Z-score the data
dat.z <- zscore(dat.log)

#______________________________________________________________________________#
#--------- SECTION 3: READ IN AND PREP COVARIATE DATA ----------
#______________________________________________________________________________#

# Read data
covar <- read_csv(file.path(nitrate, "model_data/MARSS_CLIMcovar_data.csv"))

# Interpolate a few missing values 
covariates <-covar %>%
  group_by(HEELCode) %>%
  arrange(HEELCode, time) %>%
  mutate(snowmelt = na.approx(snowmelt, maxgap=4, rule =2),
         water_temp = na.approx(water_temp, maxgap=4, rule =2))

#______________________________________________________________________________#
# SECTION 4: SET UP LOWER C MATRICES --------------------------------------
#______________________________________________________________________________#

#Transpose each column into 13x24 matrix, river site by month
covar.name <- colnames(covariates[,-c(1:6)])
all.cast <- list()
for (i in 1:length(covar.name)) {
  all.cast[[i]] <- acast(covariates, HEELCode ~ time, value.var = covar.name[[i]])
}
names(all.cast) <- covar.name

# zscore each temporal matrix
all.cast.z <- 
  map(all.cast, ~.x %>% zscore())

#### STOP COMMENT/UNCOMMENT BELOW per response variable 

# Lower c matrices for NO3 yield (doesn't include discharge metrics)
c.matrices <- all.cast.z[c(1:4)]

# Lower c matrices for NO3mgl, d15N, d18O
#c.matrices <- all.cast.z

#______________________________________________________________________________#
# SECTION 5: SET UP UPPER C MATRICES --------------------------------------
#______________________________________________________________________________#

watershed <- c("A","B","C","E","F","G","H","I","J","K","L","M","Z")
mountain <- c("cascades","cascades", "cascades",
              'cascades',"cascades", "cascades","cascades", 
              "olympics", "cascades","cascades",
              "olympics","cascades", "olympics")
lulc <- c("urban","urban", "ag",
          'ag', "forest", "forest",
          "forest", "forest", "forest",
          "ag","forest", "ag", "forest")
geomorph <- c("midlow","midlow", "midlow",
              'midlow', "large", "midlow",
              "large", "elevslope", "large",
              "midlow","elevslope", "midlow", "elevslope")

#### STOP COMMENT/UNCOMMENT per response variable

# YIELD
C.names <- covar.name[1:4]

# NO3mgL, d15N, d18O
#C.names <- covar.name

# Create matrices
basin.list <- watershed.list <- mountain.list <- lulc.list <- geomorph.list <- list()
for(i in 1:length(C.names)){
  
  # create matrices
  basin.list[[i]] <- matrix(list(0), nrow = 13, ncol = 13)
  watershed.list[[i]] <- matrix(list(0), nrow = 13, ncol = 13)
  mountain.list[[i]] <- matrix(list(0),13,13)
  lulc.list[[i]] <- matrix(list(0),13,13)
  geomorph.list[[i]] <- matrix(list(0),13,13)
  
  # change the diag
  diag(basin.list[[i]]) <- C.names[i]
  diag(watershed.list[[i]]) <- paste0(watershed,'.',C.names[i])
  diag(mountain.list[[i]]) <- paste0(mountain,".",C.names[i])
  diag(lulc.list[[i]]) <- paste0(lulc,".",C.names[i])
  diag(geomorph.list[[i]]) <- paste0(geomorph,".",C.names[i])
  
}

# name each matrix
names(basin.list) <- C.names
names(watershed.list) <- paste0(C.names,"_u")
names(mountain.list) <- paste0(C.names,"_mtn")
names(lulc.list) <- paste0(C.names,"_lulc")
names(geomorph.list) <- paste0(C.names,"_geo")

C.matrices <- c(basin.list, watershed.list, mountain.list, lulc.list, geomorph.list)

#______________________________________________________________________________#
# SECTION 6: MARSS MODEL MATRICES --------------------------------------------
#______________________________________________________________________________#

#Add a null model (no covariates)
null <- list(null = "zero")

# Rep c matrices 5 times for all scales
c.models <- c(null, c.matrices, c.matrices, c.matrices, c.matrices, c.matrices)

# C matrices
C.models <- c(null, C.matrices)

mod.list = list(
  B = "identity",
  U = "zero",
  Z = "identity",
  A = "zero",
  R = "diagonal and equal",
  x0 = "unequal",
  tinitx = 0 )

Q.models2 <- c('diagonal and unequal', "diagonal and equal","equalvarcov","unconstrained")

#______________________________________________________________________________#
#--------- SECTION 7: RUN MARSS COVARIATE MODELS ----------
#______________________________________________________________________________#

# This will take awhile.

out.tab <- NULL
fits <- list()
for(i in 1:length(C.models)){
  
  for(Q.model in Q.models2){
    fit.model = c(list(C=C.models[[i]], c = c.models[[i]],Q=Q.model), mod.list)
    fit = MARSS(dat.z, model=fit.model,
                silent=TRUE, control=list(maxit=5000))
    out=data.frame(mod_covar=names(C.models)[i], 
                   #effect_num=length(unique(diag(C.models[[i]]))),
                   Q=Q.model,
                   logLik=fit$logLik, 
                   AICc=fit$AICc, 
                   num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence,
                   stringsAsFactors = FALSE)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
    
  }
}

min.AICc <- order(out.tab$AICc)
out.tab.1 <- out.tab[min.AICc, ]
out.tab.1 <- cbind(out.tab.1,delta.AICc = out.tab.1$AICc - out.tab.1$AICc[1])
out.tab.1


#______________________________________________________________________________#
#--------- SECTION 8: CREATE AND WRITE MODEL OUTPUTS ----------
#______________________________________________________________________________#

# Export supplementary AIC table
write_csv(out.tab.1 %>% select(-c(num.iter, converged)), 
          here("data/model_output/hydroclimate_covariates/NO3yield_C_AIC.csv"))# Change this per response variable

# Tidy and export best model states for Figure 3

## best model 
best_i <- which(out.tab[,"AICc"] == min(out.tab[,"AICc"]))
best_fit <- fits[[best_i]]
 
## check acf
plot(best_fit, plot.type = 'acf.std.model.resids.ytt1')

## Tidy coefficients for Figure 4/5
names <- make_clean_names(paste(out.tab$mod_covar,out.tab$Q), replace = c('equalvarcov' = 'eqvc'))
names(fits) <- names

fig_fits <- fits[c("ppt_mmmtn_eqvc", "water_tempmtn_eqvc",
                   "spi_12mtn_eqvc", "snowmeltmtn_eqvc",
                   "ppt_mmu_eqvc")]

output_coef <- here("data/model_output/hydroclimate_covariates/covar_coef")
tidyCI <- list()
for (i in 1:length(fig_fits)) {
  tidyCI[[i]] <- MARSSparamCIs(fig_fits[[i]], method='hessian')
  write_csv(tidy(tidyCI[[i]]), file.path(output_coef, paste0(names(fig_fits)[i],"_YIELD.csv")))
}


