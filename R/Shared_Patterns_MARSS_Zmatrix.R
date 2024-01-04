#--------- MARSS SHARED SEASONAL PATTERNS Time Series Analysis ----------

# This script tests for coherence among response variables using multivariate autoregressive state space models. This work relies heavily on the the MARSS package written by E. E. Holmes, M. D. Scheuerell, and E. J. Ward. Please find links to the package user guide and online book below.

# MARSS User guide: https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf

# MARSS Online book: https://atsa-es.github.io/atsa-labs/chap-mss.html

#--------- SECTION 1: SETTING UP R AND LOADING NECESSARY PACKAGES ----------

# install and/or load the necessary R packages

if("pacman" %in% installed.packages() == FALSE){install.packages("pacman")}

pacman::p_load(MARSS, here, tidyverse, reshape2, zoo, broom, janitor) 

#--------- SECTION 2: READ IN AND PREP RESPONSE VARIABLE  ----------

# Read data
data <- read_csv(file.path(here("data", "model_data", "MARSS_NO3response_data.csv")))

# Select response variable 
yy <- data %>% select(time,
                      NO3yield_km2_composite,#Change this per response variable
                     HEELCode)

# Transposing and zscore matrix
dat <- reshape2::acast(yy, HEELCode ~ time, value.var = "NO3yield_km2_composite")#Change this per response variable

#Change this per response variable
dat.log <- log(dat) #NO3mgL/NO3yield
#dat.log <- dat #d15N/d18O 

#Z-score the data
dat.z <- zscore(dat.log)


# SECTION 3: SET UP MARSS MODEL MATRICES -------------------------------------

# Here we evaluated the data support for the following hypotheses about PS nitrate river trends

# Each Z model is a hypothesis
Z.models <- list(
  #one hidden state/trend
  H1 = matrix(1,13,1),
  #trends are defined by watershed
  H2 = factor(c("cedar","green", "nooksack",
                'samish', "skagit", "still",
                "snoho", "duck", "puyallup",
                "nisqually","skoko", "deschutes", "elwha")), 
  # trends defined by mountain range
  H3 = factor(c("cascades","cascades", "cascades",
                'cascades',"cascades", "cascades","cascades", 
                "olympics", "cascades","cascades",
                "olympics","cascades", "olympics")),
  # trends defined by lulc
  H4 = factor(c("urban","urban", "crop",
                'ag', "forest", "forest",
                "forest", "forest", "forest",
                "ag","forest", "ag", "forest")),
  # trends defined by soil/geomorphic characteristics
  H5 = factor(c("midlow","midlow", "midlow",
                'midlow', "large", "midlow",
                "large", "elevslope", "large",
                "midlow","elevslope", "midlow", "elevslope"))) 

names(Z.models) <- c("basin","watershed",'mountain',"landuse","soils_geo")


mod.list = list(
  B = "identity",
  U = "zero",
  A = "zero",
  R = "diagonal and equal",
  x0 = "unequal",
  tinitx = 0 )

Q.models2 <- c("diagonal and equal", "diagonal and unequal","equalvarcov","unconstrained")


# SECTION 4: RUN MARSS MODELS ---------------------------------------------

out.tab <- NULL
fits <- list()
for(i in 1:length(Z.models)){
  
  for(Q.model in Q.models2){
    fit.model = c(list(Z=Z.models[[i]], Q=Q.model), mod.list)
    fit = MARSS(dat.z, model=fit.model,
                silent=TRUE, control=list(maxit=6000))
    out=data.frame(model_Z=names(Z.models)[i], 
                   trend_num=length(unique(Z.models[[i]])),
                   Q=Q.model,
                   logLik=fit$logLik, AICc=fit$AICc, 
                   num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence,
                   stringsAsFactors = FALSE)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
    
  }
}

min.AICc <- order(out.tab$AICc)
out.tab.1 <- out.tab[min.AICc, ]%>%
  cbind(.,delta.AICc = out.tab.1$AICc - out.tab.1$AICc[1])
out.tab.1


# SECTION 6: EXPORT MODEL DATA --------------------------------------------

# Export supplementary AIC table
write_csv(out.tab.1 %>% select(-c(num.iter, converged)), 
          here("data/model_output/shared_patterns/NO3yield_Z_AIC.csv"))# Change this per response variable

# Tidy and export best model states for Figure 3

## best model 
best_i <- which(out.tab[,"AICc"] == min(out.tab[,"AICc"]))
best_fit <- fits[[best_i]]

years = c(1:24)
par(mfrow = c(1, 2))
for (i in 1:2) {
  plot(years, best_fit$states[i, ], ylab = "best fit", 
       xlab = "", type = "l")
  lines(years, best_fit$states[i, ] - 1.96 * best_fit$states.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  lines(years, best_fit$states[i, ] + 1.96 * best_fit$states.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  title(rownames(best_fit$states)[i])
}

## check acf
plot(best_fit, plot.type = 'acf.std.model.resids.ytt1')

## Tidy states for Figure 3
mtn_states <- bind_cols(as_tibble(t(best_fit$states)), as_tibble(t(best_fit$states.se)))%>%
  rename(cascades_fit = cascades...1, olympics_fit = olympics...2, cascades_se = cascades...3, olympics_se = olympics...4)%>%
  mutate(date = seq.Date(as.Date("2014-10-01"),as.Date("2016-09-01"),"month"))

## EXPORT STATES FOR FIGURE 3
write_csv(mtn_states, file.path(here("data/model_output/shared_patterns/NO3yield_mtn_states.csv")))# Change this per response variable

