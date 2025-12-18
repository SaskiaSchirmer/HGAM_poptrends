########################################
### 0. load packages
########################################
library(tidyverse)
library(data.table) # read in simulated data from directories
library(ggplot2)
library(poptrend)
library(gratia)
library(mgcv)

########################################
### 1. Simulation study
########################################

#####################################
## 1.1 Simulate artificial data sets
#####################################

# packages
library(tidyverse)

# name of working directory
dir <- "~/simulation_study/simulated_data/"


# set parameters
oT <- 30  # observation time
n_loc = 5 # number of locations
n_per_loc <- 20 # number of locations per group
# number of subpopulations: n_loc * n_per_loc
N0 <- 1000 # population size at start
niter <- 250 # number of iterations per scenario

# possible scenarios
setting <- expand.grid(nonlinear = c(FALSE, TRUE),
            group = c(FALSE, TRUE),
            negativebinomial = c(FALSE, TRUE))

# load functions

# subpopulation specific part of the growth rate
lambda_func <- function(loc, oT, time) {
  if(loc == 1) {
    lambda <- 0
  } else if(loc == 2) {
    lambda <- -0.05
  } else if(loc == 3) {
    lambda <- 0.05
  } else if(loc == 4) {
    lambda <- -9/oT^2*exp(-(9*(oT-2*time)^2)/(4*oT^2))*(2*time-oT)
  } else if(loc == 5) {
    lambda <- 9/oT^2*exp(-(9*(oT-2*time)^2)/(4*oT^2))*(2*time-oT)
  }
  return(lambda)
}

# common part of the growth rate
lambda_overall <- function(oT, time, type = "inc") {
  # log-linear increasing common part of the growth rate
  if(type == "inc")  lambda <- 0.1 
  # non-linear common part of the growth rate
  if(type == "high-low") {
    lambda <- -9/oT^2/2*exp(-(9*(time-oT/4)^2)/(oT^2))*(4*time-oT) 
  } 
  
  return(lambda)
}

# calculate deterministic population size from growth rate
mu_func <- function(loc, oT, N0, gr = group_structure, type = type) {
  N <- numeric(oT)
  N[1] <- N0
  
  if(gr) {
    for(i in 2:oT) {
      N[i] <- N[i-1]*(1 + lambda_func(loc, oT, i-1) + 
                        lambda_overall(oT, i-1, type = type))
    }
  } else {
    for(i in 2:oT) {
      N[i] <- N[i-1]*(1 + lambda_overall(oT, i-1, type = type))
    }
  }
  return(N)
}

# generate stochastic population size with the deterministic population size as mean
count_func <- function(oT, n_loc, n_per_loc, N0, family = "nb", size = 100, 
                       group_structure, type = type) {
  
  if(family == "nb") {
    sz <- function(geo1_id, oT, N0, gr = group_structure, type = type){
      rep(size, oT)
    }
  } else {
    sz <- function(geo1_id, oT, N0, gr = group_structure, type = type){
      mu_func(geo1_id, oT, N0, gr = group_structure, type = type)/size
    }
  }
  
  empirical <-
    tibble(
      year = rep(1:oT, n_loc*n_per_loc),
      geo1_id = rep(rep(sample(1:n_loc), each = oT), n_per_loc),
      iter = rep(1:n_per_loc, each = n_loc*oT) # iteration per location
    ) %>% 
    group_by(geo1_id, year, iter) %>% 
    mutate(
      count = stats::rnbinom(n = 1, mu = mu_func(geo1_id, oT, N0, 
                                                 gr = group_structure, 
                                                 type = type)[year], 
                             size = sz(geo1_id, oT, N0, gr = group_structure, type = type)[year]) # Generate count data that actually depends on time
    ) |> 
    mutate(location_id_fac = paste(geo1_id, iter, collapse = "_"))
  
  empirical$geo1_id <- factor(empirical$geo1_id)
  empirical$location_id_fac <- factor(empirical$location_id_fac)
  
  empirical
}

# simulate data for all scenarios

for(i in 1:nrow(setting)) {
  set <- setting[i,]
  
  if(set$nonlinear) {
    type <- "high-low"
  } else {
    type <- "inc"
  }
  if(set$group) {
    group_structure <- TRUE
  } else {
    group_structure <- FALSE
  }
  if(set$negativebinomial) {
    fam_data <- "nb"
  }else {
    fam_data <- "quasipoisson"
  }
  
  dir_name <- paste(dir, paste(setting[i,], collapse = "_"), sep = "")
  dir.create(dir_name)
  
  for(k in 1:niter) {
    df <- count_func(oT, n_loc, n_per_loc, N0, fam_data, size = 100,
                     group_structure = group_structure, type = type)
    
    write.table(df, file = paste(dir_name, "/df_", k, sep = ""), row.names = FALSE)
  }
}


###############################################
### 1.2 Overview on artificial data (Figure A1)
###############################################

tmp <- NULL

for(i in 1:nrow(setting)) {
  dir_name <-paste(dir, paste(setting[i,], collapse = "_"), sep = "")
  type_name <- paste(str_sub(setting[i,], end = 1), collapse = "")
  
  tbl_fread <- 
    list.files(dir_name, full.names = TRUE) %>% 
    map_df(~fread(.))
  
  tmp <- tmp %>% rbind(tbl_fread %>% 
                         group_by(year) %>% 
                         summarise(me = mean(count),
                                   medi = median(count),
                                   sd = sd(count), 
                                   lq = quantile(count, 0.025), 
                                   uq = quantile(count, 0.975)) %>% 
                         mutate(type = type_name))
}

tmp <- tmp %>% separate(type, c("nonlinear", "group", "negbin"), 
                        sep = 1:3, remove = F)


group_labs <- c("no group", "group")
names(group_labs) <- c("F", "T")

negbin_labs <- c("quasi-poisson", "negative binomial")
names(negbin_labs) <- c("F", "T")

nonlin_labs <- c("log-linear", "non-linear")
names(nonlin_labs) <- c("F", "T")

oT <- 30
ggplot(data = tmp %>% filter(type %in% c("FFF", "FFT", "FTF", "FTT"))) + 
  geom_line(data = data.frame(x = 1:oT, y= mu_func(1, oT, 1000, FALSE, "inc")), 
            aes(x = x, y = y, linetype = "true")) +
  geom_pointrange(aes(x = year, y = me, ymin = lq, ymax = uq, 
                      color = factor(type), shape = "mean")) +
  geom_point(aes(x= year, y = medi, color = factor(type), shape = "median")) +
  facet_grid(group~negbin, scales = "free", 
             labeller = labeller(group = group_labs, negbin = negbin_labs)) +
  scale_shape_manual(values = c(1,2), labels = c("mean", "median")) +
  scale_linetype_manual(values = 1, labels = "true") +
  labs(y = "average of simulated count data", color = "model", 
       shape = "Central tendency", linetype = "true count", x = "year")  +
  theme(text = element_text(size = 20))

ggplot(data = tmp %>% filter(type %in% c("TFF", "TFT", "TTF", "TTT"))) + 
  geom_line(data = data.frame(x = 1:oT, 
                              y= mu_func(1, oT, 1000, FALSE, "high-low")), 
            aes(x = x, y = y, linetype = "true")) +
  geom_pointrange(aes(x = year, y = me, ymin = lq, ymax = uq, 
                      color = factor(type), shape = "mean")) +
  geom_point(aes(x= year, y = medi, color = factor(type), shape = "median")) +
  facet_grid(group~negbin, scales = "free", 
             labeller = labeller(group = group_labs, negbin = negbin_labs)) +
  scale_shape_manual(values = c(1,2), labels = c("mean", "median")) +
  scale_linetype_manual(values = 1, labels = "true") +
  labs(y = "average of simulated count data", color = "model", 
       shape = "Central tendency", linetype = "true count", x = "year") +
  theme(text = element_text(size = 20))


####################################################
### 1.3 simulate models for all simulated data sets
####################################################

list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "kableExtra"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

#################################
# load functions
#################################

hgam <- function(df, fam) {
  poptrend::ptrend(count ~ trend(year, k = 8, tempRE = F, type = "smooth") +
                     s(year, location_id_fac, bs = "fs", k = 3, m=1),
                   family = fam, data = df, engine = "bam", discrete = TRUE, nthreads = 8)
}

poptrend <- function(df, family) {
  poptrend::ptrend(count ~ trend(year, type = 'smooth') + s(location_id_fac, bs = "re"),
                   family = family,
                   data = df, engine = "bam", discrete = TRUE, nthreads = 8)
}

pop_linear <- function(df, family) {
  poptrend::ptrend(count ~ trend(year, type = 'loglinear') + s(location_id_fac, bs = "re"),
                   family = family,
                   data = df, engine = "bam", discrete = TRUE, nthreads = 8)
}


lin_group <- function(df, family) {
  poptrend::ptrend(count ~ trend(year, type = "loglinear") + s(location_id_fac,year, bs = "fs", m = 1),
                   family = family,
                   data = df, engine = "bam", discrete = TRUE, nthreads = 8)
}



run_mod <- function(model_func, df, family) {
  warns <- list()
  
  trend <- NULL
  withCallingHandlers({trend <- model_func(df, family)}, 
                      warning = function(warn) {warns <<- c(warns, list(warn$message))
                      trend <- trend
                      })
  list(trend = trend, warns = warns)
}



warn_handle <- function(trend, iter) {
  # gam_mod <- trend$trend$gam
  gam_boot <- as.data.frame(trend$trend$bootTrend) %>% 
    mutate(year = trend$trend$trendFrame$year) %>% 
    pivot_longer(names_to = "b_iter", values_to = "value", 
                 cols = grep("V", colnames(.)), names_prefix = "V") %>% 
    mutate(iter = i)
  
  warns <- trend$warns
  
  if(length(warns) == 1) {
    if(grepl("Modell hat 1-d-Glättungen derselben Variable wiederholt",warns)) {
      # syear <- smooth_estimates(gam_mod) %>% mutate(iter = i)
      gam_boot <- gam_boot %>% mutate(iter = i)
    } else if(grepl("Schritt Fehlgeschlagen in der Schätzung von theta",warns)) {
      # syear <- smooth_estimates(gam_mod) %>% mutate(iter = i) %>% mutate(
      #   est = NA, se = NA
      # )
      
      gam_boot <- gam_boot %>% mutate(value = NA)
    }
  } else if(length(warns) > 1) {
    if(sum(grepl("Schritt Fehlgeschlagen in der Schätzung von theta",warns))> 0) {
      # syear <- smooth_estimates(gam_mod) %>% mutate(iter = i) %>% mutate(
      #   est = NA, se = NA
      # )
      
      gam_boot <- gam_boot %>% mutate(value = NA)
    } else {
      print("Have a look at the warning message text.")
    }
  }
  
  gam_boot
}

run_one_mod <- function(data, oT, n_loc, n_per_loc, N0, 
                        fam_model, 
                        size, model_func, i) {
  
  
  
  trend_hg_nb <- run_mod(model_func, data, fam_model)
  
  syear <- warn_handle(trend_hg_nb, i)
  
  syear
}

##########################
## set parallel settings
##########################



n.cores <- 8

model_name <- "hglm"
model_function <- lin_group
fam <- "quasipoisson"

# set parameters
oT <- 30
n_loc = 5
n_per_loc <- 20 # number of locations per group
N0 <- 1000
niter <- 250

# for every data type in 1:8 do
#for(i in 1:nrow(setting)){
i<-8

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
set <- setting[i,]

path <- paste("simulated_data/", paste(setting[i,], collapse = "_"), sep = "")
ls_files <- list.files(path, full.names = TRUE)

x <- foreach(
  i = 1:niter, 
  .combine = 'rbind'
) %dopar% {
  library(tidyverse)
  df <- read.table(list.files(path, full.names = TRUE)[i], header = TRUE)
  df$location_id_fac <- factor(df$location_id_fac)
  run_one_mod(df, oT, n_loc, n_per_loc, N0, fam, 100,
              model_function, i)
}
x


dir_path <- paste("sim_results/",paste(setting[i,], collapse = "_"), sep = "")
if(!dir.exists(dir_path)) {
  dir.create(dir_path)
}

save(x, file = paste(dir_path,
                     "/", model_name, "_", fam,".Rdata", sep = ""))

est_year_gs_qu <- x %>% group_by(iter, b_iter) %>% 
  mutate(lambda = exp((log(value)-lag(log(value), order_by = year))/(year - lag(year, order_by = year))),
         N = value) %>% 
  mutate(N = N/first(N)) %>% 
  group_by(year) %>% filter(!is.na(lambda)) %>% 
  summarise(me_lambda = mean(lambda), 
            lq_lambda = quantile(lambda, prob = 0.025),
            uq_lambda = quantile(lambda, prob = 0.975),
            sd_lambda = sd(lambda),
            me_N = mean(N), 
            lq_N = quantile(N, prob = 0.025),
            uq_N = quantile(N, prob = 0.975),
            sd_N = sd(N)) 

save(est_year_gs_qu, file = paste(dir_path,
                                  "/", model_name, "_", fam,"summary.Rdata", sep = ""))

parallel::stopCluster(cl = my.cluster)


###############################################
### 1.4 Simulation results (Figure 1)
###############################################

plot_models <- function(summary_path = "~/simulation_study/sim_results_summary/",
                        data_path = "FALSE_TRUE_FALSE",
                        true_lambda = lambda_overall, 
                        type, 
                        title = "",
                        oT,
                        labels = c("TRIM",
                                   "TRIM_cp",
                                   "FFF",
                                   "TFF",
                                   "FFT",
                                   "TFT",
                                   
                                   "FTF",
                                   "TTF",
                                   "FTT",
                                   "TTT")) {
  ls <- list.files(paste(summary_path, data_path, sep = ""), recursive = TRUE, full.names = TRUE)
  
  res <- data.frame(matrix(NA, ncol = 10, nrow = 1))
  colnames(res) <- c("year", "me_lambda", "lq_lambda", "uq_lambda", "sd_lambda", 
                     "me_N", "lq_N", "uq_N", "sd_N", "model")
  for(file in ls){
    load(file)
    
    res <- rbind(res, cbind(est_year_gs_qu, model = file))
    
  }
  
  print(levels(factor(res$model)))
  res <- res %>% slice(-1) %>% mutate(model = (factor(model, labels = labels)))
  
  
  
  
  pl_lambda <- ggplot(res) +
    geom_line(aes(x = year, y = (me_lambda), color = model), show.legend = TRUE) +
    geom_ribbon(aes(x = year, y = (me_lambda), ymin = (lq_lambda), ymax = (uq_lambda), color = model, fill = model),alpha = 0.1, show.legend = FALSE) +
    geom_line(aes(x = year, y = 1+true_lambda(max(year), year, type))) +
    labs(y = "population growth rate", 
         title = title)
  
  pl_N <- ggplot() +
    geom_line(data = res, aes(x = year, y = (me_N),  color = model), show.legend = TRUE) +
    geom_ribbon(data = res, aes(x = year, y = (me_N), ymin = (lq_N), ymax = (uq_N), 
                                color = model, fill = model), alpha = 0.1, 
                show.legend = FALSE, linetype = 2) +
    geom_line(data = data.frame(x = 1:oT, y= mu_func(1, oT, 1, FALSE, type)), aes(x = x, y = y)) +
    labs(y = "index of population size", 
         title = title)
  list(pl_lambda, pl_N)
}



pl0 <- plot_models(data_path = "FALSE_FALSE_FALSE",
                   title = "FFF",
                   type = "inc", oT = 30)

pl0[[1]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))# +
ylim(1.08,1.13)
ggsave("FFF.png", width = 30, height = 4.5, units = "cm")
pl0[[2]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(0,25)
ggsave("FFF_n.png", width = 30, height = 4.5, units = "cm")

pl <- plot_models(data_path = "FALSE_TRUE_FALSE",
                  title = "subpopulation structure / log-linear / quasi-poisson",
                  type = "inc", oT = 30)

pl[[1]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 20))#+
ylim(1.08,1.13)
ggsave("FTF.png", width = 30, height = 4.5, units = "cm")
pl[[2]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(0,25)
ggsave("FTF_n.png", width = 30, height = 4.5, units = "cm")

pl2 <- plot_models(data_path = "FALSE_FALSE_TRUE",
                   title = "FFT",
                   type = "inc", oT = 30)

pl2[[1]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(1.08,1.13)
ggsave("FFT.png", width = 30, height = 4.5, units = "cm")
pl2[[2]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(0,25)
ggsave("FFT_n.png", width = 30, height = 4.5, units = "cm")

pl3 <- plot_models(data_path = "FALSE_TRUE_TRUE",
                   title = "FTT",
                   type = "inc", oT = 30)

pl3[[1]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(1.08,1.13)
ggsave("FTT.png", width = 30, height = 4.5, units = "cm")
pl3[[2]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(0,25)
ggsave("FTT_n.png", width = 30, height = 4.5, units = "cm")

pl4 <- plot_models(data_path = "TRUE_FALSE_FALSE",
                   title = "TFF",
                   type = "high-low", oT = 30)

pl4[[1]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(.85,1.1)
ggsave("TFF.png", width = 30, height = 4.5, units = "cm")
pl4[[2]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(0.3,1.7)
ggsave("TFF_n.png", width = 30, height = 4.5, units = "cm")

pl5 <- plot_models(data_path = "TRUE_FALSE_TRUE",
                   title = "TFT",
                   type = "high-low", oT = 30)

pl5[[1]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(.85,1.1)
ggsave("TFT.png", width = 30, height = 4.5, units = "cm")
pl5[[2]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(0.3,1.7)
ggsave("TFT_n.png", width = 30, height = 4.5, units = "cm")

pl6 <- plot_models(data_path = "TRUE_TRUE_FALSE",
                   title = "TTF",
                   type = "high-low", oT = 30)

pl6[[1]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(.85,1.1)
ggsave("TTF.png", width = 30, height = 4.5, units = "cm")
pl6[[2]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(0.3,1.7)
ggsave("TTF_n.png", width = 30, height = 4.5, units = "cm")


pl7 <- plot_models(data_path = "TRUE_TRUE_TRUE",
                   title = "TTT",
                   type = "high-low", oT = 30)

pl7[[1]] + facet_grid(. ~ model)+ theme(legend.position = "none", text = element_text(size = 8))#+
ylim(.85,1.1)
ggsave("TTT.png", width = 30, height = 4.5, units = "cm")

pl7[[2]] + facet_grid(. ~ model) + theme(legend.position = "none", text = element_text(size = 8))#+
ylim(0.3,1.7)
ggsave("TTT_n.png", width = 30, height = 4.5, units = "cm")


#######################################################
### 1.5 goodness of fit plots (Figure B2)
#######################################################

df <- read.table("~/df_1",
                 header = T)

df$location_id_fac <- factor(df$location_id_fac)

# set parameters
oT <- 30
n_loc = 5
n_per_loc <- 20 # number of locations per group
N0 <- 1000
niter <- 1

ls_diagPlot <- list()

model_func <- c("hgam", "poptrend", "pop_linear", "lin_group")
famil <- c("quasipoisson", "nb")
k <- 1
for(i in model_func) {
  for(j in famil) {
    
    res <- run_mod(get(i), df, j)
    
    ls_diagPlot[[k]] <- appraise(res$trend$gam)[[2]] + ggplot2::labs(title = paste(i, j, sep = "_"))
    
    k <- k+1
  }
  
  
}

pdf("sim_diagPlots_TTF.pdf")

ls_diagPlot[[5]] + ggplot2::ylim(-40,40) + ls_diagPlot[[6]]+ ggplot2::ylim(-40,40) +ls_diagPlot[[7]]+ ggplot2::ylim(-40,40) + ls_diagPlot[[8]] + ggplot2::ylim(-40,40)

ls_diagPlot[[3]]+ ggplot2::ylim(-40,40) + ls_diagPlot[[4]]+ ggplot2::ylim(-40,40) + ls_diagPlot[[1]] + ggplot2::ylim(-40,40) + ls_diagPlot[[2]]+ ggplot2::ylim(-40,40)
dev.off()
