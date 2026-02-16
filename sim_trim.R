library(tidyverse)
#install.packages("poptrend")
library(poptrend)
library(mgcv)

#install.packages("Downloads/rtrim_2.1.1.tar.gz", repos = NULL, type ="source")
library(rtrim)

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

## compare models

setting <- expand.grid(nonlinear = c(FALSE, TRUE),
                       group = c(FALSE, TRUE),
                       negativebinomial = c(FALSE, TRUE))


trim <- function(df, cp, sc = TRUE) {
    change <- 1
    if(cp) change <- "auto"
    rtrim::trim(count ~ location_id_fac + year, data = df, model = 2, overdisp = TRUE,
                changepoints=change, serialcor = sc)
}




run_mod <- function(model_func = trim, df, cp) {
    warns <- list()
    
    trend <- NULL
    
    
    t <- tryCatch({trend <- model_func(df, cp)
    idx <- rtrim::index(trend, level = 0.95, which = "fitted")
    idx$iter <- i
    idx
    }, 
    warning = function(warn) {
        if(grepl("Serial correlation(.*?)consider disabling it.", warn$message)){
            print("Serial correlation disabled")
            trend <- model_func(df, cp, FALSE)
            idx <- rtrim::index(trend, level = 0.95, which = "fitted")
            idx$iter <- i
            idx
        }
        
    })
    t
}



run_one_mod <- function(data, oT, n_loc, n_per_loc, N0, 
                        cp, 
                        size, model_func, i) {
    
    
    
    trend_hg_nb <- run_mod(model_func, data, cp)
    
    
    
    trend_hg_nb
}


n.cores <- 6

model_name <- "trim_cp"
model_function <- trim
cp <- TRUE


# set parameters
oT <- 30
n_loc = 5
n_per_loc <- 20 # number of locations per group
N0 <- 1000
start_i <- 1
niter <- 250

# for every data type in 1:8 do
#for(i in 1:nrow(setting)){
i <- 3
my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
set <- setting[i,]

path <- paste("simulated_data/", paste(setting[i,], collapse = "_"), sep = "")
ls_files <- list.files(path, full.names = TRUE)
# init run


x <- foreach(
    i = start_i:niter, 
    .combine = 'rbind'
) %dopar% {
    library(tidyverse)
    df <- read.table(list.files(path, full.names = TRUE)[i], header = TRUE)
    df$location_id_fac <- factor(df$location_id_fac)
    run_one_mod(df, oT, n_loc, n_per_loc, N0, cp, 100,
                model_function, i)
}
x

dir_path <- paste("sim_results/",paste(setting[i,], collapse = "_"), sep = "")
if(!dir.exists(dir_path)) {
    dir.create(dir_path)
}

save(x, file = paste(dir_path,
                     "/", model_name, "_", cp,"_", start_i, "_", niter, ".Rdata", sep = ""))

est_year_gs_qu <- x %>% group_by(iter) %>% 
    mutate(lambda = exp((log(fitted)-lag(log(fitted), order_by = time))/(time - lag(time, order_by = time))),
           N = fitted) %>% 
    group_by(time) %>% filter(!is.na(lambda)) %>% 
    rename(year = time) %>%
    summarise(me_lambda = mean(lambda), 
              lq_lambda = quantile(lambda, prob = 0.025),
              uq_lambda = quantile(lambda, prob = 0.975),
              sd_lambda = sd(lambda),
              me_N = mean(N), 
              lq_N = quantile(N, prob = 0.025),
              uq_N = quantile(N, prob = 0.975),
              sd_N = sd(N)) 

save(est_year_gs_qu, file = paste(dir_path,
                                  "/", model_name, "_", cp,"summary.Rdata", sep = ""))

parallel::stopCluster(cl = my.cluster)
#}
