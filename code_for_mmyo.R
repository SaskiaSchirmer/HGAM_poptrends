### load packages

library(tidyverse)
#install.packages("poptrend")
library(poptrend)
library(mgcv)
#install.packages("rtrim")
library(rtrim)
library(patchwork)


######################################################
### 2.1 Estimate population trends for Myotis myotis 
###     with different models - example
######################################################

df <- read.table("mmyo_1991_2019.txt")
df$location_id_fac <- factor(df$location_id)

# use this model when estimating a hierarchical trend:

# HGAM_nb: tempRE = F, family = "nb", type = "smooth"
# HGAM_qp: tempRE = F, family = "quasipoisson", type = "smooth"
# HGLM_nb: tempRE = F, family = "nb", type = "loglinear"
# HGLM_qp: tempRE = F, family = "quasipoisson", type = "loglinear"
# HGAM_nb_tre: tempRE = T, family = "nb", type = "smooth"
# HGAM_qp_tre: tempRE = T, family = "quasipoisson", type = "smooth"
# HGLM_nb_tre: tempRE = T, family = "nb", type = "loglinear"
# HGLM_qp_tre: tempRE = T, family = "quasipoisson", type = "loglinear"

out <- poptrend::ptrend(count ~ trend(year, tempRE = F, type = "smooth") +
                          s(year, location_id_fac, bs = "fs", m=1),
                        family = "nb", data = df, engine = "bam", discrete = TRUE, nthreads = 8)

save(out, file = "mmyo_hgam_nb.Rdata")

# use this mode when estimating a non-hierarchical trend:

# GAM_nb: tempRE = F, family = "nb", type = "smooth"
# GAM_qp: tempRE = F, family = "quasipoisson", type = "smooth"
# GLM_nb: tempRE = F, family = "nb", type = "loglinear"
# GLM_qp: tempRE = F, family = "quasipoisson", type = "loglinear"
# GAM_nb_tre: tempRE = T, family = "nb", type = "smooth"
# GAM_qp_tre: tempRE = T, family = "quasipoisson", type = "smooth"
# GLM_nb_tre: tempRE = T, family = "nb", type = "loglinear"
# GLM_qp_tre: tempRE = T, family = "quasipoisson", type = "loglinear"


out3 <- ptrend(count ~ trend(year, type = 'smooth') + s(location_id_fac, bs = "re"),
               family = "nb",
               data = df, engine = "bam", discrete = TRUE, nthreads = 8)

# TRIM plots
library(rtrim)

out17 <- trim(count~location_id_fac + year, data= df,
              model=2, serialcor=T,overdisp=T,max_iter = 1000)

out18 <- trim(count ~ location_id_fac + year, data = df, model = 3, serialcor=T,overdisp=T,max_iter = 1000)
idx1 <- index(out17)
idx2 <- index(out18)
plot(idx1, idx2, names = c("model 2", "model 3"))



##############################################################
### 2.2 Overviews over the Myotis myotis data (Figure 2)
##############################################################

df %>% 
ggplot(aes(x = year, y = count, color = location_id_fac)) +
  geom_line(linewidth = 1.5) +
  theme(legend.position = "none",
        text = element_text(size = 20)) -> pl1

ggplot(df %>% group_by( year) %>% summarise(me = mean(count, na.rm = TRUE), med = median(count, na.rm = TRUE),
                                            lq = quantile(count, 0.025, na.rm = TRUE), 
                                            uq = quantile(count, 0.975, na.rm = TRUE))) + 
  geom_pointrange(aes(x = year, y = me, ymin = lq, ymax = uq,
                      shape = "mean")) +
  geom_point(aes(x= year, y = med, shape = "median")) +
  scale_shape_manual(values = c(1,2), labels = c("mean", "median")) +
  scale_linetype_manual(values = 1, labels = "true") +
  labs(y = "average population size", color = "model", shape = "Central tendency",
       linetype = "true count", x = "year")  +
  theme(text = element_text(size = 20)) -> pl2


df %>% group_by(year) %>% summarise(mean = mean(count, na.rm = TRUE), 
                                    var = var(count, na.rm = TRUE)) %>% 
  ggplot(aes(x = mean, y = var)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, col = "#E69F00")+
  labs(y = "annual variance of counts over all locations", x = "annual mean count over all locations") +
  theme(text = element_text(size = 20)) -> pl3


pl1 + pl2 + pl3  + plot_annotation(tag_levels= "A")


ggsave(filename = "Figure2.pdf", width = 20, height = 10)






########################################################
### 2.3 plot population trend estimates (Figure 3) - example
########################################################

pdf(file = "mmyo_GLMs.pdf")
par(mfrow = c(2,2))
load("~/mmyo_glm_qp.Rdata")

plot(out8, ylim = c(1,3)) + title("GLM qp")
rm(out8)
load("~/mmyo_glm_nb.Rdata")

plot(out7, ylim = c(1,3)) + title("GLM nb")
rm(out7)

load("~/mmyo_hglm_qp.Rdata")
plot(out6, ylim = c(1,3)) + title("HGLM qp")
rm(out6)

load("~/mmyo_hglm_nb.Rdata")

plot(out5, ylim = c(1,3)) + title("HGLM nb")
rm(out5)

dev.off()



###################################################################
### 2.4 plot diagnostics (Figure 4) - example
###################################################################

tmp <- appraise(out6$gam)
pl6 <- tmp[[2]] + ylim(-20,20) + labs(title = "HGLM qp")
rm(out6)


load("~/mmyo_hglm_nb.Rdata")
tmp <- appraise(out5$gam)
pl5 <- tmp[[2]] + ylim(-20,20) + labs(title = "HGLM nb")
rm(out5)
load("~/mmyo_glm_qp.Rdata")
tmp <- appraise(out8$gam)
pl8 <- tmp[[2]] + ylim(-20,20) + labs(title = "GLM qp")
rm(out8)

load("~/mmyo_glm_nb.Rdata")
tmp <- appraise(out7$gam)
pl7 <- tmp[[2]] + ylim(-20,20) + labs(title = "GLM nb")
rm(out7)


pdf(file = "mmyo_GLMs_diagPlots.pdf")
pl8 + pl7 + pl6 + pl5
dev.off()
