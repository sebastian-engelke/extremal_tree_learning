library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(latex2exp)
library(graphicalExtremes)
theme_set(theme_bw() +
            theme(plot.background = element_blank(),
                  legend.background = element_blank(),
                  strip.background = element_rect(fill = "white")))

my_palette <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9",
                "#CC79A7")

## Plot variogram against chi 

dirichletTheta = function(a1, a2){
  n = 50000
  la1 = length(a1)
  la2 = length(a2)
  stopifnot(la1 == la2)
  if(la1 == 1)
    mean(pmax(rgamma(n, shape=a1, rate=a1), rgamma(n, shape=a2, rate=a2)))
  else
    rowMeans(pmax(matrix(rgamma(n*la1, shape=matrix(a1, nrow=la1, ncol=n), rate=matrix(a1, nrow=la1, ncol=n)), nrow=la1, ncol=n), 
                  matrix(rgamma(n*la2, shape=matrix(a2, nrow=la2, ncol=n), rate=matrix(a2, nrow=la2, ncol=n)), nrow=la2, ncol=n)))
} 

dirichletThetaSym = function(a) dirichletTheta(a,a)

chi.vec <-  seq(0.1,0.9, len=25)
log.par.vec <-  log(2-chi.vec)/log(2)
neglog.par.vec <-  -log(2)/log(chi.vec)
set.seed(331256)
dirichlet.vec <-  sapply(chi.vec, FUN=function(chi) uniroot(function(a) dirichletThetaSym(a) - (2-chi), interval = c(0.01,200))$root)


Gamma.vec = c(chi2Gamma(chi.vec),0)
G.dirichlet = c(trigamma(dirichlet.vec) + trigamma(dirichlet.vec + 1),0)
G.logistic = c(log.par.vec^2 * (trigamma(1-log.par.vec) + pi^2/6),0)

dat <- tibble(chi.vec=c(chi.vec,1), G.log = G.logistic, G.HR = Gamma.vec, G.dirichlet = G.dirichlet) %>% 
  gather(key = "variable", value = "value", -chi.vec)

ggplot(dat, 
       aes(x = chi.vec, y = value, col = variable)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  xlab(expression(paste("Extremal correlation ", chi, sep="")))+
  ylab(expression(paste("Extremal variogram ", Gamma, sep=""))) +
  coord_cartesian(xlim=c(0,1),ylim = c(0,10))+
  scale_color_manual(values = rep(my_palette[1:3]))
ggsave(filename = "Gamma_Chi.pdf", height = 7, width=7)



###################################
##### Simulation Study 1: IID Noise
###################################


##### Simulation Study 1: IID Noise, Huesler--Reiss

dat <- read_rds("data/sim_study_1-HR_dirichlet_nsim300.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="HR") %>% 
  filter(category=="err") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))

ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Wrong edge rate") +
  scale_color_manual(values = rep(my_palette[1:4], 3)) +
  ylim(0,.5)
ggsave(filename = "HR_err_d20.pdf", height = 7, width=7)



dat <- read_rds("data/sim_study_1-HR_dirichlet_nsim300.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="HR") %>% 
  filter(category=="srr") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  scale_color_manual(values = rep(my_palette[1:4], 3)) +
  ylim(0,1)
ggsave(filename = "HR_srr_d20.pdf", height = 7, width=7)



##### Simulation Study 1: IID Noise, Dirichlet


dat <- read_rds("data/sim_study_1-HR_dirichlet_nsim300.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="dirichlet") %>% 
  filter(category=="err") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))

ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Wrong edge rate") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,.5)
ggsave(filename = "dirichlet_err_d20.pdf", height = 7, width=7)



dat <- read_rds("data/sim_study_1-HR_dirichlet_nsim300.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="dirichlet") %>% 
  filter(category=="srr") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "dirichlet_srr_d20.pdf", height = 7, width=7)



###################################
##### Simulation Study 1: Tree Noise
###################################



##### Simulation Study 1: Tree Noise, Huesler--Reiss



dat <- read_rds("data/sim_study_1-HR_dirichlet_nsim300_treeNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="HR") %>% 
  filter(category=="err") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))

ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Wrong edge rate") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,.5)
ggsave(filename = "HR_err_d20_treeNoise.pdf", height = 7, width=7)

dat <- read_rds("data/sim_study_1-HR_dirichlet_nsim300_treeNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="HR") %>% 
  filter(category=="srr") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))

ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "HR_srr_d20_treeNoise.pdf", height = 7, width=7)




##### Simulation Study 1: Tree Noise, Dirichlet

dat <- read_rds("data/sim_study_1-HR_dirichlet_nsim300_treeNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="dirichlet") %>% 
  filter(category=="err") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))

ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Wrong edge rate") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,.5)
ggsave(filename = "dirichlet_err_d20_treeNoise.pdf", height = 7, width=7)


dat <- read_rds("data/sim_study_1-HR_dirichlet_nsim300_treeNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="dirichlet") %>% 
  filter(category=="srr") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))

ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "dirichlet_srr_d20_treeNoise.pdf", height = 7, width=7)



###################################
##### Simulation Study 1: Time
###################################


dat <- read_rds("data/sim_study_1-HR_dirichlet_nsim300_treeNoise_time.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="HR") %>% 
  filter(category=="time") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))

ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) + 
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  #facet_wrap(vars(category), ncol = 1, scales = "free") +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab(TeX("Computaiton time in seconds ($log_{10}$ scale)")) +
  scale_y_log10(breaks=c(0.01, 0.1, 1, 10, 100), labels=c("0.01", "0.1", "1", "10", "100")) +
  scale_color_manual(values = rep(my_palette[1:4], 3))
ggsave(filename = "HR_time_d20_treeNoise.pdf", height = 7, width=7)



###################################
##### Simulation Study 2: IID Noise
###################################

dat <- read_rds("data/sim_study_2-HR_nsim300.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(category=="srr") %>% 
  group_by(category, type, chi.edge) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = chi.edge, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  xlab(expression(paste("Extremal correlation ", chi, sep=""))) +
  coord_cartesian(xlim=c(0,1))+
  scale_color_manual(values = rep(my_palette[1:4], 3))
ggsave(filename = "HR_srr_d20_chiedge.pdf", height = 7, width=7)



###################################
##### Simulation Study 2: Tree Noise
###################################

dat <- read_rds("data/sim_study_2-HR_nsim300_treeNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  mutate(value=if_else(category=="time", log10(value), value)) %>% 
  filter(category=="srr") %>% 
  group_by(category, type, chi.edge) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = chi.edge, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  xlab(expression(paste("Extremal correlation ", chi, sep=""))) +
  coord_cartesian(xlim=c(0,1))+
  scale_color_manual(values = rep(my_palette[1:4], 3))
ggsave(filename = "HR_srr_d20_chiedge_treeNoise.pdf", height = 7, width=7)




###################################
##### Simulation Study 3: IID Noise
###################################

##### Simulation Study 3: IID Noise, n=500

dat <- read_rds("data/sim_study_3-HR_nsim300_vario_CL_bothNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(n=as.character(n)) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(category=="srr") %>% 
  filter(noise=="iid") %>% 
  filter(n=="500") %>% 
  group_by(category, type, p, n) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = 1-p, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  xlab("Exceedance probability k/n") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "HR_srr_n500_d20_kn.pdf", height = 7, width=7)



##### Simulation Study 3: IID Noise, n=1000

dat <- read_rds("data/sim_study_3-HR_nsim300_vario_CL_bothNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(n=as.character(n)) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(category=="srr") %>% 
  filter(noise=="iid") %>% 
  filter(n=="1000") %>% 
  group_by(category, type, p, n) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = 1-p, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("") +
  xlab("Exceedance probability k/n") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "HR_srr_n1000_d20_kn.pdf", height = 7, width=7)





##### Simulation Study 3: IID Noise, n=2000

dat <- read_rds("data/sim_study_3-HR_nsim300_vario_CL_bothNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(n=as.character(n)) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(category=="srr") %>% 
  filter(noise=="iid") %>% 
  filter(n=="2000") %>% 
  group_by(category, type, p, n) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = 1-p, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("") +
  xlab("Exceedance probability k/n") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "HR_srr_n2000_d20_kn.pdf", height = 7, width=7)

###################################
##### Simulation Study 3: Tree Noise
###################################

##### Simulation Study 3: Tree Noise, n=500

dat <- read_rds("data/sim_study_3-HR_nsim300_vario_CL_bothNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(n=as.character(n)) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(category=="srr") %>% 
  filter(noise=="tree") %>% 
  filter(n=="500") %>% 
  group_by(category, type, p, n) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = 1-p, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  xlab("Exceedance probability k/n") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "HR_srr_n500_d20_kn_treeNoise.pdf", height = 7, width=7)



##### Simulation Study 3: Tree Noise, n=1000

dat <- read_rds("data/sim_study_3-HR_nsim300_vario_CL_bothNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(n=as.character(n)) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(category=="srr") %>% 
  filter(noise=="tree") %>% 
  filter(n=="1000") %>% 
  group_by(category, type, p, n) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = 1-p, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("") +
  xlab("Exceedance probability k/n") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "HR_srr_n1000_d20_kn_treeNoise.pdf", height = 7, width=7)




##### Simulation Study 3: Tree Noise, n=2000

dat <- read_rds("data/sim_study_3-HR_nsim300_vario_CL_bothNoise.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(n=as.character(n)) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(category=="srr") %>% 
  filter(noise=="tree") %>% 
  filter(n=="2000") %>% 
  group_by(category, type, p, n) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = 1-p, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("") +
  xlab("Exceedance probability k/n") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "HR_srr_n2000_d20_kn_treeNoise.pdf", height = 7, width=7)





###################################
##### Simulation Study 4: IID Noise
###################################

##### Simulation Study 4: IID Noise, n=1000

dat <- read_rds("data/sim_study_4-HR_nsim300_d10_20_30_50_100_ML_200_300.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="HR") %>% 
  filter(category=="srr") %>% 
  mutate(value=if_else(category=="time", log10(value), value)) %>% 
  group_by(category, type, d) %>% 
  summarise(mean_value = mean(value)) %>% 
  mutate(mean_value = if_else(d > 100 & type == "srr4", NaN, mean_value))



ggplot(dat, 
       aes(x = d, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  xlab("Dimension d") +
  scale_color_manual(values = rep(my_palette[1:4], 3))+
  ylim(0,1)
ggsave(filename = "HR_srr_n1000_dim.pdf", height = 7, width=7)




dat <- read_rds("data/sim_study_4-HR_nsim300_d10_20_30_50_100_ML_200_300.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="HR") %>% 
  filter(category=="time") %>% 
  mutate(value=if_else(category=="time", log10(value), value)) %>% 
  group_by(category, type, d) %>% 
  summarise(mean_value = mean(value)) %>% 
  mutate(mean_value = if_else(d > 100 & type == "time4", NaN, mean_value))




ggplot(dat, 
       aes(x = d, y = mean_value, col = type)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  xlab("Dimension d") +
  ylab(TeX("Computaiton time in seconds ($log_{10}$ scale)")) +
  scale_color_manual(values = rep(my_palette[1:4], 3))
ggsave(filename = "HR_srr_n1000_dim_time.pdf", height = 7, width=7)




