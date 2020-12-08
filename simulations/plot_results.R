library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)

theme_set(theme_bw() +
            theme(plot.background = element_blank(),
                  legend.background = element_blank(),
                  strip.background = element_rect(fill = "white")))

my_palette <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9",
                "#CC79A7")

# prepare to plot
dat <- read_rds("output/sim_study_1-2020-12-05_07_40_52.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  filter(model=="HR") %>% 
  mutate(value=if_else(category=="time", log10(value), value)) %>% 
  #filter(noise=="none") %>% 
  group_by(category, type, n) %>% 
  summarise(mean_value = mean(value))

ggplot(dat, 
       aes(x = n, y = mean_value, col = type)) +
  facet_wrap(vars(category), ncol = 1, scales = "free") +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("value") +
  scale_color_manual(values = rep(my_palette[1:4], 4))




# prepare to plot
dat <- read_rds("output/sim_study_2-2020-11-20_10_01_53.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  mutate(value=if_else(category=="time", log10(value), value)) %>% 
  group_by(category, type, chi.edge) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
  aes(x = chi.edge, y = mean_value, col = type)) +
  facet_wrap(vars(category), ncol = 1, scales = "free") +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("value") +
  scale_color_manual(values = rep(my_palette[1:4], 4))

# dat1 <- read_rds("output/sim_study_3-2020-11-23_09_23_48.rds")
# dat2 <- read_rds("output/sim_study_3-2020-11-23_16_28_27.rds")
# datall <- bind_rows(dat1,dat2)
# saveRDS(datall, file = "output/sim_study_3-HR_nsim300_vario_treeNoise1124.rds")


dat <- read_rds("output/sim_study_3-2020-12-03_08_38_04.rds") %>% 
  unnest(cols = c("perf")) %>% 
  mutate(n=as.character(n)) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  #filter(category=="srr") %>% 
  filter(type=="srr3") %>% 
  filter(noise=="tree") %>% 
  group_by(category, type, p, n) %>% 
  summarise(mean_value = mean(value))


ggplot(dat, 
       aes(x = 1-p, y = mean_value, col = n)) +
  #theme(legend.position = "none") +  f
  facet_wrap(vars(category), ncol = 1, scales = "free") +
  theme(axis.text=element_text(size=22), axis.title=element_text(size=22)) +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("Structure recovery rate error") +
  xlab("Exceedance probability k/n") +
  scale_color_manual(values = rep(my_palette[1:4], 4))+
  ylim(0,1)



# prepare to plot
dat <- read_rds("output/sim_study_3-2020-12-03_14_05_09.rds")  %>% 
  unnest(cols = c("perf")) %>% 
  mutate(category = stringr::str_extract(type, "[a-z]+")) %>% 
  mutate(n=as.character(n)) %>% 
  mutate(value=if_else(category=="time", log10(value), value)) %>% 
  filter(noise=="tree") %>% 
  filter(n=="2000") %>% 
  filter(category=="time") %>% 
  group_by(category, type, p, n) %>% 
  summarise(mean_value = mean(value))

ggplot(dat, 
       aes(x = 1-p, y = mean_value, col = type)) +
  #facet_wrap(vars(type), ncol = 1, scales = "free") +
  geom_line(size = 1, alpha = .75) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  ylab("value") +
  scale_color_manual(values = rep(my_palette[1:4], 4))





