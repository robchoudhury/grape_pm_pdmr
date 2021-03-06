library(tidyverse)
library(lubridate)
library(agricolae)
library(lme4)
library(multcomp)
library(emmeans)

cultivar=read_csv("data/cultivar_codes.csv")
data=read_csv("data/all_shasta_data.csv") %>%
  rename_all(. %>%tolower) %>%
  rename(leaf_severity=`leaf severity (%)`,
         leaf_incidence=`leaf incidence`,
         cluster_incidence=`cluster inc`,
         cluster_severity = `cluster sev`) %>%
  mutate(date=lubridate::mdy(date),
         cluster_incidence=as.numeric(cluster_incidence),
         cluster_severity=as.numeric(cluster_severity)*100,
         year=lubridate::year(date),
         yday=lubridate::yday(date)) %>%
  full_join(., cultivar) %>%
  mutate(cultivar_name=tolower(cultivar_name)) %>%
  filter(rater == "rac")

data_summary=data  %>%
  group_by(year, cultivar_name, yday) %>%
  summarize(mean_leaf_inc=mean(leaf_incidence, na.rm = T),
            mean_leaf_sev=mean(leaf_severity, na.rm = T),
            mean_cluster_inc=mean(cluster_incidence, na.rm = T),
            mean_cluster_sev=mean(cluster_severity, na.rm = T))

data_summary_rep=data  %>%
  group_by(year, cultivar_name, yday, rep) %>%
  summarize(mean_leaf_inc=mean(leaf_incidence, na.rm = T),
            mean_leaf_sev=mean(leaf_severity, na.rm = T),
            mean_cluster_inc=mean(cluster_incidence, na.rm = T),
            mean_cluster_sev=mean(cluster_severity, na.rm = T))

audpc_summary=data %>%
  filter(rater == "rac") %>%
  group_by(rep, year, cultivar_name, yday) %>%
  summarize(mean_leaf_inc=mean(leaf_incidence, na.rm = T),
            mean_leaf_sev=mean(leaf_severity, na.rm = T),
            mean_cluster_inc=mean(cluster_incidence, na.rm = T),
            mean_cluster_sev=mean(cluster_severity, na.rm = T)) %>%
  group_by(cultivar_name, year, rep) %>%
  drop_na() %>%
    summarize(audpc_leaf=audpc(evaluation = mean_leaf_sev, dates = yday ),
              audpc_cluster=audpc(evaluation = mean_cluster_sev, dates = yday ))

#2012
data_summary_2012=data %>%
  filter(rater == "rac" & year==2012 & yday==265) %>%
  group_by(cultivar_name, yday) %>%
  summarize(mean_leaf_inc=mean(leaf_incidence, na.rm = T),
            mean_leaf_sev=mean(leaf_severity, na.rm = T),
            mean_cluster_inc=mean(cluster_incidence, na.rm = T),
            mean_cluster_sev=mean(cluster_severity, na.rm = T))

audpc_summary_2012_leaf=data %>%
  filter(rater == "rac"& year==2012) %>%
  group_by(rep, cultivar_name, yday) %>%
  summarize(mean_leaf_inc=mean(leaf_incidence, na.rm = T),
            mean_leaf_sev=mean(leaf_severity, na.rm = T),
            mean_cluster_inc=mean(cluster_incidence, na.rm = T),
            mean_cluster_sev=mean(cluster_severity, na.rm = T)) %>%
  group_by(cultivar_name,  rep) %>%
  summarize(audpc_leaf=audpc(evaluation = mean_leaf_sev, dates = yday ))

audpc_summary_2012_cluster=data %>%
  filter(rater == "rac"& year==2012) %>%
  group_by(rep, cultivar_name, yday) %>%
  summarize(mean_leaf_inc=mean(leaf_incidence, na.rm = T),
            mean_leaf_sev=mean(leaf_severity, na.rm = T),
            mean_cluster_inc=mean(cluster_incidence, na.rm = T),
            mean_cluster_sev=mean(cluster_severity, na.rm = T)) %>%
  group_by(cultivar_name,  rep) %>%
  drop_na() %>%
  summarize(audpc_cluster=audpc(evaluation = mean_cluster_sev, dates = yday ))

ggplot(data_summary_2012, 
       aes(reorder(cultivar_name, mean_leaf_sev, FUN = mean, na.rm=T), 
           mean_leaf_sev, color=cultivar_name))+
  geom_jitter()+
  geom_boxplot(alpha=0.3)
  
data_anova_2012=data %>%
  filter(rater == "rac" & year==2012 & yday==265)

grape_aov <- aov(cluster_severity ~ rep+cultivar, data=data_anova_2012)
shapiro.test(grape_aov$resid)
qqnorm(grape_aov$resid)
qqline(grape_aov$resid)
bartlett.test(data_anova_2012$cluster_severity, data_anova_2012$cultivar)
fligner.test(data_anova_2012$cluster_severity, data_anova_2012$cultivar)

grape_aov <- aov(leaf_severity ~ rep+cultivar, data=data_anova_2012)
shapiro.test(grape_aov$resid)
qqnorm(grape_aov$resid)
qqline(grape_aov$resid)
bartlett.test(data_anova_2012$leaf_severity, data_anova_2012$cultivar)
fligner.test(data_anova_2012$leaf_severity, data_anova_2012$cultivar)
grape_aov_block <- aov(cluster_severity ~ rep, data=data_anova_2012)
summary(grape_aov_block)
grape_aov_block <- aov(leaf_severity ~ rep, data=data_anova_2012)
summary(grape_aov_block)
block_effect= data_anova_2012 %>%
  group_by(rep) %>%
  summarize(mean_leaf=mean(leaf_severity),
            se_leaf=sd(leaf_severity),
            mean_cluster=mean(cluster_severity),
            se_cluster=sd(cluster_severity))



grape_lme_aov <- lmer(leaf_severity ~ cultivar_name+(1|rep), data=data_anova_2012, REML=F)
anova(grape_lme_aov)
grape_lme_null=lmer(leaf_severity ~ (1|rep), data=data_anova_2012, REML=F) # reduced model
anova(grape_lme_null, grape_lme_aov) # anova function to compare the full and reduced models
summary(grape_lme_aov)
ranef(grape_lme_aov)
grape_glht<-summary(glht(grape_lme_aov, linfct=mcp(cultivar_name="Tukey")))
grape_glht # print all pairwise comparisons
cld(grape_glht)
par(mar=c(5,22,3,1), mgp=c(2,.7,0)) # set wide left margin 
plot(grape_glht)
lsmeans(grape_lme_aov, pairwise ~ cultivar_name, data=data_anova_2012, adjust="none")

grape_lme_aov <- lmer(cluster_severity ~ cultivar_name+(1|rep), data=data_anova_2012, REML=F)
anova(grape_lme_aov)
grape_lme_null=lmer(cluster_severity ~ (1|rep), data=data_anova_2012, REML=F) # reduced model
anova(grape_lme_null, grape_lme_aov) # anova function to compare the full and reduced models
summary(grape_lme_aov)
ranef(grape_lme_aov)
grape_glht_cluster<-summary(glht(grape_lme_aov, linfct=mcp(cultivar_name="Tukey")))
grape_glht_cluster # print all pairwise comparisons
cld(grape_glht_cluster)
par(mar=c(5,22,3,1), mgp=c(2,.7,0)) # set wide left margin 
plot(grape_glht_cluster)
lsmeans(grape_lme_aov, pairwise ~ cultivar_name, data=data_anova_2012, adjust="none")

ggplot(data_summary, aes(mean_leaf_sev, mean_leaf_inc, color=cultivar_name))+
  geom_point()

ggplot(data_summary, aes(mean_cluster_sev, mean_cluster_inc))+
  geom_point(aes(color=cultivar_name))+
  geom_smooth()

test=aov(cluster_severity ~ rep, data=data_anova_2012)
summary(test)
test=aov(leaf_severity ~ rep, data=data_anova_2012)
summary(test)

cld(grape_glht)
cld(grape_glht_cluster)

data = data %>%
  mutate(yday = yday(date))

ggplot(data, aes(yday, 
                 leaf_severity, 
                 group = cultivar_name, 
                 color = cultivar_name)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", na.rm = T, size =2)+
  stat_summary(fun.y = mean, geom = "line", na.rm = T, size =2)+
  facet_grid(year~., scales = "free") +
  scale_color_viridis_d(option = "C", 
                        guide = guide_legend(ncol = 3,title.position = "top")) +
  theme_dark() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  ylab("Leaf Disease Severity (%)") +
  xlab("Julian Date")
ggsave(filename = "figures/leaf_severity.png", width = 8, height = 6, units = "in")  

ggplot(data, aes(yday, 
                 cluster_severity, 
                 group = cultivar_name, 
                 color = cultivar_name)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", na.rm = T, size =2)+
  stat_summary(fun.y = mean, geom = "line", na.rm = T, size =2)+
  facet_grid(year~., scales = "free") +
  scale_color_viridis_d(option = "C", 
                        guide = guide_legend(ncol = 3,title.position = "top")) +
  theme_dark() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  ylab("Leaf Disease Severity (%)") +
  xlab("Julian Date")
ggsave(filename = "figures/cluster_severity.png", width = 8, height = 6, units = "in")  

#####
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}
#####
ggplot(audpc_summary, aes(reorder(cultivar_name, audpc_leaf, FUN = mean), 
                 audpc_leaf, 
                 group = cultivar_name, 
                 color = cultivar_name)) +
  geom_jitter(size = 2)+
  stat_summary(fun.data = mean_se, geom = "errorbar", na.rm = T, size =2)+
  facet_grid(year~., scales = "free") +
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)})+
  scale_color_viridis_d(option = "C", 
                        guide = guide_legend(ncol = 3,title.position = "top")) +
  theme_dark() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  ylab("Area Under Disease Progress Curve (AUDPC") +
  xlab("Cultivar Name")
ggsave(filename = "figures/audpc_leaf.png", width = 10, height = 6, units = "in")  


ggplot(data_summary_rep, aes(mean_leaf_sev, 
                             mean_cluster_sev, 
                             color = yday))+
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis_c(option = "D", 
                        guide = guide_legend(title.position = "top")) +
  theme_dark() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  ylab("Mean Cluster Severity (%)") +
  xlab("Mean Leaf Severity (%)")
ggsave(filename = "figures/leaf_cluster_concordance.png", width = 6, height = 6, units = "in")  

###
ggplot(audpc_summary, aes(audpc_leaf, 
                          audpc_cluster, 
                             color = cultivar_name))+
  geom_point(size = 2) +
  #geom_abline(slope = 1, intercept = 0) +
  scale_color_viridis_d(option = "D", 
                        guide = guide_legend(ncol = 3, title.position = "top")) +
  theme_dark() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  ylab("Cluster AUDPC") +
  xlab("Leaf AUDPC")
ggsave(filename = "figures/leaf_cluster_audpc_concordance.png", width = 6, height = 6, units = "in")  
