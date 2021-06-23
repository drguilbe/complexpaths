#Topological Measures for Identifying and Predicting the Spread of Complex Contagions
#Guilbeault & Centola, 2021, Nature Communications 
#June 2021 

rm(list=ls());gc(); options(max.print=999999); options(warn=-1)
library(dplyr);
library(ggplot2);
library(RVAideMemoire);
library(sjPlot)
library(multiwayvcov)
library(lmtest)
library(tidyr)

min_max_norm<-function(x){(x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x,na.rm=TRUE))}

#load Data
fig2_main_SI_data<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/fig2_main_SI_data.csv")
fig2ABCD_node_data<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/fig2ABCD_node_data.csv")
fig2_cascade_data<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/fig2_cascade_data.csv")
fig3<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/fig3_revised.csv")
fig4<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/fig4.csv")
figs7a_node_data<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs7a_node_data.csv")
figs7a_cascade<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs7a_cascade_data.csv")
figs7b_node_data<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs7b_node_data.csv")
figs8_s9_main_SI_data<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs8_s9_main_SI_data.csv")
figs11_to_s16_SI_data<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs11_to_s16_SI_data.csv")
fig_s13_raw<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/fig_s13_raw.csv")
figs16ab<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs16ab.csv")
figs17<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs17.csv")
figs18a_main<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs18a_main.csv")
figs18a_perc<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs18a.csv")
figs18b<-read.csv("C:/Users/dougl/Desktop/Complex_Paths_Data/figs18b.csv")

#######################
###ORGANIZE RAW DATA###
#######################

##############
#Main Figures#
##############

fig2_main_SI_data_homo<-subset(fig2_main_SI_data, grepl("homo", cond)) 
fig2_main_SI_data_hetero<-subset(fig2_main_SI_data, grepl("hetero", cond)) 

fig2_main_SI_data_homo_agg<-fig2_main_SI_data_homo %>% group_by(T, p, simmod, figure, cond) %>% 
  summarise(avg_bw_size = mean(avg_bw_size, na.rm=T),simp_deg = mean(simp_deg, na.rm=T), 
            step_size = mean(step_size, na.rm=T),num_adopt = mean(num_adopt, na.rm=T), plc = mean(bw_len, na.rm=T)) %>%
  group_by(T, figure, cond) %>% mutate(num_adopt_norm=min_max_norm(num_adopt), step_size_norm=min_max_norm(step_size), plc_norm=min_max_norm(plc))
fig2_main_SI_data_hetero_agg<-fig2_main_SI_data_hetero %>% group_by(p, simmod, figure, cond) %>% 
  summarise(avg_bw_size = mean(avg_bw_size, na.rm=T), simp_deg = mean(simp_deg, na.rm=T),
            step_size = mean(step_size, na.rm=T), num_adopt = mean(num_adopt, na.rm=T), plc = mean(bw_len, na.rm=T)) %>%
  group_by(figure, cond) %>% mutate(num_adopt_norm=min_max_norm(num_adopt), step_size_norm=min_max_norm(step_size),  plc_norm=min_max_norm(plc));
fig2_main_SI_data_hetero_agg$T<-NA
fig2_main_SI_data_agg<-rbind(fig2_main_SI_data_homo_agg,fig2_main_SI_data_hetero_agg);

fig2ABC<-subset(fig2ABCD_node_data, condition == "HomoThresh")

fig2ABC_agg1<-fig2ABC %>% group_by(N, k, T, p) %>% summarise(prop_steps_over_Wc=sum(avg_step_over_Wc)/length(avg_step_over_Wc) )
fig2ABC_agg2<-fig2ABC %>% group_by(N, k, T, p) %>% summarise_all(funs(mean(.,na.rm=T)))
fig2ABC_agg<-merge(fig2ABC_agg1, fig2ABC_agg2, by=c("N", "k", "T", "p"))
fig2ABC_final<-merge(subset(fig2_cascade_data, condition == "HomoThresh"), fig2ABC_agg, by=c("N", "k", "T", "p", "sim", "ID"))

fig2D<-subset(fig2ABCD_node_data, condition != "HomoThresh")
fig2D_agg1<-fig2D %>% group_by(N, k, p, simmod) %>% summarise(prop_steps_over_Wc=sum(avg_step_over_Wc)/length(avg_step_over_Wc) )
fig2D_agg_final<-fig2D_agg1 %>% group_by(N, k, p) %>% summarise_all(funs(mean(.,na.rm=T)))
fig2D_agg_final$prop_steps_over_Wc_norm<- min_max_norm(fig2D_agg_final$prop_steps_over_Wc)
fig2D_cascade_final<-subset(fig2_cascade_data, condition == "HeteroThresh")
fig2D_cascade_final<-fig2D_cascade_final %>% group_by(p,k,N) %>% summarise(prop_cascade = mean(prop_cascade))
fig2D_final<-merge(fig2D_cascade_final, fig2D_agg_final, by=c("N", "k", "p"))

fig3_top_perc<-fig3 %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=ci_norm, 1); fig3_top_perc$seeding.strat.org<-"Percolation"
fig3_top_core<-fig3 %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=core, 1); fig3_top_core$seeding.strat.org<-"K-core"
fig3_top_btwn<-fig3 %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=simp_btwn, 1); fig3_top_btwn$seeding.strat.org<-"Betweenness"
fig3_top_eigen<-fig3 %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=simp_eigen, 1); fig3_top_eigen$seeding.strat.org<-"Eigen."
fig3_top_deg<-fig3 %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=simp_deg, 1); fig3_top_deg$seeding.strat.org<-"Degree"
fig3_top_plcc<-fig3 %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=bw_len, 1); fig3_top_plcc$seeding.strat.org<-"Complex"
fig3_all<-rbind(fig3_top_perc, fig3_top_core, fig3_top_btwn, fig3_top_eigen, fig3_top_deg, fig3_top_plcc); fig3_all$FinalAdopt_mod<-fig3_all$num_adopt-fig3_all$num_seeds
fig3_plot<-fig3_all %>% group_by(netname, T, N) %>% mutate(PLc_norm = min_max_norm(bw_len), adopt_norm=min_max_norm(FinalAdopt_mod))
fig3_plot<-fig3_plot %>% group_by(netname, T) %>% mutate(Degree_norm = min_max_norm(simp_deg), core_norm = min_max_norm(core)); 
fig3_plot$seeding.strat.org<-as.factor(fig3_plot$seeding.strat.org); levels(fig3_plot$seeding.strat.org)<-c("Top.Btwn","Top.Complex", "Top.Deg","Top.Eigen","Top.Core", "Top.Perc.")
fig3_plot$Btwn.cent<-fig3_plot$simp_btwn; fig3_plot$Eigen.cent<-fig3_plot$simp_eigen
fig3_plot$seeding.strat.org <- factor(fig3_plot$seeding.strat.org, levels = c("Top.Eigen", "Top.Deg","Top.Btwn", "Top.Core", "Top.Perc.", "Top.Complex"))
fig3_plot_final<-fig3_plot %>% group_by(netname, seeding.strat.org) %>% summarise_all(mean, na.rm=T)

fig4A<-subset(fig4, takeup == 1 & hasleader == 1)
fig4A_avg<- fig4A %>% group_by(hhid, vid, netname, takeup, takeup_count_i, takeup_count, takeup_count_nonleader_i, 
                               takeup_pct, takeup_pct_nonleader, nrooms, nbeds, electricity, latrine, nroomscapita, nbedscapita, trimester, latrine_backup) %>% 
  dplyr::summarise(simp_btwn = mean(simp_btwn), percolation = mean(percolation), simp_close = mean(simp_close),
                   simp_eigen = mean(simp_eigen), simp_deg = mean(simp_deg), simp_reach = mean(simp_reach), core = mean(core),
                   bw_len = mean(bw_len), diffusion_cent = mean(diffusion_b), ci_norm=mean(ci_norm) )

fig4A_plc<-fig4A_avg %>% group_by(vid, netname) %>% top_n(wt=bw_len, 1) %>% mutate(cond="complex")
fig4A_btwn<-fig4A_avg %>% group_by(vid, netname) %>% top_n(wt=simp_btwn, 1) %>% mutate(cond="btwn")
fig4A_close<-fig4A_avg %>% group_by(vid, netname) %>% top_n(wt=simp_close, 1) %>% mutate(cond="close")
fig4A_deg<-fig4A_avg %>% group_by(vid, netname) %>% top_n(wt=simp_deg, 1) %>% mutate(cond="deg")
fig4A_percolation<-fig4A_avg %>% group_by(vid, netname) %>% top_n(wt=ci_norm, 1) %>% mutate(cond="percolation")
fig4A_eigen<-fig4A_avg %>% group_by(vid, netname) %>% top_n(wt=simp_eigen, 1) %>% mutate(cond="eigen")
fig4A_core<-fig4A_avg %>% group_by(vid, netname) %>% top_n(wt=core, 1) %>% mutate(cond="core")
fig4A_diff<-fig4A_avg %>% group_by(vid, netname) %>% top_n(wt=diffusion_cent, 1) %>% mutate(cond="diffusion_cent")
fig4A_seeds<-rbind(fig4A_btwn, fig4A_close, fig4A_deg, fig4A_eigen, fig4A_diff, fig4A_percolation)
fig4A_random<-fig4A_seeds %>% group_by(vid, netname) %>% sample_n(1) %>% mutate(cond="Random")
fig4A_all<-rbind(fig4A_plc, fig4A_btwn, fig4A_close, fig4A_deg, fig4A_percolation, fig4A_eigen, fig4A_core, fig4A_diff, fig4A_random)
fig4A_all<-subset(fig4A_all, !cond %in% c("diffusion_cent", "close"))
fig4A_all$cond<-as.factor(fig4A_all$cond)
levels(fig4A_all$cond)<-c("Betweenness", "Complex", "K-core", "Degree", "Eigen", "Percolation", "Random")

fig4B<-subset(fig4, takeup == 1)
fig4B_avg<- fig4B %>% group_by(hhid, vid, netname, takeup, takeup_count_i, takeup_count, takeup_count_nonleader_i, 
                               takeup_pct, takeup_pct_nonleader, nrooms, nbeds, electricity, latrine, nroomscapita, nbedscapita, trimester, latrine_backup) %>% 
  dplyr::summarise(simp_btwn = mean(simp_btwn), percolation = mean(percolation), simp_close = mean(simp_close),
                   simp_eigen = mean(simp_eigen), simp_deg = mean(simp_deg), simp_reach = mean(simp_reach), core = mean(core),
                   bw_len = mean(bw_len), diffusion_cent = mean(diffusion_b), ci_norm=mean(ci_norm) )

fig4B_plc<-fig4B_avg %>% group_by(vid, netname) %>% top_n(wt=bw_len, 1) %>% mutate(cond="complex")
fig4B_btwn<-fig4B_avg %>% group_by(vid, netname) %>% top_n(wt=simp_btwn, 1) %>% mutate(cond="btwn")
fig4B_close<-fig4B_avg %>% group_by(vid, netname) %>% top_n(wt=simp_close, 1) %>% mutate(cond="close")
fig4B_deg<-fig4B_avg %>% group_by(vid, netname) %>% top_n(wt=simp_deg, 1) %>% mutate(cond="deg")
fig4B_percolation<-fig4B_avg %>% group_by(vid, netname) %>% top_n(wt=ci_norm, 1) %>% mutate(cond="percolation")
fig4B_eigen<-fig4B_avg %>% group_by(vid, netname) %>% top_n(wt=simp_eigen, 1) %>% mutate(cond="eigen")
fig4B_core<-fig4B_avg %>% group_by(vid, netname) %>% top_n(wt=core, 1) %>% mutate(cond="core")
fig4B_diff<-fig4B_avg %>% group_by(vid, netname) %>% top_n(wt=diffusion_cent, 1) %>% mutate(cond="diffusion_cent")
fig4B_seeds<-rbind(fig4B_btwn, fig4B_close, fig4B_deg, fig4B_eigen, fig4B_diff, fig4B_percolation)
fig4B_random<-fig4B_seeds %>% group_by(vid, netname) %>% sample_n(1) %>% mutate(cond="Random")
fig4B_all<-rbind(fig4B_plc, fig4B_btwn, fig4B_close, fig4B_deg, fig4B_percolation, fig4B_eigen, fig4B_core, fig4B_diff, fig4B_random)
fig4B_all<-subset(fig4B_all, !cond %in% c("diffusion_cent", "close"))
fig4B_all$cond<-as.factor(fig4B_all$cond)
levels(fig4B_all$cond)<-c("Betweenness", "Complex", "K-core", "Degree", "Eigen", "Percolation", "Random")
fig4C<-fig4 %>% group_by(num_seeds, T, vid) %>% mutate(complex_norm=min_max_norm(bw_len), takeup_norm=min_max_norm(takeup_count), takeup_nonleader_norm=min_max_norm(takeup_count_nonleader))

#######################
#Supplementary Figures#
#######################
figs7a_node1<-figs7a_node_data %>% group_by(N, T) %>% dplyr::summarise(prop_steps_over_Wc=sum(avg_step_over_Wc)/length(avg_step_over_Wc) )
figs7a_node2<-figs7a_node_data %>% group_by(N, T) %>% dplyr::summarise_all(funs(mean(.,na.rm=T)))
figs7a_node_mod<-merge(figs7a_node1, figs7a_node2, by=c("N", "T")) 
figs7a_node_mod<-figs7a_node_mod[,names(figs7a_node_mod) %in% c("N", "T", "global_simp_reach", "avg_simp_btwnness", "avg_simp_closeness",
                                                                "avg_clustering", "prop_cascade", "avg_adopt", "avg_adopt_step", "avg_R_by_seed",
                                                                "max_adopt_step", "max_adopt", "bw_len", "avg_bw_size", "avg_simp_path", "prop_steps_over_Wc",
                                                                "simp_btwn","simp_close", "simp_deg", "simp_eigen", "simp_reach", "step_size", "Wc")]

figs7a<-merge(figs7a_cascade, figs7a_node_mod, by=c("N", "T"))
figs7b_agg<-figs7b_node_data %>% group_by(T, sim) %>% dplyr::summarise(avg_step_over_Wc_norm = mean(avg_step_over_Wc_norm),prop_cascade=sum(cascade)/length(cascade) ) 
figs7b<-figs7b_agg %>% group_by(T) %>% dplyr::summarise_all(funs(mean(.,na.rm=T)))

figS8_S9_SI_data_homo<-subset(figs8_s9_main_SI_data, grepl("homo", cond)) 
figS8_S9_SI_data_hetero<-subset(figs8_s9_main_SI_data, grepl("hetero", cond)) 

figS8_S9_data_homo_agg<-figS8_S9_SI_data_homo %>% group_by(T, p, simmod, figure, cond) %>% 
  summarise(avg_bw_size = mean(avg_bw_size, na.rm=T),
            simp_deg = mean(simp_deg, na.rm=T), step_size = mean(step_size, na.rm=T),
            num_adopt = mean(num_adopt, na.rm=T)) %>%
  group_by(T, figure, cond) %>% mutate(num_adopt_norm=min_max_norm(num_adopt), step_size_norm=min_max_norm(step_size) )

figS8_S9_data_hetero_agg<-figS8_S9_SI_data_hetero %>% group_by(p, simmod, figure, cond) %>% 
  summarise(avg_bw_size = mean(avg_bw_size, na.rm=T), simp_deg = mean(simp_deg, na.rm=T),
            step_size = mean(step_size, na.rm=T), num_adopt = mean(num_adopt, na.rm=T)) %>%
  group_by(figure, cond) %>% mutate(num_adopt_norm=min_max_norm(num_adopt), step_size_norm=min_max_norm(step_size));

figS8_S9_data_hetero_agg$T<-NA
figs8_s9_main_SI_data_agg<-rbind(figS8_S9_data_homo_agg,figS8_S9_data_hetero_agg);

figs11_to_s16_main<-subset(figs11_to_s16_SI_data, seed_method != "greedy");
figs11_to_s16_greedy<-subset(figs11_to_s16_SI_data, seed_method == "greedy");
figs11_to_s16_greedy_agg<-figs11_to_s16_greedy %>% group_by(num_seeds, N, k, p, sim, seed_method, ID, m, theta, condition, figure) %>%
  dplyr::summarise(num_adopt = max(num_adopt), btwn_cent=mean(btwn_cent), close_cent=mean(close_cent),
                   eigen_cent=mean(eigen_cent), deg_cent=mean(deg_cent), reach_cent=mean(reach_cent),
                   PlC=mean(PlC)); figs11_to_s16_greedy_agg$T<-"Hetero";
figs11_to_s16_org<-rbind(as.data.frame(figs11_to_s16_main), as.data.frame(figs11_to_s16_greedy_agg))

figs13a<-subset(fig_s13_raw, figure == "s13a")
figs13a_main<-subset(figs13a, seed_method != "greedy")
figs13a_greedy<-subset(figs13a, seed_method == "greedy")
figs13a_greedy_agg<-figs13a_greedy %>% group_by(num_seeds, N, k, p, T, sim, ID, m, seed_method, condition, figure) %>%
  dplyr::summarise(num_adopt = max(num_adopt), btwn_cent=mean(btwn_cent), close_cent=mean(close_cent),
                   eigen_cent=mean(eigen_cent), deg_cent=mean(deg_cent), reach_cent=mean(reach_cent),
                   PlC=mean(PlC)); figs13a_greedy_agg$theta<-NA
figs13a_org<-rbind(as.data.frame(figs13a_main), as.data.frame(figs13a_greedy_agg))
figs13b<-subset(fig_s13_raw, figure == "s13b")
figs13b_org<-figs13b %>% group_by(N, k, p, T, num_seeds, seed_method) %>% dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt, na.rm=T)),
                                                                                           CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt,conf.level = 0.95)$conf.int[2], mean(num_adopt, na.rm=T)), num_adopt = mean(num_adopt, na.rm=T))
figs16a<-subset(figs16ab, figure=="s16a")
figs16a<-subset(figs16a, seed_method %in%  c("btwn", "greedy", "PlC"))
figs16a_org<-figs16a %>% group_by(N, k, p, T, num_seeds, seed_method) %>% dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt,conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                                                                                           CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)), num_adopt = mean(num_adopt))
figs16a_plot<- figs16a_org %>% group_by(N, k, T, num_seeds, seed_method) %>%
  summarise(CI_low = mean(CI_low, na.rm=T), CI_hi = mean(CI_hi, na.rm=T), num_adopt = mean(num_adopt)) %>% group_by(N, k, num_seeds,seed_method) %>%
  summarise(CI_low = mean(CI_low, na.rm=T), CI_hi = mean(CI_hi, na.rm=T), num_adopt = mean(num_adopt,na.rm=T))

figs16b<-subset(figs11_to_s16_SI_data, figure=="s16b")
figs16b<-subset(figs16b, seed_method %in%  c("btwn", "greedy", "PlC"))
figs16b_org<-figs16b %>% group_by(N, k, p, T, num_seeds, seed_method) %>% dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt,conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                                                                                           CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)),num_adopt = mean(num_adopt))
figs16b_plot<-figs16b_org %>% group_by(N, k, num_seeds, seed_method) %>% summarise(CI_low = mean(CI_low), CI_hi = mean(CI_hi), num_adopt = mean(num_adopt))

figs17$num_adopt_mod<-figs17$num_adopt-figs17$num_seeds
figs17$simID<-paste(figs17$sim, figs17$ID, sep="_")
figs17 <- figs17 %>% group_by(num_seeds, m, simID) %>% 
  mutate(num_adopt_norm = min_max_norm(num_adopt_mod), simp_btwn_norm = min_max_norm(simp_btwn), 
         simp_close_norm = min_max_norm(simp_close), simp_deg_norm = min_max_norm(simp_deg), 
         simp_reach_norm = min_max_norm(simp_reach), plcc_norm = min_max_norm(bw_len), 
         simp_eigen_norm = min_max_norm(simp_eigen), core_norm = min_max_norm(core))
figs17_long <- gather(figs17, Centrality, measurement, simp_btwn_norm:simp_eigen_norm, factor_key=TRUE)
figs17_long$centrality<-as.factor(figs17_long$Centrality)
levels(figs17_long$Centrality)<-c("Betweenness", "Closeness", "Degree", "Reach", "Complex", "Eigenvector")
figs17_long<-figs17_long %>% group_by(Centrality, m) %>% mutate(measure_bin=ntile(measurement, 20))
figs17_plot<-figs17_long %>% group_by(Centrality, m, measure_bin) %>% 
  dplyr::summarise(
    cilow = ifelse(var(num_adopt_norm, na.rm=T) != 0, t.test(num_adopt_norm, conf.level = 0.95)$conf.int[1], mean(num_adopt_norm)),
    cihi = ifelse(var(num_adopt_norm, na.rm=T) != 0, t.test(num_adopt_norm, conf.level = 0.95)$conf.int[2], mean(num_adopt_norm)), 
    adoption = mean(num_adopt_norm, na.rm=T) )

figs18a_top_perc<-figs18a_perc %>% group_by(num_seeds, T, netname, N, d) %>% top_n(wt=ci_norm, 1)
figs18a_top_perc$seeding.strat.org<-"Percolation"
figs18a_top_core<-figs18a_main %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=core, 1)
figs18a_top_core$seeding.strat.org<-"K-core"; figs18a_top_core$d<-NA
figs18a_top_btwn<-figs18a_main %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=simp_btwn, 1)
figs18a_top_btwn$seeding.strat.org<-"Betweenness"; figs18a_top_btwn$d<-NA
figs18a_top_eigen<-figs18a_main %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=simp_eigen, 1)
figs18a_top_eigen$seeding.strat.org<-"Eigen."; figs18a_top_eigen$d<-NA
figs18a_top_deg<-figs18a_main %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=simp_deg, 1)
figs18a_top_deg$seeding.strat.org<-"Degree"; figs18a_top_deg$d<-NA
figs18a_top_plcc<-figs18a_main %>% group_by(num_seeds, T, netname, N) %>% top_n(wt=bw_len, 1)
figs18a_top_plcc$seeding.strat.org<-"Complex"; figs18a_top_plcc$d<-NA
figs18a_all<-rbind(figs18a_top_perc, figs18a_top_core, figs18a_top_btwn, figs18a_top_eigen, figs18a_top_deg, figs18a_top_plcc)
figs18a_all$FinalAdopt_mod<-figs18a_all$num_adopt-figs18a_all$num_seeds
figs18a_all<-figs18a_all %>% group_by(netname, T, N) %>% mutate(PLc_norm = min_max_norm(bw_len), adopt_norm=min_max_norm(FinalAdopt_mod))
figs18a_plot<-figs18a_all %>% group_by(seeding.strat.org,netname,N,T,d) %>% dplyr::summarise(adopt_norm = mean(adopt_norm, na.rm=T)) %>% 
  group_by(seeding.strat.org,netname,N,d) %>% dplyr::summarise(adopt_norm = mean(adopt_norm, na.rm=T)) %>% 
  group_by(seeding.strat.org,d) %>% dplyr::summarise(
    cilow=ifelse(var(adopt_norm)==0, mean(adopt_norm), t.test(adopt_norm)$conf.int[1]), 
    cihi=ifelse(var(adopt_norm)==0, mean(adopt_norm), t.test(adopt_norm)$conf.int[2]),
    adopt_norm = mean(adopt_norm, na.rm=T) )
figs18a_plot$d<-as.factor(figs18a_plot$d)

figs18b<-subset(figs18b, takeup ==1)
figs18b_avg_d<- figs18b %>% group_by(d, hhid, vid, netname, takeup, takeup_count_i, takeup_count, takeup_count_nonleader_i, takeup_pct, takeup_pct_nonleader, nrooms, nbeds, electricity, latrine, nroomscapita, nbedscapita, trimester, latrine_backup) %>% 
  dplyr::summarise(simp_btwn = mean(simp_btwn), percolation = mean(percolation), simp_close = mean(simp_close),
                   simp_eigen = mean(simp_eigen), simp_deg = mean(simp_deg), simp_reach = mean(simp_reach), core = mean(core),
                   bw_len = mean(bw_len), diffusion_cent = mean(diffusion_b), ci_norm=mean(ci_norm) )
figs18b_top_percolation<-figs18b_avg_d %>% group_by(d, vid, netname) %>% top_n(wt=ci_norm, 1) %>% mutate(cond="percolation")
figs18b_top_percolation<-figs18b_top_percolation %>% group_by(d,vid,netname) %>% sample_n(1)
figs18b_plot<-figs18b_top_percolation %>% group_by(d) %>% 
  dplyr::summarise(cilow=ifelse(var(takeup_pct)==0, mean(takeup_pct), t.test(takeup_pct)$conf.int[1])/100, 
                   cihi=ifelse(var(takeup_pct)==0, mean(takeup_pct), t.test(takeup_pct)$conf.int[2])/100,
                   takeup_pct = mean(takeup_pct, na.rm=T)/100) 
figs18b_plot$d<-as.factor(figs18b_plot$d)

##############
#PLOT FIGURES#
##############

##############
#Main Figures#
##############

#########
#Fig2ABC#
#########

#Fig2A
fig2A = subset(fig2ABC_final, T == 2)

ggplot(data=fig2A, aes(x=p)) + theme_bw() + theme_bw() +
  geom_point(size=7, aes(y=prop_cascade), shape=16, color="black")	+
  geom_line(size=5,aes(y=prop_cascade), position=position_jitter(w=0.1, h=0))+
  geom_point(size=7, aes(y=prop_steps_over_Wc), shape=16, color="red") + 
  geom_line(size=4.5,aes(y=prop_steps_over_Wc),color="red")+
  labs(y="Frequency", x="p", linetype=NULL) +
  theme(plot.title = element_blank(),axis.title.y=element_text(size = 30, face="bold"),
        axis.title.x=element_text(size = 30, face="bold"), axis.text.x=element_text(size = 40, angle = 0, face="bold"),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black", size=2)) + scale_x_continuous(trans = scales::log2_trans(),breaks = scales::trans_breaks("log2", function(x) 2^x),
                                                                                 labels = scales::trans_format("log2", scales::math_format(2^.x))) +
  theme(axis.text.y=element_text(size = 40, face="bold")) + coord_cartesian(ylim=c(0,1.05))

#fig2B
fig2B = subset(fig2ABC_final, T == 3)

ggplot(data=fig2B, aes(x=p)) + theme_bw() + theme_bw() +
  geom_point(size=7, aes(y=prop_cascade), shape=16, color="black")	+
  geom_line(size=5,aes(y=prop_cascade), position=position_jitter(w=0.1, h=0))+
  geom_point(size=7, aes(y=prop_steps_over_Wc), shape=16, color="red") + 
  geom_line(size=4.5,aes(y=prop_steps_over_Wc),color="red")+
  labs(y="Frequency", x="p", linetype=NULL) +
  theme(plot.title = element_blank(),axis.title.y=element_text(size = 30, face="bold"),
        axis.title.x=element_text(size = 30, face="bold"), axis.text.x=element_text(size = 40, angle = 0, face="bold"),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black", size=2)) +
  scale_x_continuous(trans = scales::log2_trans(),breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x))) +
  theme(axis.text.y=element_text(size = 40, face="bold")) + coord_cartesian(ylim=c(0,1.05))

#fig2C
fig2C = subset(fig2ABC_final, T == 4)

ggplot(data=fig2C, aes(x=p)) + theme_bw() + theme_bw() +
  geom_point(size=7, aes(y=prop_cascade), shape=16, color="black")	+
  geom_line(size=5,aes(y=prop_cascade), position=position_jitter(w=0.1, h=0))+
  geom_point(size=7, aes(y=prop_steps_over_Wc), shape=16, color="red") + geom_line(size=4.5,aes(y=prop_steps_over_Wc),color="red")+
  labs(y="Frequency", x="p", linetype=NULL) +
  theme(plot.title = element_blank(),axis.title.y=element_text(size = 30, face="bold"),
        axis.title.x=element_text(size = 30, face="bold"), axis.text.x=element_text(size = 40, angle = 0, face="bold"),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black", size=2)) + scale_x_continuous(trans = scales::log2_trans(),breaks =
                                                                                   scales::trans_breaks("log2", function(x) 2^x),
                                                                                 labels = scales::trans_format("log2", scales::math_format(2^.x))) +
  theme(axis.text.y=element_text(size = 40, face="bold")) + coord_cartesian(ylim=c(0,1.05))

#fig2D
ggplot(data=fig2D_final, aes(x=p)) + theme_bw() + theme_bw() +
  geom_point(size=7, aes(y=prop_cascade), shape=16, color="black")	+
  geom_line(size=5,aes(y=prop_cascade), position=position_jitter(w=0.1, h=0))+
  geom_point(size=7, aes(y=prop_steps_over_Wc_norm), shape=16, color="red") +
  geom_line(size=4.5,aes(y=prop_steps_over_Wc_norm),color="red")+ labs(y="Frequency", x="p", linetype=NULL) +
  theme(plot.title = element_blank(),axis.title.y=element_text(size = 30, face="bold"),
        axis.title.x=element_text(size = 30, face="bold"), axis.text.x=element_text(size = 40, angle = 0, face="bold"),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black", size=2)) +
  scale_x_continuous(trans = scales::log2_trans(),breaks = scales::trans_breaks("log2", function(x) 2^x),
                     labels = scales::trans_format("log2", scales::math_format(2^.x))) +
  theme(axis.text.y=element_text(size = 40, face="bold")) + coord_cartesian(ylim=c(0,1.05))

#Fig. 2E
fig2e<-subset(fig2_main_SI_data_agg, figure == "2e") #bridge

ggplot(data=fig2e, aes(x=plc_norm, group = cond)) + theme_bw() + geom_point(size=5, aes(y=num_adopt_norm, shape=cond), stroke = 1, alpha=1)	+
  scale_shape_manual(values=c(3,1)) + ylab("Proportion of Adopters") + xlab("Average Bridge Width of Graph") + 
  theme(plot.title = element_text(size=28, hjust = 0.5),axis.title.x=element_blank(),
        axis.title.y=element_blank(),axis.text.x=element_text(size = 28, angle = 0),
        axis.text.y=element_text(size = 28),legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black") )	+ 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), labels =c(0,0.2,0.4,0.6,0.8,1))

cor.test(fig2e$plc_norm, fig2e$plc_norm, method="spearman")
spearman.ci(fig2e$plc_norm, fig2e$plc_norm, conf.level = 0.95)

cor.test(subset(fig2e, cond == "homo_sf_fract")$plc_norm, subset(fig2e, cond == "homo_sf_fract")$num_adopt_norm, method="spearman")
spearman.ci(subset(fig2e, cond == "homo_sf_fract")$plc_norm, subset(fig2e, cond == "homo_sf_fract")$num_adopt_norm,conf.level = 0.95)
cor.test(subset(fig2e, cond == "hetero_sf_fract")$plc_norm, subset(fig2e, cond == "hetero_sf_fract")$num_adopt_norm, method="spearman")
spearman.ci(fig2e$plc_norm, fig2e$num_adopt_norm, conf.level = 0.95)

#Bridge width version
ggplot(data=fig2e, aes(x=step_size_norm, group = cond)) + theme_bw() + geom_point(size=5, aes(y=num_adopt_norm, shape=cond), stroke = 1, alpha=1)	+
  scale_shape_manual(values=c(3,1)) + ylab("Proportion of Adopters") + xlab("Average Bridge Width of Graph") + 
  theme(plot.title = element_text(size=28, hjust = 0.5),axis.title.x=element_blank(),
        axis.title.y=element_blank(),axis.text.x=element_text(size = 28, angle = 0),
        axis.text.y=element_text(size = 28),legend.position = "none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black") )	+ 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), labels =c(0,0.2,0.4,0.6,0.8,1))

cor.test(fig2e$step_size_norm, fig2e$num_adopt_norm, method="spearman")
spearman.ci(fig2e$step_size_norm, fig2e$num_adopt_norm, conf.level = 0.95)

cor.test(subset(fig2e, cond == "homo_sf_fract")$step_size_norm, subset(fig2e, cond == "homo_sf_fract")$num_adopt_norm, method="spearman")
spearman.ci(subset(fig2e, cond == "homo_sf_fract")$step_size_norm, subset(fig2e, cond == "homo_sf_fract")$num_adopt_norm,conf.level = 0.95)
cor.test(subset(fig2e, cond == "hetero_sf_fract")$step_size_norm, subset(fig2e, cond == "hetero_sf_fract")$num_adopt_norm, method="spearman")
spearman.ci(fig2e$step_size_norm, fig2e$num_adopt_norm, conf.level = 0.95)

#############
#Fig3ABCDEFG#
#############
pairwise.wilcox.test(fig3_all$num_adopt, fig3_all$seeding.strat.org)
fig3_plot_final$Centrality.Seeding.Strategy<-fig3_plot_final$seeding.strat.org
levels(fig3_plot_final$Centrality.Seeding.Strategy)<-c("Eigenvector", "Degree", "Betweenness", "K-core", "Percolation", "Complex")

#Fig. 3A
ggplot(data=fig3_plot_final, aes(x=PLc_norm, y=adopt_norm)) + theme_bw() +
  geom_jitter(size=9, aes(fill=Centrality.Seeding.Strategy),colour="black",pch=21, stroke = 0.5, shape=16, width=0.3, height=0.05) +
  scale_fill_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  ylab("Proportion of Adopters") + xlab("Complex Centrality") + 
  theme(plot.title = element_text(size=50, hjust =  0.5),axis.title.y=element_text(size=50),
        axis.title.x=element_text(size=70), 
        axis.text.x=element_text(size = 50, angle = 0),
        axis.text.y=element_text(size = 50), legend.text=element_text(size=30),legend.title=element_text(size=30),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black") )	+ scale_y_continuous(breaks=seq(0,1, 0.2)) + 
  scale_x_continuous(trans = scales::log2_trans(), breaks = scales::trans_breaks("log10", function(x) 10^x), labels =	scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(breaks = seq(0,1,0.2),labels = seq(0,1,0.2))

#Fig. 3B
fig3B<-fig3_plot_final %>% group_by(seeding.strat.org) %>% summarise(CI_lower = t.test(Btwn.cent, conf.level = 0.95)$conf.int[1], CI_upper = t.test(Btwn.cent, conf.level = 0.95)$conf.int[2], 
                                                                   Btwn.cent=mean(Btwn.cent,na.rm=T))
fig3B$seeding.strat.org <- factor(fig3B$seeding.strat.org, levels = c("Top.Eigen", "Top.Deg","Top.Btwn", "Top.Core", "Top.Perc.", "Top.Complex"))

ggplot(fig3B, aes(x=seeding.strat.org, y=Btwn.cent, fill=seeding.strat.org)) +
  geom_bar(stat="identity", color="black", size=2, width=1, alpha=0.8) + 
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), linetype="solid",width=0.3, size=3)+
  scale_fill_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.text.y = element_text(size=70),
                     axis.title.x = element_blank(),axis.title.y = element_text(size=70),
                     strip.text.x = element_text(size=40),
                     legend.position = "none", 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Betweenness Centrality", linetype=NULL) +
  scale_y_continuous(breaks=seq(0,0.08,0.02), labels = seq(0,0.08,0.02)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  coord_cartesian(ylim=c(0.01,0.1))

#Fig 3C
fig3C<-fig3_plot_final %>% group_by(seeding.strat.org) %>% 
  summarise(CI_lower = ifelse(var(Degree_norm)==0, mean(Degree_norm, na.rm=T), t.test(Degree_norm, conf.level = 0.95)$conf.int[1]), 
            CI_upper = ifelse(var(Degree_norm)==0, mean(Degree_norm, na.rm=T), t.test(Degree_norm, conf.level = 0.95)$conf.int[2]), 
            Degree_norm=mean(Degree_norm, na.rm=T))
fig3C$seeding.strat.org <- factor(fig3C$seeding.strat.org, levels = c("Top.Eigen", "Top.Deg","Top.Btwn", "Top.Core", "Top.Perc.", "Top.Complex"))

ggplot(fig3C, aes(x=seeding.strat.org, y=Degree_norm, fill=seeding.strat.org)) +
  geom_bar(stat="identity", color="black", size=2, width=1, alpha=0.8) + 
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), linetype="solid",width=0.3, size=3)+
  scale_fill_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=70), axis.title.x = element_blank(), axis.title.y = element_text(size=70), 
        strip.text.x = element_text(size=30), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +labs(y="Degree (Normalized)", linetype=NULL) +
  scale_y_continuous(breaks=seq(0.2,0.95,0.2), labels=seq(0.2,0.95,0.2)) + coord_cartesian(ylim=c(0.3,1)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

#Fig 3D
fig3D<-fig3_plot_final %>% group_by(seeding.strat.org) %>%  summarise( CI_lower = t.test(simp_eigen, conf.level = 0.95)$conf.int[1], 
                                                                       CI_upper = t.test(simp_eigen, conf.level = 0.95)$conf.int[2], Eigen.cent=mean(simp_eigen, na.rm=T))
fig3D$seeding.strat.org <- factor(fig3D$seeding.strat.org, levels = c("Top.Eigen", "Top.Deg","Top.Btwn", "Top.Core", "Top.Perc.", "Top.Complex"))

ggplot(fig3D, aes(x=seeding.strat.org, y=Eigen.cent, fill=seeding.strat.org)) +
  geom_bar(stat="identity", color="black", size=2, width=1, alpha=0.8) + 
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), linetype="solid", width=0.3, size=3)+
  scale_fill_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  theme_bw() + theme(axis.text.x = element_blank(), legend.position = "none",
                     axis.text.y = element_text(size=70), axis.title.x = element_blank(),
                     axis.title.y = element_text(size=70), strip.text.x = element_text(size=30), #legend.position = "none",
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Eigenvector Centrality", linetype=NULL) +
  scale_y_continuous(labels=seq(0,0.28, 0.06), breaks=seq(0,0.28, 0.06)) + 
  coord_cartesian(ylim=c(0.06,0.26)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

#Fig 3E
fig3E<-fig3_plot_final %>% group_by(seeding.strat.org) %>% summarise(CI_lower = ifelse(var(ci_norm)==0, mean(ci_norm), t.test(ci_norm, conf.level = 0.95)$conf.int[1]), 
                                                                     CI_upper = ifelse(var(ci_norm)==0, mean(ci_norm), t.test(ci_norm, conf.level = 0.95)$conf.int[2]), 
                                                                     percolation=mean(ci_norm, na.rm=T))
fig3E$seeding.strat.org <- factor(fig3E$seeding.strat.org, levels = c("Top.Eigen", "Top.Deg","Top.Btwn", "Top.Core", "Top.Perc.", "Top.Complex"))

ggplot(fig3E, aes(x=seeding.strat.org, y=percolation, fill=seeding.strat.org)) +
  geom_bar(stat="identity", color="black", size=2, width=1, alpha=0.8) + 
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), linetype="solid", width=0.3, size=3)+
  scale_fill_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  theme_bw() + theme(axis.text.x = element_blank(), legend.position = "none",
                     axis.text.y = element_text(size=70), axis.title.x = element_blank(),
                     axis.title.y = element_text(size=70), strip.text.x = element_text(size=30), #legend.position = "none",
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Percolation Centrality", linetype=NULL) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  coord_cartesian(ylim=c(0.28,1))

#Fig 3F
fig3F<-fig3_plot_final %>% group_by(seeding.strat.org) %>% 
  summarise(CI_lower = ifelse(var(core_norm, na.rm=T)==0, mean(core_norm, na.rm=T), t.test(core_norm, conf.level = 0.95)$conf.int[1]), 
            CI_upper = ifelse(var(core_norm, na.rm=T)==0, mean(core_norm, na.rm=T), t.test(core_norm, conf.level = 0.95)$conf.int[2]), 
            core_norm=mean(core_norm, na.rm=T))
fig3F$seeding.strat.org <- factor(fig3F$seeding.strat.org, levels = c("Top.Eigen", "Top.Deg","Top.Btwn", "Top.Core", "Top.Perc.", "Top.Complex"))

ggplot(fig3F, aes(x=seeding.strat.org, y=core_norm, fill=seeding.strat.org)) +
  geom_bar(stat="identity", color="black", size=2, width=1, alpha=0.8) + 
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), linetype="solid", width=0.3, size=3)+
  scale_fill_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  theme_bw() + theme(axis.text.x = element_blank(), legend.position="none",
                     axis.text.y = element_text(size=70), axis.title.x = element_blank(),
                     axis.title.y = element_text(size=70), strip.text.x = element_text(size=30), #legend.position = "none",
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="K-core (Normalized)", linetype=NULL) + 
  coord_cartesian(ylim=c(0.25,1)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

#Fig 3
fig3G<-fig3_plot_final %>% group_by(seeding.strat.org) %>% 
  summarise(CI_lower = ifelse(var(PLc_norm, na.rm=T)==0, mean(PLc_norm, na.rm=T), t.test(PLc_norm, conf.level = 0.95)$conf.int[1]), 
            CI_upper = ifelse(var(PLc_norm, na.rm=T)==0, mean(PLc_norm, na.rm=T), t.test(PLc_norm, conf.level = 0.95)$conf.int[2]), 
            PLc_norm=mean(PLc_norm, na.rm=T))
fig3G$seeding.strat.org <- factor(fig3G$seeding.strat.org, levels = c("Top.Eigen", "Top.Deg","Top.Btwn", "Top.Core", "Top.Perc.", "Top.Complex"))

ggplot(fig3G, aes(x=seeding.strat.org, y=PLc_norm, fill=seeding.strat.org)) +
  geom_bar(stat="identity", color="black", size=2, width=1, alpha=0.8) + 
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), linetype="solid", width=0.3, size=3)+
  scale_fill_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  theme_bw() + theme(axis.text.x = element_blank(),legend.position = "none",
                     axis.text.y = element_text(size=70), axis.title.x = element_blank(),
                     axis.title.y = element_text(size=70), strip.text.x = element_text(size=30), #legend.position = "none",
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Complex Centrality", linetype=NULL) + 
  coord_cartesian(ylim=c(0.28,1)) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

######
#Fig4#
######

#Fig4A#
fig4A_all$Probability.Adoption.All.Neighbors<-fig4A_all$takeup_pct
fig4A_all$Probability.Adoption.NonLeaders<-fig4A_all$takeup_pct_nonleader
fig4A_all$Seeding.Strategy<-fig4A_all$cond
fig4A_all$vid<-as.factor(fig4A_all$vid)
fig4A_all<-within(fig4A_all, Seeding.Strategy <- relevel(Seeding.Strategy, ref = c("Random")))
fig4A_all$Probability.Adoption.All.Neighbors<-fig4A_all$Probability.Adoption.All.Neighbors/100
fig4A_all_s<-fig4A_all %>% group_by(vid, netname, Seeding.Strategy) %>% sample_n(1)

fig4A_model <- lm(Probability.Adoption.All.Neighbors ~ Seeding.Strategy + nrooms + nbeds + electricity + latrine + nroomscapita+ nbedscapita + trimester + latrine_backup + vid, data=fig4A_all_s)
fig4A_coeffs<-data.frame(fig4A_model$coefficients)
fig4A_coeffs$measure<-rownames(fig4A_coeffs)
colnames(fig4A_coeffs)<-c("coeffs", "measure")
fig4A_CIs<-data.frame(confint(fig4A_model))
fig4A_CIs$measure<-rownames(fig4A_CIs)
colnames(fig4A_CIs)<-c("CI.Lower", "CI.Upper", "measure")
fig4A_lm_data<-merge(fig4A_coeffs, fig4A_CIs, by="measure")
fig4A_intercept<-subset(fig4A_lm_data, measure=="(Intercept)")$coeffs
fig4A_plot<-subset(fig4A_lm_data, measure %in% c("Seeding.StrategyComplex", "Seeding.StrategyDegree", "Seeding.StrategyEigen", "Seeding.StrategyK-core", 
                                                 "Seeding.StrategyPercolation", "Seeding.StrategyBetweenness"))
fig4A_plot$measure<-as.factor(fig4A_plot$measure)
levels(fig4A_plot$measure)<-c("Betweenness","Complex", "Degree", "Eigenvector", "K-core", "Percolation")
fig4A_plot$measure<-factor(fig4A_plot$measure, levels=c("Eigenvector", "Degree", "Betweenness", "K-core", "Percolation", "Complex"))

ggplot(fig4A_plot, aes(x=reorder(measure, coeffs), y=coeffs, fill=measure, color=measure)) +
  geom_point(size=30) + geom_errorbar(aes(ymin=CI.Lower, ymax=CI.Upper), linetype="solid",width=0, size=8)+
  scale_color_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  ylab("P(Peer Adopts)") + ggtitle("Seeding with Leader Households") + 
  theme_bw() + theme(axis.text.y = element_text(size=80),
                     #axis.text.x = element_text(size=50, angle=25, vjust=0.6), 
                     axis.text.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=80),
                     strip.text.x = element_text(size=0),
                     plot.title = element_text(size=72, hjust=0.5), 
                     legend.position = "none", 
                     legend.title = element_blank(),
                     legend.text = element_text(size=20), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  geom_hline(yintercept = 0, linetype="solid", size=3) 

#Supplementary Table S2
summary(fig4A_model)
tab_model(fig4A_model)
fig4A_model_vcov <- cluster.vcov(fig4A_model, fig4A_all_s$vid) #clustered standard errors 
coeftest(fig4A_model, fig4A_model_vcov)

#Fig4B#
fig4B_all$Probability.Adoption.All.Neighbors<-fig4B_all$takeup_pct
fig4B_all$Probability.Adoption.NonLeaders<-fig4B_all$takeup_pct_nonleader
fig4B_all$Seeding.Strategy<-fig4B_all$cond
fig4B_all$vid<-as.factor(fig4B_all$vid)
fig4B_all<-within(fig4B_all, Seeding.Strategy <- relevel(Seeding.Strategy, ref = c("Random")))
fig4B_all$Probability.Adoption.All.Neighbors<-fig4B_all$Probability.Adoption.All.Neighbors/100
fig4B_all_s<-fig4B_all %>% group_by(vid, netname, Seeding.Strategy) %>% sample_n(1)
fig4B_model <- lm(Probability.Adoption.All.Neighbors ~ Seeding.Strategy + nrooms + nbeds + electricity + latrine + 
                    nroomscapita+ nbedscapita + trimester + latrine_backup + vid, data=fig4B_all_s)
fig4B_coeffs<-data.frame(fig4B_model$coefficients)
fig4B_coeffs$measure<-rownames(fig4B_coeffs)
colnames(fig4B_coeffs)<-c("coeffs", "measure")
fig4B_CIs<-data.frame(confint(fig4B_model))
fig4B_CIs$measure<-rownames(fig4B_CIs)
colnames(fig4B_CIs)<-c("CI.Lower", "CI.Upper", "measure")
fig4B_lm_data<-merge(fig4B_coeffs, fig4B_CIs, by="measure")
fig4B_intercept<-subset(fig4B_lm_data, measure=="(Intercept)")$coeffs
fig4B_plot<-subset(fig4B_lm_data, measure %in% c("Seeding.StrategyComplex", "Seeding.StrategyDegree", "Seeding.StrategyEigen", "Seeding.StrategyK-core", 
                                                 "Seeding.StrategyPercolation", "Seeding.StrategyBetweenness"))
fig4B_plot$measure<-as.factor(fig4B_plot$measure)
levels(fig4B_plot$measure)<-c("Betweenness","Complex", "Degree", "Eigenvector", "K-core", "Percolation")
fig4B_plot$coeffs<-fig4B_plot$coeffs+fig4B_intercept
fig4B_plot$CI.Lower<-fig4B_plot$CI.Lower+fig4B_intercept
fig4B_plot$CI.Upper<-fig4B_plot$CI.Upper+fig4B_intercept
fig4B_plot$measure<-factor(fig4B_plot$measure, levels=c("Eigenvector", "Degree", "Betweenness", "K-core", "Percolation", "Complex"))

ggplot(fig4B_plot, aes(x=reorder(measure, coeffs), y=coeffs, fill=measure, color=measure)) +
  geom_point(size=30) + geom_errorbar(aes(ymin=CI.Lower, ymax=CI.Upper), linetype="solid",width=0, size=8)+
  scale_color_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  ylab("P(Peer Adopts)") + ggtitle("Seeding with Leader Households") + 
  theme_bw() + theme(axis.text.y = element_text(size=80),
                     #axis.text.x = element_text(size=50, angle=25, vjust=0.6), 
                     axis.text.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size=80),
                     strip.text.x = element_text(size=0),
                     plot.title = element_text(size=72, hjust=0.5), 
                     legend.position = "none", 
                     legend.title = element_blank(),
                     legend.text = element_text(size=20), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  geom_hline(yintercept = fig4B_intercept, linetype="solid", size=3) 

#Supplementary Table S3
summary(fig4B_model)
tab_model(fig4B_model)
fig4B_model_vcov <- cluster.vcov(fig4B_model, fig4B_all_s$vid) #clustered standard errors 
coeftest(fig4B_model, fig4B_model_vcov)

#Fig. 4C#
fig4C_plot<- fig4C %>% group_by(vid) %>% dplyr::summarise(takeup_norm = mean(takeup_norm,na.rm=T), takeup_nonleader_norm = mean(takeup_nonleader_norm, na.rm=T),avg_PLcc = mean(complex_norm,na.rm=T))
cor.test(fig4C$bw_len, fig4C$complex_norm)
cor.test(fig4C_plot$takeup_norm, fig4C_plot$avg_PLcc)
cor.test(fig4C_plot$takeup_nonleader_norm, fig4C_plot$avg_PLcc)

fig4C_plot$Complex.Cent.Bin<-ntile(fig4C_plot$avg_PLcc, 5)
fig4C_plot_bin<-fig4C_plot %>% group_by(Complex.Cent.Bin) %>% dplyr::summarise(avg.frac.adopt=mean(takeup_norm, na.rm=T))

ggplot(fig4C_plot_bin, aes(x=Complex.Cent.Bin, y=avg.frac.adopt)) +
  geom_point(size=13, alpha=0.8) + 
  xlab("Average Complex Path Length\n of Village (Binned)") + 
  ylab("Proportion of Adopters") + ggtitle("Seeding with any Village Member") + 
  theme_bw() + theme(axis.text.x = element_text(size=40, vjust=0.6), 
                     axis.text.y = element_text(size=40),
                     axis.title.x = element_text(size=40),
                     axis.title.y = element_text(size=40),
                     strip.text.x = element_text(size=40),
                     plot.title = element_blank(), 
                     legend.position = c(0.75,0.15), 
                     legend.title = element_blank(),
                     legend.text = element_text(size=20),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))


#######################
#Supplementary Figures#
#######################

#Fig s7a
ggplot(data=figs7a, aes(x=T, group=1)) + theme_bw() + theme_bw() + geom_point(size=5, aes(y=prop_cascade), shape=16, color="black") +
  geom_line(size=3,aes(y=prop_cascade))+ geom_point(size=5, aes(y=prop_steps_over_Wc), shape=16, color="red") + 
  geom_line(size=3,aes(y=prop_steps_over_Wc), color = "red")+
  ylab("Frequency") + xlab("Threshold") + ggtitle(paste("Global Diffusion", sep="")) +
  theme(plot.title = element_blank(), axis.title.y=element_text(size = 30, face="bold"),
        axis.title.x=element_text(size = 30, angle = 0, face="bold"), axis.text.x=element_text(size = 32, angle = 0, face="bold"),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black", size=2)) + theme(axis.text.y=element_text(size = 32, face="bold")) +
  scale_x_continuous(breaks = seq(1,7,1))

#Fig s7b
ggplot(data=figs7b, aes(x=T, group=1)) + theme_bw() + theme_bw() + geom_point(size=5, aes(y=prop_cascade), shape=16, color="black") +
  geom_line(size=3,aes(y=prop_cascade))+ geom_point(size=5, aes(y=avg_step_over_Wc_norm), shape=16, color="red") +	
  geom_line(size=3,aes(y=avg_step_over_Wc_norm), color = "red")+ ylab("Frequency") + xlab("Threshold (Fractional)") +
  ggtitle(paste("Global Diffusion", sep="")) + 
  theme(plot.title = element_blank(), axis.title.y=element_text(size = 40, face="bold"),
        axis.title.x=element_text(size = 40, angle = 0, face="bold"), axis.text.x=element_text(size = 32, angle = 0, face="bold"),
        legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black", size=2)) + theme(axis.text.y=element_text(size = 32, face="bold")) +
  scale_x_continuous(breaks = seq(0.1,1,0.2)) + coord_cartesian(ylim=c(0,1.05))

#Fig s8
figs8<-subset(figs8_s9_main_SI_data_agg, figure == "s8")

ggplot(data=figs8, aes(x=step_size_norm, group = cond)) + theme_bw() + 
  geom_point(size=5, aes(y=num_adopt_norm, shape=cond), stroke = 1, alpha=1)	+
  scale_shape_manual(values=c(3,1)) +
  ylab("Proportion of Adopters") + 
  xlab("Average Bridge Width of Graph") + 
  theme(plot.title = element_text(size=28, hjust = 0.5),axis.title.x=element_text(size=28, hjust = 0.5), 
        axis.title.y=element_text(size=28, hjust = 0.5),axis.text.x=element_text(size = 28, angle = 0),
        axis.text.y=element_text(size = 28),legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour =	"black") ) + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), labels =c(0,0.2,0.4,0.6,0.8,1))

#Fig s9
figs9<-subset(figs8_s9_main_SI_data_agg, figure == "s9")

ggplot(data=figs9, aes(x=step_size_norm, group = cond)) + theme_bw() + 
  geom_jitter(size=5, aes(y=num_adopt_norm, shape=cond), stroke = 1, alpha=1, height=0,width=0.01)	+
  scale_shape_manual(values=c(9,5)) + 
  ylab("Proportion of Adopters") + 
  xlab("Average Bridge Width of Graph\n(Normalized)") + 
  ggtitle("k-regular Graphs (k=8)") + 
  theme(plot.title = element_text(size=28, hjust = 0.5), axis.title.x=element_text(size = 28),
        axis.title.y=element_text(size = 28),axis.text.x=element_text(size = 28, angle = 0),
        axis.text.y=element_text(size = 28),legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour =	"black"))	+ 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), labels =c(0,0.2,0.4,0.6,0.8,1))

#Fig s10
figs10<-fig3_plot %>% group_by(N, T, num_seeds, seeding.strat.org) %>% sample_n(1)
figs10_plot<-figs10 %>% group_by(T, seeding.strat.org) %>% 
  dplyr::summarise(
    adopt_avg = mean(adopt_norm, na.rm=T), 
    cilow= ifelse(var(adopt_norm, na.rm=T)==0, mean(adopt_norm, na.rm=T), t.test(adopt_norm)$conf.int[1]), 
    cihi= ifelse(var(adopt_norm, na.rm=T)==0, mean(adopt_norm, na.rm=T), t.test(adopt_norm)$conf.int[2]))

figs10_plot$T<-as.factor(figs10_plot$T)
figs10_plot$seeding.strat.org<-as.factor(figs10_plot$seeding.strat.org)
figs10_plot$Threshold<-figs10_plot$T

ggplot(figs10_plot, 
       aes(x=Threshold, y=adopt_avg, fill=seeding.strat.org)) +
  geom_bar(stat="identity", color="black", size=1, width=0.8, alpha=0.8, position=position_dodge(0.8)) + 
  geom_errorbar(aes(ymin=cilow, ymax=cihi), linetype="solid",width=0.3, size=1, position=position_dodge(0.8))+
  scale_fill_manual(values = c("indianred3", "olivedrab1","orange", "grey", "darkslateblue", "cadetblue2")) +
  theme_bw() + theme(axis.text.x = element_text(size=50),
                     axis.text.y = element_text(size=50),
                     axis.title.x = element_text(size=50),
                     axis.title.y = element_text(size=50),
                     strip.text.x = element_text(size=40),
                     legend.text = element_text(size=20),
                     legend.position = c(0.85,0.85), 
                     legend.title=element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Proportion of Adopters", linetype=NULL)

#Table S1
levels(figs10$seeding.strat.org)<-c("Eigenvector", "Degree", "Betweenness", "K-core", "Percolation", "Complex")
figs10$Seeding.Strategy<-figs10$seeding.strat.org
figs10$Threshold<-figs10$T

figs10$Proportion.Adopters<-figs10$num_adopt/figs10$N
#figs10$Proportion.Adopters<-figs10$adopt_norm
table_S1<-lm(Proportion.Adopters ~ Threshold + Seeding.Strategy, data = figs10)
summary(table_S1)
table_S1_vcov <- cluster.vcov(table_S1, figs10$net)
coeftest(table_S1, table_S1_vcov)
#library(sjPlot)
#tab_model(table_S1)

#Fig S11A
figS11A_plot<-subset(figs11_to_s16_org, figure == "s11a") %>% group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, 
                                   t.test(num_adopt, conf.level = 0.95)$conf.int[1], 
                                   mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, 
                                  t.test(num_adopt, 
                                         conf.level = 0.95)$conf.int[2], 
                                  mean(num_adopt)),num_adopt = mean(num_adopt))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){ 
  new_row<-data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0,seed_method=currmethod , figure="2a")
  figS11A_plot<-rbind(as.data.frame(figS11A_plot),new_row)
}

figS11A_plot$num_seeds_perc<-as.numeric(as.character(figS11A_plot$num_seeds/figS11A_plot$N))
figS11A_plot$seed_method<-as.factor(figS11A_plot$seed_method)
levels(figS11A_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(subset(figS11A_plot, num_seeds<= 30), 
       aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000, ymax=CI_hi/1000, shape =	seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+ scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values =c("black", "black", "black", "black", "red")) +
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.8,0.2),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.03,0.01), labels=scales::percent_format(trim=TRUE,accuracy=1)) + 
  scale_y_continuous(breaks=c(0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(0,0.7))

#Fig S11B
figS11B_plot<-subset(figs11_to_s16_org, figure == "s11b") %>% 
  group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)),
                   num_adopt = mean(num_adopt))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){ 
  new_row<-data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0, seed_method=currmethod, figure="2b")
  figS11B_plot<-rbind(as.data.frame(figS11B_plot),new_row)
}

figS11B_plot$num_seeds_perc<-as.numeric(as.character(figS11B_plot$num_seeds/figS11B_plot$N))
figS11B_plot$targets<-figS11B_plot$N-figS11B_plot$num_seeds
figS11B_plot$prop_adopt<-as.numeric(as.character(figS11B_plot$num_adopt/figS11B_plot$targets))
figS11B_plot$num_seeds_perc<-as.numeric(as.character(figS11B_plot$num_seeds/figS11B_plot$N))
figS11B_plot$CI_low_prop_adopt<-figS11B_plot$CI_low/figS11B_plot$targets 
figS11B_plot$CI_hi_prop_adopt<-figS11B_plot$CI_hi/figS11B_plot$targets 
figS11B_plot<-as.data.frame(figS11B_plot)

figS11B_plot$seed_method<-as.factor(figS11B_plot$seed_method)
levels(figS11B_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(figS11B_plot, aes(x = num_seeds_perc, y = prop_adopt, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4,aes(ymin=CI_low_prop_adopt,ymax=CI_hi_prop_adopt, shape = seed_method,color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+
  scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values=c("black", "black", "black", "black", "red"))+ 
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.8,0.2),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.03,0.01),labels=scales::percent_format(trim=TRUE,accuracy=1)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), labels=c(0,0.2,0.4,0.6,0.8,1)) + 
  coord_cartesian(ylim=c(0,1))

#Fig 11C
figS11C_plot<-subset(figs11_to_s16_org, figure == "s11c") %>% group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.99)$conf.int[1], mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.99)$conf.int[2], mean(num_adopt)),
                   num_adopt = mean(num_adopt))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){
  new_row<-data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0, seed_method=currmethod, figure="2c")
  figS11C_plot<-rbind(as.data.frame(figS11C_plot),new_row)
}

figS11C_plot$num_seeds_perc<-as.numeric(as.character(figS11C_plot$num_seeds/figS11C_plot$N))
figS11C_plot$seed_method<-as.factor(figS11C_plot$seed_method)
levels(figS11C_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(subset(figS11C_plot, num_seeds<= 20),
       aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000, ymax=CI_hi/1000, shape = seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+
  scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values =c("black", "black", "black", "black", "red")) +
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.8,0.2),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.03,0.01),labels=scales::percent_format(trim=TRUE,accuracy=1)) 

#Fig 11D
figS11D_plot<-subset(figs11_to_s16_org, figure == "s11d") %>% group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)),
                   num_adopt = mean(num_adopt))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){ 
  new_row<-data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0,seed_method=currmethod,	figure="2d")
  figS11D_plot<-rbind(as.data.frame(figS11D_plot),new_row)
}

figS11D_plot$num_seeds_perc<-as.numeric(as.character(figS11D_plot$num_seeds/figS11D_plot$N))
figS11D_plot$seed_method<-as.factor(figS11D_plot$seed_method)
levels(figS11D_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(subset(figS11D_plot, num_seeds<= 30),
       aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000, ymax=CI_hi/1000, shape = seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method)) +
  scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values =c("black", "black", "black", "black", "red")) +
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.8,0.2),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.03,0.01), labels=scales::percent_format(trim=TRUE,accuracy=1))

#Fig s12a
figs12a_plot<-subset(figs11_to_s16_org, figure == "s12a") %>% group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)),
                   num_adopt = mean(num_adopt))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){
  new_row<-data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0, seed_method=currmethod, figure="s6a")
  figs12a_plot<-rbind(as.data.frame(figs12a_plot),new_row)
}

figs12a_plot$num_seeds_perc<-as.numeric(as.character(figs12a_plot$num_seeds/figs12a_plot$N))
figs12a_plot$seed_method<-as.factor(figs12a_plot$seed_method)
levels(figs12a_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(subset(figs12a_plot, num_seeds<= 20),
       aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000,ymax=CI_hi/1000, shape = seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+ scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values=c("black", "black", "black", "black", "red"))+ xlab("Proportion of Seeds") + 
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.8,0.2),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.03,0.01), labels=scales::percent_format(trim=TRUE,accuracy=1)) 

#Fig s12b
figs11_to_s16_org$num_adopt_mod<-figs11_to_s16_org$num_adopt-figs11_to_s16_org$num_seeds
figs12b_plot<-subset(figs11_to_s16_org, figure == "s12b") %>% group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt_mod) != 0, t.test(num_adopt_mod, conf.level = 0.95)$conf.int[1], mean(num_adopt_mod)),
                   CI_hi = ifelse(var(num_adopt_mod) != 0, t.test(num_adopt_mod, conf.level = 0.95)$conf.int[2], mean(num_adopt_mod)),
                   num_adopt = mean(num_adopt_mod))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){
  new_row<-data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0, seed_method=currmethod, figure="s6b")
  figs12b_plot<-rbind(as.data.frame(figs12b_plot),new_row)
}

figs12b_plot$num_seeds_perc<-as.numeric(as.character(figs12b_plot$num_seeds/figs12b_plot$N))
figs12b_plot$seed_method<-as.factor(figs12b_plot$seed_method)
levels(figs12b_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(subset(figs12b_plot, num_seeds<= 20), aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000,ymax=CI_hi/1000, shape = seed_method, color = seed_method)) + 
  geom_line(size=1.2, aes(color=seed_method))+ scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values=c("black", "black", "black", "black", "red"))+ 
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.8,0.2),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.03,0.01), labels=scales::percent_format(trim=TRUE,accuracy=1))

#fig s13a
figs13a_plot<-figs13a_org %>% group_by(N, k, p, T, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)),
                   num_adopt = mean(num_adopt))

for(currT in c(2, 3, 4)){
  for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){
    new_row<-data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,T = currT, num_adopt=0, seed_method=currmethod, figure="s7a")
    figs13a_plot<-rbind(as.data.frame(figs13a_plot),new_row)
  }
}

figs13a_plot<-figs13a_plot %>% group_by(N, k, p, num_seeds, seed_method) %>% 
  dplyr::summarise(CI_low = mean(CI_low, na.rm=T), CI_hi = mean(CI_hi, na.rm=T), num_adopt = mean(num_adopt, na.rm=T))

figs13a_plot$num_seeds_perc<- as.numeric(as.character(figs13a_plot$num_seeds/figs13a_plot$N))
figs13a_plot$seed_method<-as.factor(figs13a_plot$seed_method)
levels(figs13a_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(subset(figs13a_plot, num_seeds<= 12), 
       aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000,ymax=CI_hi/1000, shape = seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+ scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values=c("black", "black", "black", "black", "red"))+ 
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.85,0.2),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.012,0.002), labels=scales::percent_format(trim=TRUE,accuracy=0.1))

#fig s13b
figs13b_plot<-figs13b_org %>% group_by(N, k, p, num_seeds, seed_method) %>% 
  dplyr::summarise(CI_low = mean(CI_low), CI_hi = mean(CI_hi), num_adopt = mean(num_adopt))

figs13b_plot$seed_method<-as.factor(figs13b_plot$seed_method) 
figs13b_plot$num_seeds_perc<-figs13b_plot$num_seeds/figs13b_plot$N
figs13b_plot<-as.data.frame(figs13b_plot)

figs13b_plot$seed_method<-as.factor(figs13b_plot$seed_method)
levels(figs13b_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(figs13b_plot,aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000,ymax=CI_hi/1000, shape = seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+ scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values=c("black", "black", "black", "black", "red"))+ 
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.2,0.8),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 25, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.03,0.003), labels=scales::percent_format(trim=TRUE,accuracy=0.1))

#Fig s14a
figs14a_plot<-subset(figs11_to_s16_org, figure == "s14a") %>% group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)),
                   num_adopt = mean(num_adopt))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){
  new_row<-data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0,seed_method=currmethod, figure="s8a")
  figs14a_plot<-rbind(as.data.frame(figs14a_plot),new_row)
}

figs14a_plot$num_seeds_perc<-as.numeric(as.character(figs14a_plot$num_seeds/figs14a_plot$N))
figs14a_plot$seed_method<-as.factor(figs14a_plot$seed_method)
levels(figs14a_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(subset(figs14a_plot, num_seeds<= 15),aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000, ymax=CI_hi/1000, shape =	seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+
  scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values =c("black", "black", "black", "black", "red")) +
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.8,0.2),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.015,0.002), labels=scales::percent_format(trim=TRUE,accuracy=.1)) + 
  coord_cartesian(ylim=c(525/1000,670/1000))

#Fig s14b
figs14b_plot<-subset(figs11_to_s16_org, figure == "s14b") %>% group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)),
                   num_adopt = mean(num_adopt))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){ 
  new_row<- data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0, seed_method=currmethod, figure="s8b")
  figs14b_plot<-rbind(as.data.frame(figs14b_plot),new_row)
}

figs14b_plot$num_seeds_perc<-as.numeric(as.character(figs14b_plot$num_seeds/figs14b_plot$N))
figs14b_plot$seed_method<-as.factor(figs14b_plot$seed_method)
levels(figs14b_plot$seed_method)<-c("Betweenness", "Degree", "Eigenvector", "Greedy", "Complex")

ggplot(subset(figs14b_plot, num_seeds<= 15), aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000, ymax=CI_hi/1000, shape =	seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+
  scale_shape_manual(values =c(15, 1, 17, 18, 19)) +
  scale_color_manual(values =c("black", "black", "black", "black", "red")) +
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.15,0.83),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.015,0.005),labels=scales::percent_format(trim=TRUE,accuracy=.1)) + 
  coord_cartesian(ylim=c(800/1000,860/1000))

#Figs15a
figs15a_plot<-subset(figs11_to_s16_org, figure == "s15a") %>% group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)),
                   num_adopt = mean(num_adopt))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){ 
  new_row<- data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0,seed_method=currmethod,	figure="s9a")
  figs15a_plot<-rbind(as.data.frame(figs15a_plot),new_row)
}

figs15a_plot$num_seeds_perc<-as.numeric(as.character(figs15a_plot$num_seeds/figs15a_plot$N))
figs15a_plot$seed_method<-as.factor(figs15a_plot$seed_method)
levels(figs15a_plot$seed_method)<-c("Betweenness", "Closeness", "Degree", "Eigenvector", "Greedy", "Complex", "Reach")

ggplot(subset(figs15a_plot, num_seeds<= 30),aes(x = num_seeds_perc, y = num_adopt/1000, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low/1000, ymax=CI_hi/1000, shape =	seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+
  scale_shape_manual(values =c(15, 5, 1, 17, 18, 19, 22)) +
  scale_color_manual(values =c("black", "black", "black", "black", "black", "red", "black")) +
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.8,0.25),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.03,0.01),labels=scales::percent_format(trim=TRUE,accuracy=1)) + 
  coord_cartesian(ylim=c(0,0.7))

#Figs15b
figs15b_plot<-subset(figs11_to_s16_org, figure == "s15b") %>% group_by(N, k, p, num_seeds, seed_method, figure) %>%
  dplyr::summarise(CI_low = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                   CI_hi = ifelse(var(num_adopt) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)), num_adopt = mean(num_adopt))

for(currmethod in c("btwn", "deg", "eigen", "greedy", "PlC")){ 
  new_row<-data.frame(N=1000,k=8,p=0.5,num_seeds=0,CI_low=0,CI_hi=0,num_adopt=0, seed_method=currmethod, figure="s9b")
  figs15b_plot<-rbind(as.data.frame(figs15b_plot),new_row)
}

figs15b_plot$num_seeds_perc<-as.numeric(as.character(figs15b_plot$num_seeds/figs15b_plot$N))
figs15b_plot$targets<-figs15b_plot$N-figs15b_plot$num_seeds
figs15b_plot$prop_adopt<-as.numeric(as.character(figs15b_plot$num_adopt/figs15b_plot$targets))
figs15b_plot$num_seeds_perc<-as.numeric(as.character(figs15b_plot$num_seeds/figs15b_plot$N))
figs15b_plot$CI_low_prop_adopt<-figs15b_plot$CI_low/figs15b_plot$targets
figs15b_plot$CI_hi_prop_adopt<-figs15b_plot$CI_hi/figs15b_plot$targets
figs15b_plot<-as.data.frame(figs15b_plot)

figs15b_plot$seed_method<-as.factor(figs15b_plot$seed_method)
levels(figs15b_plot$seed_method)<-c("Betweenness", "Closeness", "Degree", "Eigenvector", "Greedy", "Complex", "Reach")

ggplot(figs15b_plot, aes(x = num_seeds_perc, y = prop_adopt, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low_prop_adopt,ymax=CI_hi_prop_adopt, shape = seed_method,color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+
  scale_shape_manual(values =c(15, 5, 1, 17, 18, 19, 22)) +
  scale_color_manual(values =c("black", "black", "black", "black", "black", "red", "black")) +
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.8,0.2),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.03,0.01), labels=scales::percent_format(trim=TRUE,accuracy=1)) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), labels =c(0,0.2,0.4,0.6,0.8,1)) + 
  coord_cartesian(ylim=c(0,1))

#Figs16a
figs16a_plot$num_seeds_perc<-figs16a_plot$num_seeds/figs16a_plot$N
figs16a_plot<-as.data.frame(figs16a_plot)
figs16a_plot<-figs16a_plot[complete.cases(figs16a_plot),]
figs16a_plot<-subset(figs16a_plot, seed_method %in%  c("btwn", "greedy", "PlC"))

for(currmethod in c("btwn", "greedy", "PlC")){ 
  new_row<-data.frame(N=1000,k=8,num_seeds=0,seed_method=currmethod,CI_low=0,CI_hi=0,num_adopt=0, num_seeds_perc=0)
  figs16a_plot<-rbind(as.data.frame(figs16a_plot),new_row)
}

figs16a_plot$num_adopt_norm<-min_max_norm(figs16a_plot$num_adopt)
figs16a_plot$CI_low_norm<-min_max_norm(figs16a_plot$CI_low)
figs16a_plot$CI_hi_norm<-min_max_norm(figs16a_plot$CI_hi)
figs16a_plot$seed_method<-as.factor(figs16a_plot$seed_method)
levels(figs16a_plot$seed_method)<-c("Betweenness", "Greedy", "Complex")

ggplot(subset(figs16a_plot, num_seeds_perc<=0.01), aes(x = num_seeds_perc, y = num_adopt_norm, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low_norm,ymax=CI_hi_norm, shape = seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+
  scale_shape_manual(values =c(15, 18, 19)) +
  scale_color_manual(values=c("black", "black", "red"))+
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.15,0.8),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.01,0.002), labels=scales::percent_format(trim=TRUE,accuracy=0.1))
  
#fig s16b
figs16b_plot$num_seeds_perc<-figs16b_plot$num_seeds/figs16b_plot$N 
figs16b_plot<-as.data.frame(figs16b_plot)
for(currmethod in c("btwn", "greedy", "PlC")){ 
  new_row<-data.frame(N=1000,k=8,num_seeds=0,seed_method=currmethod,CI_low=0,CI_hi=0,num_adopt=0,num_seeds_perc=0)
  figs16b_plot<-rbind(as.data.frame(figs16b_plot),new_row)
}

figs16b_plot<-subset(figs16b_plot, N==1000 & num_seeds<=12)
figs16b_plot$seed_method<-as.factor(figs16b_plot$seed_method)
levels(figs16b_plot$seed_method)<-c("Betweenness", "Greedy", "Complex")

figs16b_plot$num_adopt_norm<-min_max_norm(figs16b_plot$num_adopt)
figs16b_plot$CI_low_norm<-min_max_norm(figs16b_plot$CI_low)
figs16b_plot$CI_hi_norm<-min_max_norm(figs16b_plot$CI_hi)

ggplot(figs16b_plot, aes(x = num_seeds_perc, y = num_adopt_norm, group = seed_method)) + theme_bw() +
  geom_pointrange(size = 1.4, aes(ymin=CI_low_norm,ymax=CI_hi_norm, shape = seed_method, color = seed_method)) +
  geom_line(size=1.2, aes(color=seed_method))+
  scale_shape_manual(values =c(15, 18, 19)) +
  scale_color_manual(values=c("black", "black", "red"))+
  xlab("Proportion of Network as Seeds") + ylab("Proportion of Adopters") + theme(
    axis.title.x=element_text(size=30),
    axis.title.y=element_text(size=30),
    legend.text=element_text(size=20),
    legend.title=element_blank(),
    legend.position=c(0.15,0.8),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black",fill="white", size=1.4),
    axis.text.x=element_text(size = 30, vjust=0.8),axis.text.y=element_text(size = 30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = seq(0,0.012,0.002), labels=scales::percent_format(trim=TRUE,accuracy=0.1))

#fig. S17 (Standard scale-free graph)
figs17<-figs17 %>% group_by(Centrality, m) %>% mutate(measure_bin=ntile(measurement, 20))

figs17_plot<-figs17 %>% group_by(Centrality, m, measure_bin) %>% 
  dplyr::summarise(cilow = ifelse(var(num_adopt, na.rm=T) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[1], mean(num_adopt)),
                   cihi = ifelse(var(num_adopt, na.rm=T) != 0, t.test(num_adopt, conf.level = 0.95)$conf.int[2], mean(num_adopt)), 
                   adoption = mean(num_adopt, na.rm=T) )
figs17_plot$Centrality<-factor(figs17_plot$Centrality, levels=c("Complex", "Reach", "Degree", "Betweenness", "Eigenvector", "Closeness"))

#modify m to view robustness 
ggplot(data=subset(figs17_plot, m==15 & Centrality != "Closeness"), aes(x=measure_bin, y=adoption/1000, color=Centrality, group=Centrality, shape=Centrality)) + 
  theme_bw() + geom_point(size=5) +
  geom_line(size=1.2, aes(color=Centrality))+ scale_shape_manual(values = c(16, 18, 1, 15,17)) +
  scale_color_manual(values =c("red", "black", "black", "black", "black")) +
  geom_errorbar(aes(ymin=cilow/1000, ymax=cihi/1000), linetype="solid",width=0, size=1)+
  xlab("Centrality (Binned)") + ylab("Proportion of Adopters") + 
  theme_bw() + theme(axis.text.x = element_text(size=50),axis.text.y = element_text(size=50),
                     axis.title.x = element_text(size=50),axis.title.y = element_text(size=50),
                     strip.text.x = element_text(size=50),legend.position = c(0.25,0.78), 
                     legend.text=element_text(size=20),legend.title=element_text(size=20, face="bold"),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))


#fig. S18a
ggplot(subset(figs18a_plot, seeding.strat.org=="Percolation"), aes(x=d, y=adopt_norm)) +
  geom_bar(stat="identity", color="black", size=2, width=1, alpha=0.8) + 
  geom_errorbar(aes(ymin=cilow, ymax=cihi), linetype="solid",width=0.3, size=3)+
  labs(x="d\n(Percolation Seeding)") + ggtitle("Add Health (Simulated)") + 
  theme_bw() + theme(axis.text.x = element_text(size=50),axis.text.y = element_text(size=50),
                     axis.title.x = element_text(size=50),axis.title.y = element_text(size=50),
                     strip.text.x = element_text(size=50),plot.title = element_text(size=50, hjust=0.5),legend.position = "none", 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Proportion of Adopters", linetype=NULL) + coord_cartesian(ylim=c(0,0.4))

#fig. S18b
ggplot(figs18b_plot, aes(x=d, y=takeup_pct)) +
  geom_bar(stat="identity", color="black", size=2, width=1, alpha=0.8) + 
  geom_errorbar(aes(ymin=cilow, ymax=cihi), linetype="solid",width=0.3, size=3)+
  labs(x="d\n(Percolation Seeding)") + ggtitle("BSS Program (Empirical)") + 
  theme_bw() + theme(axis.text.x = element_text(size=50),axis.text.y = element_text(size=50),
                     axis.title.x = element_text(size=50),axis.title.y = element_text(size=50),
                     strip.text.x = element_text(size=50),plot.title = element_text(size=50, hjust=0.5),
                     legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y="Proportion of Adopters", linetype=NULL) + coord_cartesian(ylim=c(0,0.4))


