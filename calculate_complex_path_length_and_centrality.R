#Topological Measures for Identifying and Predicting the Spread of Complex Contagions
#Guilbeault & Centola, 2021, Nature Communications 
#June 2021 

#This script calculates average complex path length for all nodes in a given graph, g
#which then provides a measure of complex centrality 
#parallel processing is used to improve speed performance in large networks 
#increasing the number of available cores increases speed 

rm(list=ls());gc()
library(igraph)
library(dplyr)
library(tidyr)
library(MASS)

run_sim<-function(seed, N, g, gmat, gmat_simulation, T_dist, T_type, thresholds, sim_df){
  num_seeds=thresholds[seed]-1
  seeds<-c(seed, sample(as.numeric(neighbors(g, seed, mode = "total")), num_seeds))
  activated<-replicate(N, FALSE); activated[seeds]<-TRUE
  gmat_simulation[seeds,]<-1; gmat_simulation[,seeds]<-1
  gmat_simulation_run<-gmat*gmat_simulation

  spread=TRUE
  step<-1
  
  while(spread){
    influence<-colSums(gmat_simulation_run)
    influence_activated<-influence>=thresholds
    t<-which(activated)
    t_1<-which(influence_activated)
    t_full<-union(t,t_1)
    spread<-length(t_full)>length(t)
    activated[t_full]<-TRUE
    step<-step+1
    adopters<-which(activated)
    gmat_simulation[adopters,]<-1
    gmat_simulation[,adopters]<-1
    gmat_simulation_run<-gmat*gmat_simulation
  }
  
  num_adopters<-sum(activated)
  complex_g<-graph_from_adjacency_matrix(gmat_simulation_run)
  all_simple_paths<-get.shortest.paths(complex_g, seed, mode = "out")
  PLci<-mean(sapply(1:N, function(x) length(all_simple_paths$vpath[[x]])))
  
  sim_df<-rbind(sim_df, data.frame(seed = seed, N=N, adopters = num_adopters, Thresh = thresholds[seed], 
                                   T_dist = T_dist, T_type = T_type, num_seeds = num_seeds, PLci = PLci))
  return(sim_df)
}


########################
#Setup Model Parameters#
########################

N<-500
g <- sample_smallworld(1, N, 5, 0.1)
gmat<-as.matrix(as_adjacency_matrix(g))
gmat_simulation<-gmat * 0
T_dist<-"homo"
T_type<-"abs"
thresholds<-replicate(N, 4)
hist(degree(g))

centrality_df<-data.frame(seed=1:N, degree = as.numeric(degree(g)), betweenness = as.numeric(betweenness(g)), 
                          closeness = as.numeric(closeness(g)), eigen = as.numeric(eigen_centrality(g)$vector),
                          pagerank = as.numeric(page_rank(g)$vector) )


#############################
#########Run Model########### 
#Without Parallel Processing#
#############################

model_output<-data.frame()

start_time <- proc.time()
for(x in 1:N){model_output<-run_sim(x, N, g, gmat, gmat_simulation, T_dist, T_type, thresholds, model_output)}
print(proc.time() - start_time)

model_output_full<-merge(model_output, centrality_df, by="seed") #map to seed centrality 
cor.test(model_output_full$PLci, model_output_full$adopters)
cor.test(model_output_full$degree, model_output_full$adopters)
cor.test(model_output_full$betweenness, model_output_full$adopters)
cor.test(model_output_full$closeness, model_output_full$adopters)
cor.test(model_output_full$eigen, model_output_full$adopters)
cor.test(model_output_full$pagerank, model_output_full$adopters)

##########################
#######Run Model########## 
#With Parallel Processing#
##########################

library(doParallel)
library(parallel)

#setup cluster
detectCores(logical = TRUE)
cluster <- makeCluster(8)  
registerDoParallel(cluster)  
clusterExport(cluster,list('neighbors', 'graph_from_adjacency_matrix', 'get.shortest.paths'))

model_output_parallel_raw<-data.frame()
start_time <- proc.time()
model_output_parallel_raw <- clusterMap(cluster, run_sim, seed=1:N, 
                                        MoreArgs=list(N=N, g=g,gmat=gmat,
                                                      gmat_simulation=gmat_simulation,T_dist=T_dist,
                                                      T_type=T_type,thresholds=thresholds,sim_df=model_output_parallel_raw))
print(proc.time() - start_time)

model_output_parallel <- do.call('rbind', model_output_parallel_raw)

model_output_parallel_full<-merge(model_output_parallel, centrality_df, by="seed") #map to seed centrality 
cor.test(model_output_parallel_full$PLci, model_output_parallel_full$adopters)
cor.test(model_output_parallel_full$degree, model_output_parallel_full$adopters)
cor.test(model_output_parallel_full$betweenness, model_output_parallel_full$adopters)
cor.test(model_output_parallel_full$closeness, model_output_parallel_full$adopters)
cor.test(model_output_parallel_full$eigen, model_output_parallel_full$adopters)
cor.test(model_output_parallel_full$pagerank, model_output_parallel_full$adopters)







