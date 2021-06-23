#Topological Measures for Identifying and Predicting the Spread of Complex Contagions
#Guilbeault & Centola, 2021, Nature Communications 
#June 2021 

#This script calculates average complex path length for all nodes in a given graph, g
#which provides a measure of complex centrality 
#methods with and without parallel processing are provided 

rm(list=ls());gc()
library(igraph)
library(dplyr)
library(tidyr)
library(influential)
library(doParallel)
library(parallel)

#Model Functions
min_max_norm<-function(x){(x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x,na.rm=TRUE))}

run_sim<-function(seed, N, g, gmat, thresholds, num_seeds_to_add, model_output_list){
  gmat_simulation<-matrix(nrow=N,ncol=N,0)
  num_seeds_i<-num_seeds_to_add[seed]
  seeds_to_add<-numeric(num_seeds_i)
  
  if(num_seeds_i > 0){
    possible_seeds<-neighbors(g, seed, mode = "total")  
    need_seeds<-num_seeds_i > length(possible_seeds) #clustered Seeding
    if(need_seeds){
      possible_seeds<-unique(unlist(sapply(possible_seeds, function(x) neighbors(g, x, mode = "total"))))
      need_seeds<-num_seeds_i > length(possible_seeds) }
    seeds_to_add<-sample(possible_seeds, num_seeds_i) }
  
  seeds<-c(seed, seeds_to_add)
  activated<-logical(N); activated[seeds]<-TRUE
  gmat_simulation[seeds,]<-1; gmat_simulation[,seeds]<-1
  gmat_simulation_run<-gmat*gmat_simulation
  
  spread=TRUE
  while(spread){
    influence<-colSums(gmat_simulation_run)
    influence_activated<-influence>=thresholds
    t<-which(activated)
    t_1<-which(influence_activated)
    t_full<-union(t,t_1)
    spread<-length(t_full)>length(t)
    activated[t_full]<-TRUE
    adopters<-which(activated)
    gmat_simulation[adopters,]<-1
    gmat_simulation[,adopters]<-1
    gmat_simulation_run<-gmat*gmat_simulation
  }

  num_adopters<-sum(activated)
  complex_g<-graph_from_adjacency_matrix(gmat_simulation_run)
  PLci<-mean_distance(complex_g, directed = FALSE)

  model_output_list<-c(seed, N, num_seeds_i, num_adopters, PLci)

  return(model_output_list)
}

########################
#Setup Model Parameters#
########################
N<-300
g <- sample_smallworld(1, N, 4, 0.1)
gmat<-as.matrix(as_adjacency_matrix(g))

#T_dist<-"homo"
#T_type<-"abs"
#thresholds<-replicate(N, 6)
#num_seeds_to_add<-thresholds-1

T_dist<-"hetero"
T_type<-"frac" 
thresholds<-replicate(N, runif(1, 0.1,0.5)) 
if(T_type == "frac"){thresholds = round(unlist(lapply(1:N, function(x) thresholds[x] * length(neighbors(g, x, mode = "total")))))}

num_seeds_to_add<-thresholds-1
num_seeds_to_add[num_seeds_to_add<0]=0
thresholds[thresholds<=0]=1

#Get extant centrality measures for each node 
centrality_df<-data.frame(seed=1:N, degree = as.numeric(degree(g)), betweenness = as.numeric(betweenness(g)), 
                          closeness = as.numeric(closeness(g)), eigen = as.numeric(eigen_centrality(g)$vector),
                          pagerank = as.numeric(page_rank(g)$vector) )

percolation<-collective.influence(graph=g, vertices = V(g), mode="all", d=3)
centrality_df$percolation<-as.numeric(percolation); centrality_df$percolation_norm<-min_max_norm(percolation)

#############################
#########Run Model########### 
#Without Parallel Processing#
#############################

model_output_list <- vector(mode = "list", length = N)
start_time <- proc.time()
model_output_list<-lapply(1:N, function(x) run_sim(x, N, g, gmat, thresholds, num_seeds_to_add, model_output_list)) #Run model; change 1:N to any subset of specific seeds to narrow search
print(proc.time() - start_time) #view runtime 
model_output_df <- as.data.frame(Reduce(rbind, model_output_list))
colnames(model_output_df)<-c("seed","N","num_seeds", "num_adopters","PLci")
model_output_df$PLci_norm<-min_max_norm(model_output_df$PLci)
model_output_full<-merge(model_output_df, centrality_df, by="seed") #map to seed centrality 

#Basic Stats 
cor.test(model_output_full$PLci_norm, model_output_full$num_adopters)
cor.test(model_output_full$degree, model_output_full$num_adopters)
cor.test(model_output_full$eigen, model_output_full$num_adopters)
cor.test(model_output_full$betweenness, model_output_full$num_adopters)
cor.test(model_output_full$percolation_norm, model_output_full$num_adopters)

##########################
#######Run Model########## 
#With Parallel Processing#
##########################

#Initialize clusters
detectCores(logical = TRUE)
cluster <- makeCluster(8)  
registerDoParallel(cluster)  
clusterExport(cluster,list('neighbors', 'graph_from_adjacency_matrix', 'get.shortest.paths', 'mean_distance'))

#Run Parallelized Model
model_output_list <- vector(mode = "list", length = N)
start_time <- proc.time()
model_output_list <- clusterMap(cluster, run_sim, seed=1:N, MoreArgs=list(N=N, g=g,gmat=gmat,thresholds=thresholds,num_seeds_to_add=num_seeds_to_add,
                                                                          model_output_list=model_output_list)) #change 1:N to any subset of specific seeds to narrow search
print(proc.time() - start_time) #view runtime 

model_output_list_df <- as.data.frame(do.call('rbind', model_output_list))
colnames(model_output_list_df)<-c( "seed","N","num_seeds", "num_adopters","PLci")
model_output_list_df$PLci_norm<-min_max_norm(model_output_list_df$PLci)
model_output_parallel_full<-merge(model_output_list_df, centrality_df, by="seed") #map to seed centrality 

#Basic Stats 
cor.test(model_output_parallel_full$PLci_norm, model_output_parallel_full$num_adopters)
cor.test(model_output_parallel_full$degree, model_output_parallel_full$num_adopters)
cor.test(model_output_parallel_full$eigen, model_output_parallel_full$num_adopters)
cor.test(model_output_parallel_full$betweenness, model_output_parallel_full$num_adopters)
cor.test(model_output_parallel_full$percolation_norm, model_output_parallel_full$num_adopters)





