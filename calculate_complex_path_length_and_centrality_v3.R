#Topological Measures for Identifying and Predicting the Spread of Complex Contagions
#Guilbeault & Centola, 2021, Nature Communications 
#June 2021 

#This script calculates average complex path length for all nodes in a given graph, g
#which provides a measure of complex centrality; methods w/wo parallel processing are provided 
#Pipeline is provided for replicated main analyses from Guilbeault & Centola, 2021
#Using both simulated graphs and empirical Addhealth networks 

rm(list=ls());gc()
library(dplyr)
library(tidyr)
library(influential)
library(fastnet)
library(igraph)
library(doParallel)
library(parallel)

#Model Functions
min_max_norm<-function(x){(x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x,na.rm=TRUE))}

clustered_seeding<-function(seeds, g, seeds_needed){
  possible_seeds<-unique(unlist(sapply(seeds, function(x) neighbors(g, x, mode = "total"))))
  seeds_to_add<-possible_seeds[!possible_seeds %in% seeds]
  need_seeds<-seeds_needed > length(seeds_to_add)
  if(need_seeds){return(possible_seeds)}
  else{final_seeds<-c(seeds, sample(seeds_to_add, seeds_needed))
  return(final_seeds)
  }
}

get_simple_centralities<-function(g){
  centrality_df<-data.frame(seed=as.numeric(V(g)), degree = as.numeric(degree(g)), betweenness = as.numeric(betweenness(g)), eigen = as.numeric(eigen_centrality(g)$vector) ) 
  percolation<-collective.influence(graph=g, vertices = V(g), mode="all", d=3)
  centrality_df$percolation<-as.numeric(percolation); 
  return(centrality_df)
}

get_complex<-function(seed, N, g, gmat, thresholds, num_seeds_to_add, model_output_list){
  gmat_simulation<-matrix(nrow=N,ncol=N,0)
  num_seeds_i<-num_seeds_to_add[seed]
  seeds_to_add<-numeric(num_seeds_i)
  
  if(num_seeds_i > 0){
    seeds<-as.numeric(neighbors(g, seed, mode = "total")) 
    num_seeds_needed<-num_seeds_i - length(seeds)
    need_seeds<-num_seeds_needed > 0 
    if(need_seeds){#initiate clustered Seeding
      seeds<-clustered_seeding(seeds, g, num_seeds_needed)
      num_seeds_needed<-num_seeds_i - length(seeds)
      need_seeds<-num_seeds_needed > 0}else{seeds<-c(seed, sample(seeds,num_seeds_i))}
  }else{seeds=seed}
  
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
  all_distances<-distances(complex_g, seed, V(complex_g), mode="out")
  all_distances[all_distances == "Inf"]<-0
  PLci<-mean(all_distances)
  
  model_output_list<-c(seed, N, num_seeds_i, num_adopters, PLci)

  return(model_output_list)
}

###########################
#####Replicate Results#####
###With Parallel Methods###
###########################
detectCores(logical = TRUE)
cluster <- makeCluster(8)  
registerDoParallel(cluster)  
clusterExport(cluster,list('neighbors', 'graph_from_adjacency_matrix', 'distances', 'V', 'clustered_seeding'))

#Replicate with simulated graphs 
#num_graphs<-20; N<-200
#graphs<-c(); for(i in 1:num_graphs){graphs[[i]]<-to.igraph(net.holme.kim(N,4,0.5))} #use Holme & Kim SF algorithm
#gmats<-c(); for(i in 1:num_graphs){gmats[[i]]<-as.matrix(as_adjacency_matrix(graphs[[i]]))}

#Replicate with empirical Addhealth networks 
addhealth_mats<-list.files("C:/Users/dougl/Desktop/AddHealth_Networks_Largest_Components/") #DL zipped folder from github and set directory 
gmats<-c(); for(i in 1:length(addhealth_mats)){gmats[[i]]<-read.csv(paste("C:/Users/dougl/Desktop/AddHealth_Networks_Largest_Components/", addhealth_mats[i], sep=""))}
for(i in 1:length(gmats)){colnames(gmats[[i]])<-1:length(gmats[[i]]);rownames(gmats[[i]])<-1:length(gmats[[i]])}  
for(i in 1:length(gmats)){gmats[[i]]<-as.matrix(gmats[[i]])}  
graphs<-c(); for(i in 1:length(addhealth_mats)){graphs[[i]]<-graph_from_adjacency_matrix(gmats[[i]])} 

#Set threshold type
T_dist<-c("homo", "hetero")[1]
T_type<-c("abs","frac")[1]

comparison_df<-data.frame()

for(i in 1:length(graphs)){ #
  print(paste("Graph: ", i, sep=""))
  g<-graphs[[i]]
  gmat<-gmats[[i]]
  N<-nrow(gmat)
  
  if(length(V(g))>10){ #block errors from empty networks 
    thresholds<-replicate(N, 3) #thresholds<-replicate(N, runif(1, 0.1,0.5)) #for Fractional Thresholds
    if(T_type == "frac"){thresholds = round(unlist(lapply(1:N, function(x) thresholds[x] * length(neighbors(g, x, mode = "total")))))}
    num_seeds_to_add<-thresholds-1
    num_seeds_to_add[num_seeds_to_add<0]=0
    thresholds[thresholds<=0]=1
    simple_centralities_df<-get_simple_centralities(g)
    
    model_output_list <- vector(mode = "list", length = N)
    start_time <- proc.time() #run Model
    model_output_list <- clusterMap(cluster, get_complex, seed=1:N, MoreArgs=list(N=N, g=g,gmat=gmat,thresholds=thresholds,num_seeds_to_add=num_seeds_to_add,model_output_list=model_output_list)) 
    print(proc.time() - start_time) #view runtime 
    model_output_list_df <- as.data.frame(do.call('rbind', model_output_list))
    colnames(model_output_list_df)<-c( "seed","N","num_seeds", "num_adopters","PLci")
    model_output_list_df$PLci_norm<-min_max_norm(model_output_list_df$PLci)
    model_output_parallel_full<-merge(model_output_list_df, simple_centralities_df, by="seed") #map to simple centralities of seeds
    model_output_parallel_full$threshold<-thresholds; model_output_parallel_full$T_dist<-T_dist; model_output_parallel_full$T_type<-T_type
    
    top_plc<-model_output_parallel_full %>% top_n(1, wt = PLci_norm) %>% mutate(top="Complex") %>% sample_n(1)
    top_btwn<-model_output_parallel_full %>% top_n(1, wt = betweenness) %>% mutate(top="betweenness") %>% sample_n(1)
    top_degree<-model_output_parallel_full %>% top_n(1, wt = degree) %>% mutate(top="degree") %>% sample_n(1)
    top_eigen<-model_output_parallel_full %>% top_n(1, wt = eigen) %>% mutate(top="eigen") %>% sample_n(1)
    top_percolation<-model_output_parallel_full %>% top_n(1, wt = percolation) %>% mutate(top="percolation") %>% sample_n(1)
    top_final<-rbind(top_plc, top_btwn, top_degree, top_eigen, top_percolation)
    top_final$graph<-i
    
    comparison_df<-rbind(comparison_df, top_final)
  }
}

pairwise.wilcox.test(comparison_df$num_adopters, comparison_df$top, p.adjust.method ="none", paired=T)
mean(subset(comparison_df, top=="Complex")$num_adopters)
mean(subset(comparison_df, top=="betweenness")$num_adopters)
mean(subset(comparison_df, top=="degree")$num_adopters)
mean(subset(comparison_df, top=="eigen")$num_adopters)
mean(subset(comparison_df, top=="percolation")$num_adopters)


#####################################
##Model without Parallel Processing##
#####################################
N<-300
g<-net.holme.kim(N,4,0.5)
g<-to.igraph(g)
gmat<-as.matrix(as_adjacency_matrix(g))
T_dist<-c("homo", "hetero")[2]
T_type<-c("abs","frac")[2]
thresholds<-replicate(N, runif(1, 0.1,0.5)) #thresholds<-replicate(N, 8)
if(T_type == "frac"){thresholds = round(unlist(lapply(1:N, function(x) thresholds[x] * length(neighbors(g, x, mode = "total")))))}
num_seeds_to_add<-thresholds-1
num_seeds_to_add[num_seeds_to_add<0]=0
thresholds[thresholds<=0]=1
simple_centralities_df<-get_simple_centralities(g)

#Run Model
model_output_list <- vector(mode = "list", length = N)
start_time <- proc.time()
model_output_list<-lapply(1:N, function(x) get_complex(x, N, g, gmat, thresholds, num_seeds_to_add, model_output_list)) #Run model; change 1:N to any subset of specific seeds to narrow search
print(proc.time() - start_time) #view runtime 
model_output_df <- as.data.frame(Reduce(rbind, model_output_list))
colnames(model_output_df)<-c("seed","N","num_seeds", "num_adopters","PLci")
model_output_df$PLci_norm<-min_max_norm(model_output_df$PLci)
model_output_full<-merge(model_output_df, simple_centralities_df, by="seed") #map to simple centralities of seed
