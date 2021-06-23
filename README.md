# complexpaths
Guilbeault, D. & D. Centola. 2021. Topological Measures for Identifying and Predicting the Spread of Complex Contagions. Nature Communications. 

This repository contains: 
-	Raw edgelists of Addhealth topologies (82 adolescent friend networks).
-	The largest connected component extracted from each Addhealth topology. 
-	An R script that calculates complex path length and complex centrality for all nodes in a given graph g (this includes methods for parallelization). 
-	An R script that replicates all of the figures from the manuscript and supplementary material using the raw simulation data generated for this experiment (data provided below).

Features in development: 
- A Python code base for reproducing all of the measures defined in Guilbeault & Centola (2021). In the long term, these measures may be released as a package in Python and/or R.

Data: 
- Here is a link to all of the raw data used to generate the main and supplementary figures for this study: https://drive.google.com/file/d/1M9sO51Viv9h6e-p83Q6FThOiYZZpdUhJ/view?usp=sharing
- The topologies and household-level statistics for the Banerjee et al. (2013) analysis were downloaded from the following data repository (https://osf.io/qu57w/); this accessed repository used the Banerjee et al. (2013) for “Study 2” in the following article: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0245476 

Any questions regarding this code or data can be sent to douglas.guilbeault[at]haas.berkeley.edu. 
