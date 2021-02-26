# flowMapper


## Installation
It's easiest to install flowMapper with the help of the `devtools` package.

`install.packages("devtools") # if necessary`  
`devtools::install_github("matei-ionita/flowMapper")`

## Summary
Mapper is an algorithm inspired by topological data analysis, originally due to [[1]](#1). It constructs a graph which approximates the (hypothetical) space from which data points are sampled. This implementation makes design choices which are suitable for flow cytometry data: for example, only using subroutines that run in linear time, so that millions of data points can be processed in a few minutes.

Each node in the Mapper graph represents multiple data points, resulting in a compact visualization of the dataset. But the graph is more than the visualization. For example, we use the Leiden method[[2]](#2) to detect communities in the graph, which can then be interpreted as cell phenotypes.

![Leiden figure](leidenCommunities.png)


## References
<a id="1">[1]</a> 
Singh, Gurjeet, Facundo Mémoli, and Gunnar E. Carlsson. "Topological methods for the analysis of high dimensional data sets and 3d object recognition." SPBG 91 (2007): 100.

<a id="2">[2]</a> 
Traag, Vincent A., Ludo Waltman, and Nees Jan Van Eck. "From Louvain to Leiden: guaranteeing well-connected communities." Scientific reports 9.1 (2019): 1-12.
