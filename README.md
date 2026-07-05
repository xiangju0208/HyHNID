# HyHNID: Hypergraph-coupled heterogeneous network impulsive dynamics for predicting disease-related genes

Identifying disease-associated genes accurately is essential for exploring disease mechanisms and precision medicine. Most traditional network-based methods only capture pairwise interactions between biological nodes, ignoring their high-order synergistic relationships, resulting in unsatisfactory prediction performance. Therefore, this study proposes a novel computational method based on hypergraph-coupled heterogeneous network impulsive dynamics (HyHNID) for disease gene prediction. It first constructs four types of multi-view hyperedges based on biological networks of genes and diseases, and a cross-layer similarity calculation strategy is used to optimize network topology and reduce structural noises. Unlike conventional methods that rely solely on pairwise diffusion, HyHNID builds a high-order coupled dynamic system by integrating four types of hyperedge-coupling dynamics into network impulsive dynamical equations of genes and diseases, respectively. This system can simultaneously extract basic pairwise features and high-order synergistic features but also realize the cross-layer signal propagation in heterogeneous biological networks. We conduct comprehensive evaluations including five-fold cross-validation, novel disease testing and independent testing. Experimental results substantiate the positive contributions of hyperedge-coupling dynamics and demonstrate that HyHNID outperforms both pairwise diffusion models and other state-of-the-art algorithms. Consequently, the proposed method offers an effective and robust framework for disease-gene prediction and helps better understand the synergistic regulation of complex biological networks. 

<div align="center">
  <img src="Fig.1-HyHNID.png" width="80%" />
</div>

## Requirements
Matlab 2016 or above   


## Codes 
#main_HyHNID.m: cross-validation code.   
This code allows parallel execution.   <br>
 <br>
#A_HyHNID.m: the recommended algorithm in the study. <br>
% Input:  <br>
% AdjGfG: associatins between Genes (G) and Genes (G)   <br>
% AdjGfD: associatins between Genes (G) and Diseases (D)  <br>
% AdjDfD: associatins between Diseases (D) and Diseases (G)  <br>
% P0_G: initials in Gene network  <br>
% P0_D: initials in Disease network  <br> 
% Ouput: <br>
% TableScores: a table whos variable record the scores of genes. <br>
  

## Dataset 
This dataset includes: <br>
 disease-gene associations, disease-disease associations and gene-gene associations.  <br>
 


## Results 
The results will be automatically saved into the directory: results.  

 


 
