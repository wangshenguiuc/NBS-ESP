# Network-based stratification of tumor mutations (NBS2.0)

## About
This directory contains the code and resources of the following paper:

<i>"Differentially mutated networks reveal complex mechanisms for cancer heterogeneity
". Under review. </i>

1. NBS2 is a network-based tumor stratification algorithm. It takes a somatic mutation profile as input. It then clusters tumors by integrating molecular networks with this somatic mutation profile.
2. The default network is a recently published PPI network [1]. This software also supported user-defined network.
3. Our experimental results are based on TCGA data. A collection of pre-processed TCGA somatic mutation profiles can be obtained here (/data/ME_network/)
3. Please contact Sheng Wang swang91@stanford.edu if you have issues using this software.

## Overview of the Model
We introduce NBS2 algorithm to address the contamination of frequently mutated genes in network-based patient stratification. Our solution to alleviate the noise from frequently mutated genes is only using selective key gene interactions rather than all human gene interactions in network propagation. 

<p align="center">
<img  src="figure/framework.png" width="800" height="400" > 
</p>


### Step 1. Construct a small and restrictive gene network
Briefly, we selected gene interactions according to their mutually exclusive patterns in cancer because these genes are most likely to be evolved in cancer development. Since frequently mutated genes are unlikely to exhibit mutually exclusive patterns with their neighbors, their interactions are discarded by our method, thus diminishing the influence of frequently mutated genes on other genes. 

### Step 2. Network propogation based on the constructed network
Somatic mutations for each patient are represented as a binary (1,0) states on genes, in which a ‘1’ indicates a gene for which mutation (a single-nucleotide base change or the insertion or deletion of bases) has occurred in the tumor relative to germ line. For each patient, we used network propagation to smooth the influence of each mutation over its neighborhood in the small and restrictive network.

### Step 3. Patient stratification based on the 'network-smoothed' profile
We then clustered the resulting matrix of ‘network-smoothed’ profile into a predefined number of subtypes (k=2,3,...,6) (Fig. 2b). Specifically, a compact representation of each patient was first calculated by applying singular value decomposition to ‘network-smoothed’ profile. K-means++ was then used to cluster patients into subtypes based on those compact vectors. 

### Step 4. Detecting differentially mutated subnetworks
We then detected signature subnetworks by comparing the smoothed mutation profiles of different subtypes and validated these subnetworks using various of experimental datasets. 

For further details, see Online Methods of our paper. 

## Sub-directories
  - [src] contains impelmentation of the rationale model used for the beer review data. ``main.m`` is the main function to run NBS2.0.
  - [data] contains the pre-processed TCGA data which can be used to reproduce our results. 
  - [output] contains the KM-plot and clustering assignments.

<br>

## Data
  - **Pre-processed TCGA somatic mutation profile:** We provide pre-processed TCGA somatic mutation profile of 22 cancer types at (/data/TCGA_survival/). This should be sufficient for reproducing our results.   
   - **Sample network:** We use inBioMap as the sample network in this paper. The original network can be obtained at (/data/network/InBio-Map_Entrez.sif). The constructed small and restrictive network of 22 cancer types can be accessed at (/data/ME_network/)
   
   All the TCGA data is obtained from firehose.
  
**Important Note:** all data is for research-purpose only.

<br>

## Code Usage

To run the code, you need matlab (> 2015 I used) installed. Older matlab verision may not have an implementation of K-means++. Next:
  1. Clone the NBS2 repo
  2. Run `matlab main cancer_type='OV'` to run for OV cancer. For other cancer types, please change OV to the corresponding TCGA cancer name (e.g., BLCA, HNSC)

## License
NBS2 is licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0


[1]. Li, Taibo, et al. "A scored human protein-protein interaction network to catalyze genomic interpretation." Nature methods 14.1 (2017): 61-64.

