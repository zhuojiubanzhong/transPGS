# TLegene：Polygenic prediction for underrepresented populations through transfer learning by utilizing shared genetic similarity shared with European populations

# Introduction
The past two decades have witnessed remarkable advance of genome-wide association studies (GWASs) in identifying associated loci (mainly single-nucleotide polymorphisms [SNPs]) for traits and diseases. Because 
most human phenotypes are affected by hundreds or thousands of genetic variants, a single variant typically exerts quite weak impact compared to traditional non-genetic clinical factors and thereby only explains 
a very small proportion of phenotypic variation. However, the combination of multiple SNPs weighted by their effect sizes by creating a polygenic score (PGS) usually better reflects the genetic susceptibility to 
a disease. Such score represents an independent risk factor, which is as equally strong or much stronger than many established clinical risk factors, and has gained popularity to quantify individual's diseases risk.
It is now widely recognized that PGS together with clinical and environmental data can improve the possibility for risk stratification and early disease detection and even pave a road towards personalized intervention. 
As a result, PGS has been extensively utilized to many diseases such as cardiometabolic diseases.

However, current GWASs have been predominantly conducted in individuals of European (EUR) ancestry, with 94.6% in the EUR population and only 3.7% in the East Asian (EAS) population and 0.2% in the African (AFR) population.
Due to this underrepresentation, the performance of PGS behaves poorly in non-EUR populations, particularly in populations of AFR ancestry. For example, the PGS accuracy reduced approximately 78% across multiple traits in individuals of AFR ancestry
relative to those of EUR ancestry; similarly, the accuracy of PGS across traits was on average 40% lower in individuals of South Asian ancestry and 5% lower in individuals of EAS ancestry compared to that in those of EUR ancestry. The poor transferability 
of PGS derived from EUR ancestry data to non-EUR populations leads to great concern in health disparities . Therefore, there is an urgent need to develop novel PGS methods which can exploit data across diverse populations to better perform genetic risk prediction.

Increasing sample sizes in non-EUR GWASs for the understanding of genetic architecture underlying complex phenotypes is a necessary road for understudied populations such as EAS and AFR; but this requires plenty of expense and time. 
Alternatively, integrating existing knowledge available from EURs into non-EURs by novel approaches may be another promising strategy to improve the portability of PGS. Actually, there is a deal of evidence that significant genetic
similarity exists between the EUR and non-EUR populations at both SNP and gene levels. Such genetic similarity provides theoretical and biological support for trans-ethnic leveraging of EUR information into non-EUR studies.
Currently, there are a range of trans-ethnic statistical methods that help enhance the transferability of PGS across distinct ancestral groups; however, how to optimally integrate EUR information into non-EUR genetic research remains unknown.

Recently, transfer learning has been applied in various machine learning fields for knowledge transfer from informative auxiliary samples into target samples to improve learning ability in the target task. 
By borrowing this idea, we here propose transPGS, a novel transfer learning genetic prediction method with the P+T approach as the baseline model. Using the effect sizes of the baseline model as initial values, 
transPGS leverages trans-ethnic genetic similarity shared with the EUR population (i.e., auxiliary samples) to adapt the effect sizes in the non-EUR populations (i.e., target samples) such as AFR or EAS.
To illustrate the effectiveness of transPGS, we conduct extensive simulations and confirm that the predictive ability of transPGS is enhanced in non-EUR as the increased of trans-ethnic similarity. 
Then, we apply it to ten phenotypes with individual-level data from the UK Biobank (UKB) and the Kaiser Permanente/UCSF Genetic Epidemiology of Adult Health and Ageing Study.Overall, through simulations and real data applications,
we demonstrate that transPGS represents a flexible and effective polygenic score method, which can improve genetic prediction capability for individuals of non-EUR ancestry.

# Example
```ruby
library(data.table)
library(SKAT)
library(glmnet)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(doParallel)
sourceCpp("lmm_PXEM.cpp")
source("transPGS.R")


######T is the GWAS summary statistics for the target and auxiliary populations, including marginal effects as well as standard errors.
######G1 is the target population genotype data (matched to 1000 Genomes Project).
######G2 is the auxiliary population genotype data (matched to 1000 Genomes Project).

T <- data.frame(fread("/public/home/yiyangzhu/qianyi/lipid/function/data.txt"))
G1 <- data.frame(fread("/public/home/yiyangzhu/qianyi/lipid/function/target_geno.txt"))
G2 <- data.frame(fread("/public/home/yiyangzhu/qianyi/lipid/function/auxiliary_geno.txt"))
a1 <- transPGS(T,G1,G2)

head(a1)
  origina_beta       tl_beta
1 -0.006349755  0.0000894057
2  0.076720675  0.0008204309
3 -0.060972709 -0.0009179851
4 -0.057363847 -0.0054121625
5  0.180453372  0.0069451656
6  0.020002770  0.0009623524
        
```
  
# Cite
Yiyang Zhu<sup>$</sup>, Wenying Chen<sup>$</sup> , Wenying Chen<sup>$</sup> and Ping Zeng<sup>#</sup> (2024). Polygenic prediction for underrepresented populations through transfer learning by utilizing shared genetic similarity shared with European populations.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact Ping Zeng via zpstat@xzhmu.edu.cn.
