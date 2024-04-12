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
library(harmonicmeanp)
source("TLegene.R")
source("null_model_fit.R")
source("numerical_approximation.R")
data = read.table("data.txt",head=T)
weights = read.table("weights.txt",head=T)
data<-as.matrix(data)
weights<-as.matrix(weights)
p=dim(data)[2]-3
result=TLegene(data,d=2,p,R=1,
               outcome_type="Continuous",
               weight_method= "User",
               user_weight = weights)

po=result$pvalue[3]
pa=result$pvalue[4]
pf=result$pvalue[5]
ph<-cbind(po,pa,pf)
phmp<-as.vector(c(p.hmp(ph,L=length(ph))))
$pvalue

  pvalue.TLegene-oScore   4.93183049954382e-10 

  pvalue.TLegene-aScore   2.86890693215321e-16

  pvalue.TLegene-fScore  2.94312008307861e-15

  pvalue.TLegene-HMP     7.8422646713986e-16
                             
```
  
# Cite
Shuo Zhang<sup>$</sup>, Zhou Jiang<sup>$</sup> and Ping Zeng<sup>#</sup> (2022). Incorporating genetic similarity of auxiliary samples into eGene identification under the transfer learning framework.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact Ping Zeng via zpstat@xzhmu.edu.cn.
