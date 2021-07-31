# Local tree methods for classification

This repository contains the R codes for the reproducibility of the results associated to the manuscript in preparation **Local tree methods for classification** by *L. Raynal, A. Cleynen and J.-M. Marin*. 

## Description

This project focuses on the study and development of local random forest methodologies for the classification task. Here, *local* means that only one test data is of interest, and these random forest methods take into consideration this unique observation to potentially improve the predictive accuracy, compared to the classic (eager) Breiman's random forest (Breiman, 2001).

Usually, the splits on the covariate space induced by the construction of a regular classification tree can be irrelevant for the predictions of certain test data, and even have a negative impact. For example, on the two-class classification problem below, the first cut would occur in X1=0.5, to separate the two green and red high density blocs. But a more relevant one for you observed data (black *) would have been in X2=0.25. Hence, the development of such local approaches aims at circumventing these irrelevant splits by incorporating the additional information provided by a specific test data.

<img src="/Images/FragIssue-1.png" alt="Fragmentation Graph," style="zoom:38%;" />

## Methods

The studied/proposed algorithms fall in four categories depending on how the test data is used to modify the usual random forest to turn it into a local one. Different R files provide the source code for the different methods when built from scratch.

1. **Modification of the random forest splitting rule**:
   - The lazy decision tree introduced by Friedman et al. (1997) is here augmented into a forest and is implemented in `LazyDRF.R`.
   - We propose to incorporate the observed information thanks to kernels giving higher weights to training data close to the test instance. We can use either one kernel per dimension of the covariate space (which is implemented in `UniDKernelRF.R`), or one multidimensional kernel (implemented in `MultiDKernelRF.R`).
2. **Weighting of the training instances based on similarity with the test data**:
   * The case specific random forest proposed by Xu et al. (2016) is used by means of the `ranger` R package that provides an implementation.
   * As suggested by Fulton et al. (1996), we can pre-select a subset of the training data thanks by nearest neighborhood, to then train an eager random forest on this smaller training set. This is directly done  inside the scripts associated to the three examples described below.
3. **Weighting of the covariates based on their local importance with the test data**:
   * Thanks to a first Breiman's random forest, we propose to extract the importance of the different covariates with regards to the test data we are interested in, to then build another one taking profit of these importance during the covariate samplings. The weight computation function is implemented in `LocalVarImpRF.R` and given to the `ranger` function.
4. **Modification of the tree aggregation scheme**:
   * Tsymbal et al. (2006) propose to weight each trees of a random forest based on their capacity to  correctly predict the training data that are similar to the observed data. These weights are then used into a weighted aggregation scheme. This strategy employs the forest proximity matrix and the margin function, it is implemented in the file `DynamicVotingWithSelectionRF.R`.
   * We propose an alternative using kernels to measure the similarity between training set and the observation, available in `KernelVotingRF.R`.

## Comparison examples

The different R files mentioned above are sourced for the study and comparisons of the methodologies on three different examples trying to put in decline the usual random forest algorithm.

1. **A four-class balanced Gaussian mixture example** (file `Ex1-GaussianExampleBalanced.R`).
2. **A four-class unbalanced Gaussian mixture example** (file `Ex2-GaussianExampleUnbalanced.R`).
3. **A three-class population genetics example** (file `Ex3-SnpExample.R`).

The two first examples data are simulated directly into their corresponding R files. The third example uses simulations generated thanks to the DIYABC (Cornuet et al., 2014) software and saved in `reftable.bin`. The `header.txt` file is also required.

### Prerequisites

The packages required to run these examples are:

`abcrf`, `discretization`, `doParallel`, `MASS`, `mvtnorm`, `ranger`, `Rcpp`.

Run the following R command to install them if needed:

```R
install.packages(c("abcrf", "discretization", "doParallel", "MASS", "mvtnorm", "ranger", "Rcpp"))
```

### Running the examples

The three examples can be reproduced by simply running the associated R scripts.

## Usage example

If one is interested in toying around with the different methods, here is a simple example usable on the `iris` data set. We display below an example for the local random forest using multidimensional kernels in the splitting rule.

#### Load the method file

```R
source("MultiDKernelRF.R")
source("LazyDRF.R")
```

#### Data formatting

Most methods require the specification of the training covariance matrix, the response values associated to it, and the observed covariate we want the prediction for, as well as additional method specific arguments. 

```R
# Training covariates (we keep the first instance as test)
x <- iris[-1,-5]
# Response values
y <- factor(iris[-1,5])
# Vector of observed covariates
obs <- iris[1,-5]
```

#### Construction and prediction

```R
# For the forest with multidimensional kernel
res1 <- forestLocalMultiDim(x, y, obs, multiMinNodeSize = 1, alpha=0.9, ntree = 20, bootstrap = TRUE, whichKernel="MultiGauss", hfixe=TRUE)
res1$prediction # displays the prediction for the observed data

# For the forest of lazy decision trees
res2 <- LazyOptionForest(x, y, obs, Nmin=5, pGain=0.9, bootstrap=TRUE, ntree=20)
res2$prediction
```

## Authors

* **Louis Raynal** - *Harvard T.H. Chan School of Public Health*
* **Alice Cleynen** - *University of Montpellier*
* **Jean-Michel Marin** - *University of Montpellier*

## Acknowledgments

Many thanks to Erwan Scornet and Gérard Biau for some insightful discussions concerning the initial development of local tree methods.

## References

- Breiman L. (2001) Random forests. *Machine Learning*, 45:5-32.
- Cornuet J.-M., Pudlo P., Veyssier J., Dehne-Garcia A., Gautier M., Leblois R., Marin J.-M. and Estoup A.  (2014) DIYABC v2.0: a software to make approximate Bayesian computation inferences about population history using single nucleotide polymorphism, DNA sequence and microsatellite data. *Bioinformatics*, 30(8):1187–1189.
- Friedman J. H., Kohavi R. and Yun Y. (1997) Lazy decision trees. *Proceedings of the 13th National Conference on AAAI*, 717-724.
- Fulton T., Kasif S., Salzberg S. and Waltz D. L. (1996) Local induction of decision trees: towards interactive data mining. *Proceedings of the Second International Conference on Knowledge Discovery and Data Mining*, *AAAI Press*, 14-19.
- Tsymbal A., Pechenizkiy M. and Cunningham P. (2006) Dynamic integration with random forests. *In Fürnkranz J., Scheffer T. and Spiliopoulou M. (eds.), Machine Learning: ECML 2006.* *Lecture Notes in Computer Science*, 801-808,  Springer Berlin, Heidelberg.
- Xu R., Nettleton D. and Nordman D. J. (2016) Case-specific random forests. *Journal of Computational and Graphical Statistics*, 25(1):49-65.
