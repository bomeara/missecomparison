Tip rates, phylogenies and diversification: what are we estimating, and how good are the estimates? 

Pascal O Title, Dan L Rabosky

Dryad Repository README
=======================

The included data files and R scripts should make it possible to reproduce all analyses and figures. 


LIST OF FILES
-------------

******************
*** Data files ***
******************

** Trees.tar.gz
This directory contains all phylogenies, as well as BAMM analysis files, used in this study, and as described in Table 1 of the main text. File structure varies depending on data source. Most straightforward is to pair these phylogenies with the data tables described below. 
If you use these trees, please cite the appropriate source publications. 

** dataFiles/treeSummary.csv
Table containing simulation parameters and model descriptions for each phylogeny. 
	headers:
		- treeName: code name unique to each phylogeny
		- setname: code name for the simulation set the phylogeny belongs to
		- source: source publication
		- description: brief description of the simulation model
		- nTips_all: total number of tips in the phylogeny
		- processType: general phylogenetic simulation type
		- regimeMRCA: the MRCA of the rate regime
		- regimeTime: the start time of the rate regime, as measured above the root
		- regime_nTips: the number of tips in this rate regime
		- K: if a diversity-dependent model, the true carrying capacity K
		- sigma: if a simulation under the evolving rates model, the diffusion parameter sigma
		- lambda0: initial speciation rate at the start of the rate regime
		- lambda1: the rate decay parameter for speciation, 0 if under a constant rate model
		- mu0: initial extinction rate at the start of the rate regime
		- mu1: the rate decay parameter for extinction, 0 if under a constant rate model


** dataFiles/trueTipRates.csv
Table containing, for each phylogeny, the true tip rates used for simulation in the original source studies.
	headers:
		- treeName: code name unique to each phylogeny
		- setname: code name for the simulation set the phylogeny belongs to
		- tipName: numeric naming of the phylogeny tips
		- regimeID: assignment of each tip to true rate regimes
		- tipLambda: true tip speciation rate
		- tipMu: true tip extinction rate
		- tipNetDiv: true tip net diversification rate


** dataFiles/estimatedTipRates.csv
Table containing, for each phylogeny, tip rates as estimated using different tip rate metrics. Headers treeName, setname and tipname match with those in trueTipRates.csv

** dataFiles/estimatedCRBD.csv
Table containing, for each phylogeny, the estimated speciation and extinction rates from a constant-rate birth-death model. 


** dataFiles/rateErrorMetrics.csv
Table containing, for each phylogeny, accuracy of tip rates as summarized by a number of different error metrics. 

** dataFiles/regimeSummary.csv
Table containing mean tip rates by rate regime, for constant-rate multi-regime trees.

** dataFiles/regimeSummary_dd.csv
Table containing mean tip rates by rate regime, for diversity-dependent trees.

** dataFiles/regimeSummary_evolvingRates.csv
Table containing mean tip rates by rate regime, for evolving rates trees.

** dataFiles/allRatesList.rds
R object that is a merging of true and estimated tip rates tables. The object is a list where each list element is a table pertaining to a single phylogeny. Each table contains the true and estimated rates for each tip. This object is also generated in script 4. 


*****************
*** R scripts ***
*****************

** 1.simulateTreesAcrossEpsilon.R
Simulation procedure used to generate constant rate trees for Figure 1. Produces Figure S1. 

** 2.simulateEvolvingRates.R
Simulation procedure for generating the evolving rates trees.

** 3.summarizeByRegime.R
Script to calculate average tip rate of tips as grouped by true rate regime. 

** 4.calcErrorMetrics.R
Script to calculate a number of error metrics based on true and estimated tip rates.

** 5.error_relativeExtinction.R
Script to examine error as a function of relative extinction, produces Figures 1, S2, S3, S4. 

** 6.error_nRegimes.R
Script to examine error as a function of the number of the amount of rate heterogeneity in the tree. Produces Figures 3, S8, S9.

** 7.ratesComparedByRegime.R
Script to examine error in tip rates when summarized by true rate regime. Produces Figures 4, 5, S10.

** 8.tipRatesCompared.R
Script to plot true tip rates against estimated tip rates. Produces Figures 2, S5, S6, S7.

** 9.conceptualFig.R
Script to plot true and estimated rates alongside a rate shift tree and an evolving rates tree. Produces Figures 6 and 7. 

** sourceScripts/tipRateFunctions.R
R functions to calculate the various tip rates estimated for this study. 

** sourceScripts/sourceFxns.R
Miscellaneous helper functions needed for other R scripts. 

** sourceScripts/evolveRateTree_JMB_sqT.R
Evolving Rates simulation functions, as provided by Beaulieu and O'Meara 2015
https://datadryad.org/resource/doi:10.5061/dryad.33p91




