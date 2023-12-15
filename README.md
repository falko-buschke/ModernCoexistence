# The effect of plant resource richness on coexistence of a mammal herbivore community

## Description

This repository contains all the code and data needed to replicate a modern coexistence analysis for mammal herbiveres from three African protected areas.

The information was correct as of 15 December 2023. For any queries, contact Falko Buschke `falko.buschke@gmail.com`

## Study Area
This study relies on previously published data from three protected areas in Africa:

1. Gorongosa National Park, Mozambique
2. Serengeti National Park, Tanzania
3. Laikipia National Park, Kenya

<img src="https://github.com/falko-buschke/ModernCoexistence/blob/main/Study%20Area%20Map.png" alt="Study Area" width="500"/>

## Repository Stucture
The repository is made up of two sets of R-scripts and four sub-directories. 

The first set of R-scripts is for the multispecies coexistence analysis and the second set of scripts simulates the incremental removal of plant species to assess how tropic cascades affect herbivore coexistence. 

The four sub-directories contain (1) raw input data for analyses, (2) intermediate data for validation purposes, (3) processed output data from simulations, and (4) scripts to replicate figures.

### Scripts: Coexistence analysis

There are three R-scripts for the coexistence analysis; one for each protected area:

* `Coexistence_Gorongosa.R`
* `Coexistence_Serengeti.R`
* `Coexistence_Laikipia.R`

The scripts identify every possible combination of species and uses the MacArthur Consumer-Resource Model to assess whether the combination of species can coexist stably. Herbivore population synamics are calculated as:

${1 \over N_i } {dN_i \over dt} = b_i (\Sigma_l u_{il} w_l R_l - m_i)$

And the dynamics of the plant species (resources) are modelled as:

${1 \over R_l } {dR_l \over dt} = r_l (K_l - R_l) - \Sigma_i u_{il} N_i$

A combination of species is considered stable if it meets both of the following conditions:

1. All species in the assemblage have positive equilibrium densities ($N_i^*$, the population density when growth rates are zero).
2. None of the other species, those not in the assemblage, are able to invade (i.e. they all have negative invasion growth rates).

If these two conditions are met for multiple combinations of species, the combination with the highest species richness is seleceted to calcualte multispecies niche ($\mathcal{N}$) and fitness ($\mathcal{F}$) differnences according to the method outline by:

* Spaak & De Laender (2020) [Intuitive and broadly applicable definitions of niche and fitness differences.](https://doi.org/10.1111/ele.13511) *Ecology Letters*, 23, 1117 - 1128.

Niche differences are calculated as the ratio between the the differnece between the *invasion growth rate* and the *no-niche growth rate* and the differnece between the *maximum growth rate* and the *no-niche growth rate*:

$\mathcal{N_i} = {{f_i(0,\mathbf{N^{-i,* }}) - {f_i(\Sigma_{j \neq i}c_{ij} N_j^{-i,* },\mathbf{0})} }  \over f_i(0,\mathbf{0}) - f_i(\Sigma_{j \neq i}c_{ij} N_j^{-i,* },\mathbf{0}) }$

Fitness differences are the ratio between the *no-niche growth rate* and the *maximum growth rate*:

$\mathcal{F_i} = {{{f_i(\Sigma_{j \neq i}c_{ij} N_j^{-i,* },\mathbf{0})} }  \over f_i(0,\mathbf{0}) }$

The *no-niche growth rate* is the rate at which a population would grow if all its competitors consumed exactly the same resources (i.e. no niche differences), but continues to consume the samebulk amount of food. This calcualtion require as a conversion factor, $c_{ij}$ that translates the consumption ratio of species *j* into units of species *i*, so that $c_{ij} = {1 \over c_{ji}}$. 

### Scripts: Trophic cascade simulation

There are three R-scripts to simulate the incremental removal of plant resources and evaluate the maximum richness of stable communities:

* `Cascade_Gorongosa.R`
* `Cascade_Serengeti.R`
* `Cascade_Laikipia.R`

This simulation sample randomly (without replacement) plant species and then uses the same approach described above to determine the maximum species richness of the stable community of herbivores. Coexistence is established if:

1. All species in the assemblage have positive equilibrium densities ($N_i^*$, the population denisty when growth rates are zero).
2. None of the other species, those not in the assemblage, are able to invade (i.e. they all have negative invasion growth rates).

This whole process is iterated 100 times for each level of plant resources richness. 
 
In addition to the maximum richness of each interation, the presence of individual species is tracked. The persistence probability of each species at each level of plant reource richness is estimated by dividing the number of interation in which the species can persist stably, by the total number of iterations.

Results for richness and species' persistence probabilities are written to file and save in sub-directory `Processed_data`.

### Sub-directory: Raw Data

This directory includes three `.txt` files for each of the protected areas:

1. `Gorongosa_data.txt`: data for 11 mammal herbivore species and 144 plant resource species.
2. `Serengeti_data.txt`: data for 8 mammal herbivore species and 91 plant resource species.
3. `Laikipia_data.txt`: data for 12 mammal herbivore species and 121 plant resource species.

Each dataset has the following columns:

* **Column 1** `Species`: The common name of the herbivore species
* **Column 2** `BM`: The herbivore species bodymass, in kg.
* **Column 3** `Rmax`: The maximum annual rate of population increase, based on life-history traits
* **Column 4** `M`: The natural rate of attrition (mortality), which is the inverse of the average lifespan.
* **Column 5- onwards** `Sp1... Spn`: The proprion of plant species *1* to *n* in the diet of the hebivore species (row sums equal 1).


### Sub-directory: Intermediate data

The purpose of this sub-directory is to store intermediate results, which don't affect any other analyses, but are useful for method-validation and sense checks.

The directory includes three `.txt` files containing a $S \times S$ matrix, with the  consumption rates scaled by the efficiency at which plant biomass is converted into herbivore population growth rates: $b_i \Sigma_l u_{il} u_{jl}$. These matrices are used to calculate equilibrium densities:

* `U_matrix_Gorongosa.txt`
* `U_matrix_Serengeti.txt`
* `U_matrix_Laikipia.txt`

 Next, there are three matrices in `.txt` format, which include the conversion facotrs $c_{ij}$. The factors are calcualted as $c_{ij} = \sqrt{ \Sigma_l u_{jl}^2 \over \Sigma_l u_{il}^2}$. These matrices are included in files:

* `c_conversion_Gorongosa.txt`
* `c_conversion_Serengeti.txt`
* `c_conversion_Laikipia.txt`

For comparison, the directory also includes three impage in `.png` format, which plots the conversion factors $c_{ij}$ calcualtes as above against the ratio of minimum dietary requirements based on metabolic scaling: $m_i = 0.05. M_i^{0.77}$. This is just a sense check to confirm that the code is correct because conversion factors and consumption ratios should be storngly correlated. Images saved as files:

* `C-factors_Gorogosa.png`
* `C-factors_Serengeti.png`
* `C-factors_Laikipia.png`

Example:

<img src="https://github.com/falko-buschke/ModernCoexistence/blob/main/Intermediate_data/C-factors_Laikipia.png" alt="Conversion factors" width="400"/>

### Sub-directory: Processed data

The outputs from the coexistence analyses and tropic cascade simulations are saved in this directory.

There is a `.csv` file, which includes the estimates for $\mathcal{N}$, $\mathcal{F}$, *invasion growth rate*, *no-niche growth rate*, and *maximum growth rate* for each species in the three protected areas. this file is used to reproduce the composite figures for all protected areas:

* `MCT_combined_data`

Biplots of niche and fitness differences are saved in three `.png` files:

* `Gorongosa_coexistence.png`
* `Serengeti_coexistence.png`
* `Laikipia_coexistence.png`

Example:

<img src="https://github.com/falko-buschke/ModernCoexistence/blob/main/Processed_data/Laikipia_coexistence.png" alt="CoexistenceLaikipia" width="400"/>

The similation outputs include three `.txt` files that show the maximum species richness of of a stable community for each level of plant resources richness and for each of the 100 iterations:

* `Gorongosa_SpRichIter100.txt`
* `Serengeti_SpRichIter100.txt`
* `Laikipia_SpRichIter100.txt`

Lastly, there a three `.txt` files with the **persistence probability** for each speices for each level of plant resource richness:

* `Gorongosa_SpProbIter100.txt`
* `Serengeti_SpProbIter100.txt`
* `Laikipia_SpProbIter100.txt`  

### Sub-directory: Figures

This directory include the R-scripts used to reproduce the publication-ready composite figures of all the protected areas. **Both scripts rely on input files from the `Processed data` directory.**

First, the script `Combined_Coexistence.R` is used to produce the figure `Combined Coexistence.png`.

Example:

<img src="https://github.com/falko-buschke/ModernCoexistence/blob/main/Figures/Combined Coexistence.png" alt="Coexistence" width="1200"/>

Second, the script `Cascade_plots.R` is used to produce the figure `Complete_cascade.png`.

Example:

<img src="https://github.com/falko-buschke/ModernCoexistence/blob/main/Figures/Complete_cascade.png" alt="Coexistence" width="1200"/>
