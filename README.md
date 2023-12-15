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

### Coexistence analysis

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

### Trophic cascade simulation

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

### Raw Data

### Intermediate data

### Processed data

### Figures
