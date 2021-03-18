# CSH-CH-CO2 test

This test aims to simulate a crack in a C<sub>3</sub>S paste filled with water and CO<sub>2</sub> at atmospheric saturation level. This should induce the  dissolution of Ca(OH)<sub>2</sub> to form CaCO<sub>3</sub>. The concentration of CO<sub>2</sub> in solution should remain constant during the simulation. Each particle of a different phase should represent one molecule.

The test consists in 4 steps:

1. Creating Ca(OH)<sub>2</sub> crystals at given volume fraction;
2. Fill the remaining space with amorphous C--S--H (this requires another separate test to parametrise amorphous C--S--H dissolution);
3. Crack the C--S--H - Ca(OH)<sub>2</sub>  system
4. Simulate dissolution-mineralisation in the crack

The details for each step are given below: 

## Step 1: creating Ca(OH)<sub>2</sub> crystals

This first simulation is to create two crystallites of Ca(OH)<sub>2</sub>, targetting a desired volume fraction $\eta$ =  28%. This $\eta$ assumes that the paste is pure C$_3$S and is fully hydrated, hence leading only to C--S--H gel and Ca(OH)$_2$ as hydration products. The C--S--H gel is assumed to have an internal gel porosity of approximately 34%, whereas the Ca(OH)$_2$ will be a solid crystal. The volume fractions occupied by the two minerals are then obtained from the stoichiometry and molar volumes in *Masoero et al., J. Am. Ceram. Soc. 2013*. There will also be the assumption that the uncracked paste will feature no capillary pores. The calculations to obtain $\eta$ for the two phases are given in the supporting file *Vol_fracs.txt*

The simulation starts with creating 2 small FCC nuclei of Ca(OH)$_2$ with different orientations, and setting the solution to a high concentration of Ca and OH, so that further precipitation will occur. The simulation will eventually lead to $\eta >  28\%$ for Ca(OH)$_2$, but snapshots will be saved during the simulation so that later we will be able to retrieve a configuration with $\eta \approx  28\%$.

The simulation box is set to 7 x 7 x 10 nm, which will contain approximately 5,000 particles (molecules). This is a very small system, but the simulation is meant to run reasonably fast on 1 processor only. A cluster version could scale up by a factor $\sim 5$ in all directions.

To precipitate Ca(OH)$_2$ using net rates, the solution is initially set to high supersaturation with respect to Ca(OH)$_2$ precipitation, i.e. $\beta = 390$, by setting [Ca] = 0.15 M and [OH] = 0.3 M in the *sol_start* command.

