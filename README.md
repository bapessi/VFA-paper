# VFA-paper
To determine the fluxes of the metabolites and draw it in the svg file this are the following steps:
## Getting the stoich matrix 
1. Format the excel file as copied from the word file: ==get_stoichiometry.py==
2. Remove not feasible metabolites (only entrance or only output): ==correct_excel_files.py==
3. 

## After getting the stoichiometric matrix
1. Calculate the entrance fluxes: ==get_flux.py==
	1. Calculate the fluxes that will be used as the boundaries in the FBA calculation. We calculate the flux for BUT,ACE, GLYC and GLU. We have the use the fluxes before GAP because using was is consumed to the biomass is not the same.
	2. You can get the y value from ==flux_t.py==
2. With the entrance fluxes, calculate FBA: ==FBA.py==
	1. There is one problem with the main_metabolism which is the boundaries of the reactions. The way to determine the correct one is by manipulating the boundary of the reactions until we reacth the point where the input substrates are the limiting substrates. In a perfect situation the flux of biomass should be the same as what we expect from the reduced model, but I couldn't get this.
	2. Calculate the input fluxes for the three subnetworks
4. Order the fluxes for the image position: ==get_flux_for_chart.py== :
	1. First it is necessary to determine a standard width and flux in the base FBA file (preferenceably that for the max flux the width is not too big).
5. Copy the values to the svg file with Inkscape.

## Observations

- The fluxes are equal to the value of $\alpha$ that means the kinetic vector before multiplying it with the stoichiometric matrix.
- ==efm_chart.py== calculate the EFMS of a file
- ==metabolite_dependancies.py== Determine if it is possible to determine which metabolites are dependent from previous one, zero being the input metabolites.
- 
