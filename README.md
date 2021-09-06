# polytropic-fits
Python tools to fit projected polytropes (propols) to mass surface density of galaxies (as described by Sanchez Almeida et al. 2021, ApJ, in press) 
-----------

        This repository contains python scripts to fit projected
	polytropes (propols) to mass surface density profiles
	of galaxies.

	The tool is described in detail in the paper

	"Physically motivated fit to mass surface density profiles
	observed in galaxies"

	by Sanchez Almeida et al., 2021, ApJ, in press
        (it is not in arXiv yet, but an early version of
	this paper is incuded in the directory)
	
=========================


Modification history:
--------------------

	- Sept 2021. First attemp to set up all the materials needed
	   to distribute the code.

Procedure:
----------

	Download and unzip the GitHub repository 'polytropic_fit'
	

Scripts
-------

	- 'lane_emden_grid.py'
	   computes the grid to be used in the fits.
	   Just included for reference. It does not need to be run since the
	   grids described in Table 1 of the paper and used to fit are alredy
	   provided in './lane_emden_grid/' 
     
        - 'lane_emden_fit.py'
	  Shows how to carry out the fits. It has to be modified
	  for each specific application.

	  The examples selected in this script correspond to some of the
	  profiles analyzed in the paper. From Trujillo et al. 2020a,b.

	  It is prepared to be run in python interactive mode, i.e.

	  linux> ipython3
	  shell command % ipython3                
	  Python 3.9.1 (v3.9.1:1e5d33e9b9, Dec  7 2020, 12:44:01) 
	  Type 'copyright', 'credits' or 'license' for more information
	  IPython 7.19.0 -- An enhanced Interactive Python. Type '?' for help.
	  
	  In [1]: run lane_emden_fit 


Data for testing at
-------------------

        ./data/galaxies_edges_wnew_mass.profiles

Plots at
--------

	./data/plots/

Tables of propols at
--------------------

	./lane_emden_grid/

	   contains the pre-computed grid of propols
	   described in Table 1 of the paper.

Contact:
-------

To report problems or request clarifications, please contact 

Jorge Sanchez Almeida (jos@iac.es) 
