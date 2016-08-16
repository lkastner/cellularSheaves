# Tropical Homology - an Extension for polymake

## Purpose
In this extension for [polymake] (https://www.polymake.org) was developed to compute tropical homology. 

We defined sheaves and cosheaves for a polyhedral complex. For those we can associate 
a chain complex and determine its Betti numbers. 

This extension currently works with polymake version 3.0. It will not work with the polymake developer version.

## Details
For a polyhedral complex, we associate to each face sigma a linear space LL(sigma). 
This information is stored as CHOSEN_BASIS. 
For tau contained in sigma there are natural projection maps LL(sigma) to LL(tau), 
which are given by the property SIMPLE_BLOCKS.
Additionally each face is given an orientation, ORIENTATION.


On a polyhedral complex one can define several sheaves and cosheaves, which are 
a property of the polyhedral complex named SHEAF resp. COSHEAF.
A (co)sheaf associates to each face of the polyhedral complex some vector space 
and maps between those spaces. This information is given as CHOSEN_BASIS (the spaces 
assigned to the faces) and BLOCKS (the maps between the spaces).
Often the spaces will not be set as this information can be part of the maps.

We have some functions that construct some relevant sheaves: the Wp and Fp sheaves.

To a (co)sheaf we associate a chain complex, from which we compute its homology, the Betti numbers.


## References
[Tropical Homology - Ilia Itenberg, Ludmil Katzarkov, Grigory Mikhalkin, Ilia Zharkov] (http://arxiv.org/abs/1604.01838)

[Superforms, Tropical Cohomology and Poincaré Duality - Philipp Jell, Kristin Shaw, Jascha Smacka] (http://arxiv.org/abs/1512.07409)


