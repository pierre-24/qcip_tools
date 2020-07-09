# Changelog

+ Extract TD excitation energies from gaussian logs 
+ Remove ESML, since there is a new version (with an API)
+ Switch to pip-tools
+ Add PDB files

## Version 0.6

"I'm not dead!"

+ Add this changelog 
+ Handle the patched version of dalton if needed 
+ Correct gamma perpendicular
+ Add excitation: `!` and `#` , and extract from Dalton (archive) and Gaussian (FCHK)
+ Upgrade numpy and scipy version, remove matplotlib
+ Correct charge transfer 
+ Upgrade to the latest version of h5py, remove warning 
+ Add the latest formula for first and second hyperpolarizabilities (decomposition into spherical invariants)
+ Use pipenv 
+ Store best values in Romberg object 
+ Inverse the sign of gamma tensor components when computed via CC methods in Dalton 
+ Extract gamma tensor from Gaussian FCHK 
+ Handle symmetry 

## Version 0.5

"Still working for nachos"

+ Change the whole way to represent electrical properties : checks that the fields sums up, fix the corresponding permutations
+ Add geometrical derivatives from dalton 
+ Handle the number of electron for ECP in FCHK 

## Version 0.4

"Working for nachos"

+ Fix thermochemistry for atom alone 
+ Add the possibility to do full tensor differentiation 
+ Add the possibility to *update* a `.dal` file 
+ Improve ESML basis sets (!27 and !29)

## Version 0.3

+ Improve permutations  and electrical derivatives
+ Add properties for Gaussian and Dalton 
+ Add the "chemistry data files" 
+ Improve FCHK reading 
+ (*dev*) Add a release script 

## Version 0.2

+ Add GAMESS  and XYZ  files
+ Add general basis sets objects and Gaussian type basis sets 
+ Allow to creat chemistry files from scratch 
+ Add binary data file 

## Version 0.1

+ Add Dalton `.dal` files 
+ Add *helpers* 
+ Add the pint unit library 