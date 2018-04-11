# Changelog

## Next version

+ Add this changelog (#38)
+ Handle the patched version of dalton if needed (!41)
+ Correct gamma perpendicular
+ Add excitation: `!` and `#` (!46), and extract from Dalton (archive) and Gaussian (FCHK)
+ Upgrade numpy and scipy version, remove matplotlib

## Version 0.5

"Still working for nachos"

+ Change the whole way to represent electrical properties (!40): checks that the fields sums up, fix the corresponding permutations
+ Add geometrical derivatives from dalton (!39)
+ Handle the number of electron for ECP in FCHK (!36)

## Version 0.4

"Working for nachos"

+ Fix thermochemistry for atom alone (!30)
+ Add the possibility to do full tensor differentiation (!33)
+ Add the possibility to *update* a `.dal` file (!31)
+ Improve ESML basis sets (!27 and !29)

## Version 0.3

+ Improve permutations (!16) and electrical derivatives (#26 and !18 for gamma THS)
+ Add properties for Gaussian and Dalton (!24)
+ Add the "chemistry data files" (!21)
+ Improve FCHK reading (!19)
+ (*dev*) Add a release script (#25)

## Version 0.2

+ Add GAMESS (#5) and XYZ (#3) files
+ Add general basis sets objects and Gaussian type basis sets (#21)
+ Allow to creat chemistry files from scratch (#19)
+ Add binary data file (#2)

## Version 0.1

+ Everything was [migrated](https://git.pierrebeaujean.net/pierre/qcip_tools/issues/1) from the previous project.
+ Add Dalton `.dal` files (#14)
+ Add *helpers* (#4)
+ Add the pint unit library (#12)