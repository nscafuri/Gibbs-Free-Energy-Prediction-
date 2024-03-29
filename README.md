*****ABOUT SNAPE*****

During my first Postdoc I have been using VASP software to simulate the adsorption of CO2 over different monoclinic ZrO2 surfaces.
Since VASP only provided the SCF energies and the vibrational modes, I decided to write this script to 
calculate and integrate the thermodynamics parameters (H, S and G) in order to have an estimation of 
the DeltaG of adsorption, at different temperatures and pressures, by calculating the molecular partition functions.(1) 
This script is very general for any kind of system. 

Snape calculates the different contributions to S and H, considering the different motions of the system:

- GAS PHASE MOLECULES (translations, rotations and vibrations).

- SOLID SURFACES (vibrations only).   

*****INPUT*****

SNAPE requires two input files:

- vibrations -- where the vibrational modes (cm**-1) are listed in one column. 

- snape.inp  -- input file where all the input parameters are specified with the following order and units:
                
                T =  Temmperature (Kelvin)
                P =  Pressure (Pascal)
                M =  Molecular/atomic mass (a.u.,if a gas phase molecule is considered)
                sigma = parameter of molecular symmetry (adimensional).
                Escf =  SCF energy (Hartree)
                I1 = moments of Inertia listed in column (Kg/m**2). If no gas phase molecule, comment or delete this lines, then M and sigma will be neglected as well.  
                --
                In

The results are written in the snape.out file.


*****REQUIREMENTS*****

The script was tested with Python 3.6.8.


*****EXECUTION*****

python3 SNAPE.py

vibrations and snape.inp file must be in the same directory of the .py file.


***DEPENDENCIES***

- Scipy

*****EXAMPLE*****

The files are here provided for the case of H2O at 298 K and 1 atm. 


*****REFERENCES*****

(1) McQuarrie, Donald A.; "Statistical mechanics", New York: Harper & Row, 1975.
