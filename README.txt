Trident-Generator project is created to calculate neutrino trident process cross-section calculation, make differential distribution of outgoing particles and create MC event simulation dataset.

This project uses the Equivalent Photon Approximation (EPA) method where photon virtuality is ignored to make the calculation simple. This approximation method affects the total cross-section calculation and overestimates the values as large as 200%. However differential cross-sections (i.e. kinematics of outgoing particles) are not affected by the approximation, thus can be used as a reliable tool for MC production and use in various neutrino experiment environments.

This project is structured in the following way:

Master directory contains several sub-directories:

1. CalcHEP_generator: 


2. CrossSection:


3. Distributions:


4. LookupTable:


5. python_modules:
