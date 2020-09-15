# Code-for-SALOD
Code coresponding to article "Evolutionary demographic models reveal the strength of purifying selection on susceptibility alleles to late-onset diseases". 

It allows to calculate selection coefficients for susceptibility alleles to late-onset diseases in a model incorporating 1) a ditribution of age at disease-onset, 2) different age-specific fertilities for men and women, 3) the fact that a child survival depends on maternal, paternal and grandmaternal care.

1.	The code ‘A_PARAMETERS’ should be run first.
a.	Chose in the code the parameters for Survival and Fertility
b.	Chose the ‘vec.sigma’ table code corresponding to all the sigma(y1,y2,y3) for the chosen sociocultural scenario. Because it is a little time consuming to build, it is advised to store them.
2.	Run the codes ‘B_Fn_AgeFather’, B_Fn_AgeMother’, ‘B_Fn_Salpha’ and ‘C_Fn_Unions’. 
a.	These includes all the functions corresponding to equations in supplementary text, part I and III.
3.	Run the code ‘D_Population_Dynamics_Solver’
a.	It solves the corresponding Euler-Lotka equation and return population dynamics characteristics.
4.	Function ‘E_Fn_W’ allow to calculate the selective values of carriers according to:
a.	A Ld vector of survival to the disease
b.	Whether the Mother, the Grandmother and the Father also carry the allele (TRUE or FALSE)
5.	Code ‘F_Calculations_SelectionGradients_FigMaintexte’ allows to calculate the selection coefficient presented in the figure 2A and 2B of the main texte.
6.	Code ‘F_Calculations_BigGraph’ allows to reproduce supplementary figure 1
7.	Code ‘F_Calculations_AdditionalRes’ allows to reproduce supplementary figures 2 and 3
8.	Code ‘F_Calculations_DemoRegime’ allows to reproduce supplementary figure 4 (do not forget to chose parameters for Sweden into the ‘A_parameters’ file.)

Code for calculating specific disease selection coefficients (figure 2C) is available on request.


Licence - Copyright (c) 2011-2017 GitHub Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy of this 'Code-for-SALOD' and associated documentation files, to deal in the code without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, but not sell, copies of the 'Code-for-SALOD', and to permit persons to whom the 'Code-for-SALOD' is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.


