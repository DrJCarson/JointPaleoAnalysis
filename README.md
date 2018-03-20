# JointPaleoAnalysis
Code for the manuscript "Quantifying Age and Model Uncertainties in Paleoclimate Data and Dynamical Climate Models with a Joint Inferential Analysis"



Code:

Code for running SMC^2 on the CR14 model with joint age estimation

The main files are SMC2_CR14_Comp_<Data> (forced model) and SMC2_UCR14_Comp_<Data> (unforced model). <Data> indicates which data was used: SS01 is the forced simulation study, USS01 the unforced simulation study, ODP677 and ODP846 are benthic sediment cores from the ocean drilling program.

FindForcing.c and insol.f include functions used in calculating the orbital parameters. These can be compiled using, for example, R CMD SHLIB insol.f

SMC_CR14_NAM.c performs one step of the particle filter. This needs to be compiled with the gsl library, i.e. R CMD SHLIB SMC_CR14_NAM.c -lgsl -lgslcblas. 



Data:

CompactedSS01.txt / CompactedUSS01.txt - Data for the forced / unforced simulation study.

Data for the ODP677 core dated as part of the LR04 stack can be obtained from http://www.lorraine-lisiecki.com/stack.html, and and as part of the H07 stack from http://www.people.fas.harvard.edu/~phuybers/Progression/. Please cite the original sources (Shackleton et al. 1990, and Mix et al. 1995), and relevant age model paper (Lisiecki and Raymo 2005, or Huybers 2007) when using this data.


Contact:

Jake Carson
Department of Statistics
University of Warwick
Coventry
CV4 7AL
UK

jake.carson@warwick.ac.uk

