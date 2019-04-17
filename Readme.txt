This module, developed by Wenjuan Zhang and Mohammed Al Kobaisi from the department of petroleum engineering in Khalifa University
of Science and Technology, implements nonlinear finite volume discretizations. It follows the development of paper 
W. Zhang, M. Al Kobaisi, Cell-centered nonlinear finite volume methods with improved robustness, SPE Journal (in print), 2019. 

The main files are as follows:

findHAP.m - compute the coordinates and weighting coefficients of harmonic averaging points
correctHAP.m - make corrections to the harmonic averaging points using the correction algorithm proposed in the mentioned paper.
findOSflux.m - compute the one-sided flux transmissiblity matrix
PicardNTPFA - solve the incompressible single-phase flow equation using nonlienar TPFA method with Picard iterations
PicardNMPFA - solve the incompressible single-phase flow equation using nonlinear MPFA method with Picard iterations

To reproduce the results shown in figure 8 and figure 9 in the paper, run the script file runHAP.m
To reproduce the convergence tests in section 4.1 of the paper, run the script file runConv.m
To reproduce the monotonicity test in section 4.2 of the paper, run the script file runMono.m
To reproduce the field test case in section 4.3 of the paper, run the script file runField.m
