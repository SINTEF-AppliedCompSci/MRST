# shale
This MATLAB package was developed to facilitate the compositional modeling of naturally-fractured reservoirs, such as shale gas/oil reservoirs. It implements the 3D Projection-based Embedded Discrete Fracture Model (pEDFM) presented in the following journal:

   O. M. Olorode, B. Wang, H. U. Rashid. “Three-Dimensional Projection-Based Embedded Discrete Fracture Model for Compositional Simulation of Fractured Reservoirs”. SPE Journal, 2020, Volume 25, Issue 04, 2143–2161 (https://doi.org/10.2118/201243-PA).

This shale module requires the opensource MATLAB Reservoir Simulation Toolkit (http://www.sintef.no/projectweb/mrst/). 

Easiest way to get started with this package is to clone and extract it into the modules folder in MRST, and then run the cases in the "examples" folder. These include:

   gangiPlot       - Plots the Gangi Permeability correction factor
   
   langmuirPlot    - Plots the Langmuir Isotherm
   
   stochasticFracs - Generates stochastic fractures
   
   compositional3D - Illustrates use of hfm + compositional module 
   
   sorption        - Illustrates effect of Sorption
   
   diffusion       - Illustrates effect of Diffusion
   
   gangi           - Illustrates effect of a pressure-dependent permeability
   
   eagleFord       - Application of pEDFM with fractured shale reservoirs 
   
   eagleFordEDFM   - Application of EDFM with fractured shale reservoirs

Please feel free to leave us any questions using the "Issues" tab above.
