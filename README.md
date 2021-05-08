# ad-micp: A module for microbially induced calcite precipitation (MICP)

## Description
MICP is a new and sustainable technology which utilizes biochemical 
processes to create barriers by calcium carbonate cementation; therefore, 
this technology has a potential to be used for sealing leakage zones in 
geological formations. We have implemented a mathematical model for MICP 
suitable for field-scale studies. The main mechanisms in the conceptual 
model are as follow: microbial solution is injected which results in 
microbes attaching themselves to the pore walls to form biofilm, growth 
solution is injected to stimulate the development of biofilm, cementation 
solution is injected which is used by the biofilm for production of 
calcite, and the calcite reduce the pore space which in turn decreases the 
rock permeability. Further information on the mathematical model can be 
found in [A,B].

The numerical examples accompanying this module are:
* micp_1Dfhs.m
* micp_2Dfhcs.m
* micp_2Dfhrs.m
* micp_2Dfls.m
* micp_2Dfvrs.m
* micp_3Dfls.m
* micp_mrst_opm.m

The ad-micp module is compatible with the (freely available) MATLAB 
Reservoir Simulation Toolbox (MRST) provided by SINTEF Digital, see
http://www.sintef.no/projectweb/mrst/. The ad-micp module was largely based 
on:
* Bao, K., Lie, K.-A., Møyner, O., Liu, M., 2017. Fully implicit simulation 
of polymer flooding with MRST. Comput. Geosci. 21 (5-6), 1219-1244.
https://doi.org/10.1007/s10596-017-9624-5.

## Release
* 2021a

## Changes respect to the first release 
* Adding the numerical example to produce the MRST figure in paper B.
* Bug fixed (adding the porosity for the attachment in the biofilm term).
* Setting the attachment rate 'ka' to 8.51e-7 1/s.
* Adding new grids which can be built using available commands GNU Octave,
while still using the original ones (in [A]) when ad-micp is run in MATLAB.
* Adding a new function 'mrsttovtk' to write the simulation results into 
files 'name'.pvd and 'name'-n.vtu for visualization in ParaView when 
ad-micp is run in GNU Octave, while still plotting the results in the 
original format (in [A]) when ad-micp is run in MATLAB.   
The previous changes result in slight deviations from prior numerical 
results and Figures in [A]; however, the main results and conclusions 
remain the same.

## Requirements
* MRST (Tested version: 2021a)
* MATLAB (Tested version: R2021a) or GNU Octave (Tested version: 6.1.0)

## MRST dependencies
* [distmesh](http://persson.berkeley.edu/distmesh/)

## Installation (MATLAB)
* Set the folder's name to 'ad-micp'.
* Move the folder inside the 'modules' folder in MRST.
* Run the 'startup.m' file in MRST.

## Installation (GNU Octave)
* Set the folder's name to 'ad-micp'.
* Move the folder inside the 'modules' folder in MRST.
* Run the 'startup.m' file in MRST.
* Open the script 'mrst-2021a/utils/uniqueStable.m' and comment the lines
104, 112, and 117 (this forces to use the function 'fall_back' in line 110 
instead of the function 'unique' in line 115). This is necessary for older
versions of GNU Octave (e.g., 5.2.0) while for recent version (e.g, 6.1.0)
not doing this would result in annoying warnings.  

## Cite
If you use ad-micp to write a scientific publication, please cite one of 
the following papers:
[A] Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021. 
Practical approaches to study microbially induced calcite precipitation 
at the field scale. Int. J. Greenh. Gas Control 106, 103256.
https://doi.org/10.1016/j.ijggc.2021.103256.
[B] Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E. 
Numerical studies of CO2 leakage remediation by micp-based plugging 
technology. Submitted.

## Contact
David Landa-Marbán (dmar@norceresearch.no).
