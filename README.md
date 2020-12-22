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
found in Landa-Marbán et al. (Submitted).

The numerical examples accompanying this module are:
* micp_1Dfhs.m
* micp_2Dfhcs.m
* micp_2Dfhrs.m
* micp_2Dfls.m
* micp_2Dfvrs.m
* micp_3Dfls.m

The ad-micp module is compatible with the (freely available) Matlab 
Reservoir Simulation Toolbox (MRST) provided by SINTEF Digital, see
http://www.sintef.no/projectweb/mrst/ The ad-micp module was largely based 
on:
* Bao, K., Lie, K.-A., Møyner, O., Liu, M., 2017. Fully implicit simulation 
of polymer flooding with MRST. Comput. Geosci. 21 (5-6), 1219-1244.
https://doi.org/10.1007/s10596-017-9624-5.

## Requirements
* MRST (Tested version: 2020b)
* MATLAB (Tested version: R2020b)

## MRST dependencies
* [distmesh](http://persson.berkeley.edu/distmesh/)

## Installation
* Set the folder's name to 'ad-micp'.
* Move the folder inside the 'modules' folder in MRST.
* Run the 'startup.m' file in MRST.

## Cite
If you use ad-micp, please cite:
* Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E. Practical approaches 
to study microbially induced calcite precipitation at the field scale.
Submitted.

## Contact
David Landa-Marbán (dmar@norceresearch.no).
