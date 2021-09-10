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
found in the two publications [A, B] in the "Cite" part of this README.

The numerical examples accompanying this module are:
* runMICP1DCase.m and
* runMICP3DCase.m.

In addition, this module includes scripts to run the studies in [A, B].

The ad-micp module is compatible with the (freely available) MATLAB
Reservoir Simulation Toolbox (MRST) provided by SINTEF Digital, see
http://www.sintef.no/projectweb/mrst/. The ad-micp module was largely based
on:
* Bao, K., Lie, K.-A., Møyner, O., Liu, M., 2017. Fully implicit simulation
of polymer flooding with MRST. Comput. Geosci. 21 (5-6), 1219-1244.
https://doi.org/10.1007/s10596-017-9624-5.

## Release
* 2021b

## Changes respect to the 2021a release
* Adding new commented numerical examples to show complete workflow for
creating, running, and analyzing 1D and 3D flow systems using this module.
* Deleting the 'addWellMICP' script from the utility folder (the scripts
have been modified accordingly to use the 'addWell' script in MRST).
* Deleting the 'co2_2Dfls' and 'co2_3Dfls' scripts from the utility folder
(these have been included inside the corresponding scripts in the
publications/paper_A folder).
* Deleting the 'simulateScheduleADMICP' script from the utility folder (the
scripts have been modified accordingly to use the 'simulateScheduleAD'
script in MRST, and a new function 'checkCloggingMicp' has been added to
check the clogging criterium that was checked before in
'simulateScheduleADMICP').
* Adding new functions 'getPlotAfterStepCO2' and 'getPlotAfterStepMICP' to
live plot the numerical results while the simulator is running.
* Moving the 'CO2Model', 'equationsCO2', 'getFluxAndPropsCO2', and
'getFluxAndPropsWater' scripts into the co2_assessment folder in the
publications/paper_A folder (these scripts are kept to run the
corresponding scripts described in publication [A], where a very simple
CO2-water model is used, while in the new added numerical example
('runMICP3DCase') the 'TwoPhaseWaterGasModel' script in the MRST co2lab
module, a more comprehensive CO2-water model, is used to asses the CO2
leakage prior and after MICP treatment).
* Improving the code presentation of all scripts in the publications folder.

## Changes respect to the first release 
* Adding the numerical example to produce the MRST figure in paper B.
* Bug fixed (adding the porosity for the attachment in the biofilm term).
* Setting the attachment rate 'ka' to 8.51e-7 1/s.
* Adding new grids which can be built using available commands in GNU Octave,
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
* MATLAB (Tested version: R2021a) or GNU Octave (Tested version: 6.3.0)

Regarding GNU Octave, all scripts have been successfully run in GNU/Linux and
Windows (tested in Ubuntu 20.04 and Windows 10 using Octave version 6.3.0).

In macOs using GNU Octave, it is possible to run some of the examples, but
the latest version (6.3.0) fails to run the scripts where distmesh is used
to build the grids (in Octave blogs it is written this issue in macs might
be fixed soon) while an older version (6.1.0) can sucessfully run distmesh
but it fails to run the CO2 assessment using the 'TwoPhaseWaterGasModel'.
Then, for running ad-micp in macOs using Octave, for now we recommend to
use the version 6.1.0, and for assessing the CO2 distribution using the
simple 'CO2Model' in the co2_assessment folder (see/run e.g.,
publications/paper_A/micp_3Dfls).

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
104, 112, 115, and 117 (this forces to use the function 'fall_back' in
line 110 instead of the function 'unique' in line 115).
* Open the scripts 'mrst-2021a/models/co2lab/utils/co2props/CO2props.m' and
'mrst-2021a/models/co2lab/utils/co2props/SampleProp2D.m' and change the
name of the function 'fields' to 'fieldnames' in lines 120 and 484
respectively.

## Cite
If you use ad-micp to write a scientific publication, please cite one of
the following papers:
* [A] Landa-Marbán, D., Tveit, S., Kumar, K., Gasda, S.E., 2021.
Practical approaches to study microbially induced calcite precipitation
at the field scale. Int. J. Greenh. Gas Control 106, 103256.
https://doi.org/10.1016/j.ijggc.2021.103256.
* [B] Landa-Marbán, D., Kumar, K., Tveit, S., Gasda, S.E.
Numerical studies of CO2 leakage remediation by micp-based plugging
technology. In: Røkke, N.A. and Knuutila, H.K. (Eds) Short Papers from the
11th International Trondheim CCS conference, ISBN: 978-82-536-1714-5,
284-290.

## Contact
David Landa-Marbán (dmar@norceresearch.no).
