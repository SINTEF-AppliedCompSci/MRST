% AD-MICP
%
% Files
%   examples/runMICP1DCase - 1D MICP injection case
%   examples/runMICP3DCase - 3D MICP injection case and CO2 assessment including a leakage path
%   models/MICPModel - MICP type model
%   utility/checkCloggingMICP - Check if clogging as been reached in any cell
%   utility/equationsMICP - Assemble the linearized equations for the MICP system
%   utility/getDispersionAnddpWMICP - Compute the disperison coefficients and dpW on the cell centers
%   utility/getFluxAndPropsMICP - Compute the water and component fluxes and water mobility
%   utility/getPlotAfterStepCO2 - Dynamic plotting in 'simulateScheduleAD' for the CO2 assessment simulations.
%   utility/getPlotAfterStepMICP - Dynamic plotting in 'simulateScheduleAD' for the micp simulations.
%   utility/mrsttovtk - Write one or more data sets in vtk format for visualization in ParaView.

%{
Copyright 2021, NORCE Norwegian Research Centre AS, Computational
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
