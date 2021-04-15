function [pcOW, pcOG] = getCapillaryPressureBO(fluid, sW, sG)
%Compute oil-water and oil-gas capillary pressures
%
% SYNOPSIS:
%   [pcOW, pcOG] = getCapillaryPressureBO(f, sW, sG)
%
% DESCRIPTION:
%   Compute the capillary pressures relative to the oil pressure
%
% REQUIRED PARAMETERS:
%   fluid     - AD fluid model, from for example initDeckADIFluid or
%               initSimpleADIFluid.
%
%   sW        - Water saturation. Used to evaluate the water-oil capillary
%               pressure.
%
%   sG        - Gas saturation. Used to evaluate the oil-gas capillary
%               pressure.
% RETURNS:
%   pcOW - Pressure difference between the oil and water pressures. The
%          water phase pressure can be computed as
%          pW = pO - pcOW;
%
%   pcOG - Pressure difference between the oil and gas pressures. The
%          water phase pressure can be computed as
%          pG = pO + pcOG;
%
% NOTE:
%     Mind the different sign for the computation of the phase pressures
%     above!
% 

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    [pcOW, pcOG] = deal(0);
    
    % Check for capillary pressure (p_cow)
    if isfield(fluid, 'pcOW') && ~isempty(sW)
        pcOW  = fluid.pcOW(sW);
    end
    
    % Check for capillary pressure (p_cog)
    if isfield(fluid, 'pcOG') && ~isempty(sG) && nargin > 1
        pcOG  = fluid.pcOG(sG);
    end
end

