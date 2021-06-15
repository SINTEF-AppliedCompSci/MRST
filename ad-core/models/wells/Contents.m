% WELLS
%
% Files
%   computeWellContributionsNew - Setup well (residual) equations and compute corresponding source terms.
%   setupWellControlEquations   - Setup well controll (residual) equations 
%   updateConnectionDP          - Explicit update of hydrostatic pressure difference between bottom hole
%   updateSwitchedControls      - Update active well controls based on well solution structure
%   updateSwitchedWellControls  - Check for violated well limits and switch controls. 
%   WellModel                   - DEPRECATED Well model for three phase flow with black oil-style fluids

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
