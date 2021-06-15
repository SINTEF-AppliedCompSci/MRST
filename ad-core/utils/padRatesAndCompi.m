function [q_new, W, isActive] = padRatesAndCompi(q_s, W, model)
% Pad one/two/threephase values with zeros corresponding to missing phases.
%
% SYNOPSIS:
%   [q, W] = padRatesAndCompi(q, W, model);
%
% DESCRIPTION:
%   This function adds padding zeros to convert rates and wells for a
%   one/two phase model to make it appear as a three phase model with zero
%   rates for the missing phases.
%
% REQUIRED PARAMETERS:
%   q_s          - Cell array of fluxes corresponding to the number of
%                  active phases in the model.
%
%   W            - Wells compatible with the current model.
%
%   model        - Model with one or more active phases, consistent with
%                  q_s and W.
%
% RETURNS:
%
%   q_s          - 1 by 3 cell array with zero values added for missing
%                  phases.
%
%   W            - Three phase wells (.compi contains three fields, again
%                   with zeros where phases are missing in the original
%                   model).
%
%   isActive     - Indicators for which phases are present.
%
% SEE ALSO:
%   PhysicalModel, WellModel

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

    isActive = model.getActivePhases();
    dummy = zeros(size(double(q_s{1})));
    
    q_new = cell(1, 3);
    
    [q_new{:}] = deal(dummy);
    q_new(isActive) = q_s;
    
    for i = 1:numel(W)
        W(i).compi = padByIndex(W(i).compi, isActive);
    end
end

function tmp = padByIndex(val, index)
    tmp = zeros(1, 3);
    tmp(index) = val;
end
