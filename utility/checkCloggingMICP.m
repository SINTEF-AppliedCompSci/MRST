function fn = checkCloggingMICP(ok)
% Check if clogging as been reached in any cell.
% 
% This function is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/solvent/utils/getPlotAfterStepSolvent.m
%
% We refer to that function for a complete commented version of the file. 
% In this file we comment on some of the lines.
%{
Partial copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021 NORCE Norwegian Research Centre AS, Computational 
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

fn = @(model, states, reports, solver, schedule, simtime) ... 
                     afterStepFunction(model, states, reports, solver, ok);                                           
end

function [model, states, reports, solver, ok] = ...
                      afterStepFunction(model, states, reports, solver, ok)
    computed = cellfun(@(x) ~isempty(x), states);
    
    st = states(computed);
    
    current = find(computed, 1, 'last');
    
    % The simulator stops if clogging has been reached in any of the cells
    if any(model.rock.poro - st{current}.c - st{current}.b < ...
                                                          model.fluid.ptol)
        fprintf(['External function: Clogging has been reached in at ', ...
                                                      'least one cell\n']);
        ok = false; 
    end
end