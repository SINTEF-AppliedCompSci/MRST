function numGrad = computeNumericalGradient(simRes, G, S, W, rock,     ...
                                            fluid, schedule, controls, ...
                                            objectiveFunction)
% compute numerical gradient

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


epsilon = 1e-5;

numControlWells = numel( controls.well );
numSteps        = numel( schedule);
obj       = objectiveFunction(G, S, W, rock, fluid, simRes);
valInit   = obj.val;

uInit    = [controls.well.values]';
uInit    = uInit(:);
dimU     = length( uInit );

% scale epsilon
epsilon = epsilon*norm(uInit);

h = waitbar(0,'Computing numerical gradient based on perturbed controls ...');
for k = 1:dimU
    e_k       = zeros(dimU, 1); e_k(k) = 1;
    uCur      = uInit + epsilon*e_k;
    controls  = updateControls(controls, uCur);
    schedule  = updateSchedule(controls, schedule);
    simRes    = runSchedule(simRes(1).resSol, G, S, W, rock, fluid, ...
                            schedule, 'VerboseLevel', 0);
    obj       = objectiveFunction(G, S, W, rock, fluid, simRes);
    values(k) = obj.val;
    waitbar(k/dimU,h)
end
close(h)

numGrad = reshape((values' - valInit)/epsilon, numControlWells, numSteps);

function controls = updateControls(controls, u)
% update controls - strucure based on controll vector u
numControlWells = numel(controls.well);
numControls     = length(u);

U = reshape(u, numControlWells, numControls/numControlWells)';
for k = 1: numControlWells
    controls.well(k).values = U(:, k);
end
