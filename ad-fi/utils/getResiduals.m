function [meta, residuals] = getResiduals(meta, eqs, system, linsolver_diverged)
%Undocumented Utility Function

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

   if nargin < 4, solver_diverged = false; end

    % Store the residuals for debugging and convergence testing.
    residuals = cellfun(@(x) norm(x.val, 'inf'), eqs);

    if isempty(meta.res_history) 
        meta.res_history = zeros(system.nonlinear.maxIterations, numel(residuals));
    end

    meta.res_history(meta.iteration, :) = residuals;

    % Try a simple detection of oscillations, and relax the next iteration if
    % oscillations were detected.
    [oscillate stagnate] = detectNewtonOscillations(meta.res_history, system.cellwise, meta.iteration, system.nonlinear.relaxRelTol);
    if ~ linsolver_diverged,
        if oscillate
            meta.relax = max(meta.relax - system.nonlinear.relaxInc, system.nonlinear.relaxMax);
            dispif(mrstVerbose, ...
                  ['Oscillating behavior detected: Relaxation set ', ...
                   'to %.1g\n'], meta.relax);
        elseif stagnate && 0
            meta.relax = max(meta.relax + system.nonlinear.relaxInc, system.nonlinear.relaxMax);
            dispif(mrstVerbose, ...
                  ['Stagnating behavior detected: Relaxation set ', ...
                   'to %.1g\n'], meta.relax);
        end
    end
    meta.oscillate = oscillate;
    meta.stagnate  = stagnate;
    meta.linsolver_diverged = linsolver_diverged;
end
