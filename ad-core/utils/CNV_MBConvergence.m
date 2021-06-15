function [converged, values, evaluated, names] = CNV_MBConvergence(model, problem)
% Compute convergence based on total mass balance and maximum residual mass balance.
%
% SYNOPSIS:
%   [converged, values, evaluated] = CNV_MBConvergence(model, problem)
%
% DESCRIPTION:
%   Compute CNV/MB type convergence similar to what is used for black
%   oil convergence in commercial simulators.
%
% REQUIRED PARAMETERS:
%   model      - Subclass of PhysicalModel. Strongly suggested to be
%                some black oil variant, as this convergence function
%                does *not* account for general residual convergence.
%
%   problem    - LinearizedProblem class instance we want to test for
%                convergence.
%
%
% RETURNS:
%   convergence - Boolean indicating if the state used to produce the 
%                 LinearizedProblem has converged.
%                  
%   values      - 1 by 6 array containing mass balance in the first
%                 three terms followed by cnv in the last three. The
%                 phase ordering is assumed to be oil, water, gas.
%                 Phases present will return a zero in their place.
%
%   evaluated   - Logical array into problem.equations indicating which
%                 residual equations we have actually checked 
%                 convergence for.
%
%   names       - Cell array of same length as values with short names
%                 for printing/debugging.

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

    state = problem.state;
    pv = model.operators.pv;
    pvsum = sum(pv);

    % Grab tolerances
    tol_mb = model.toleranceMB;
    tol_cnv = model.toleranceCNV;

    evaluated = false(1, numeq(problem));
    
    b = model.FlowPropertyFunctions.getProperty(model, state, 'ShrinkageFactors');
    dt = problem.dt;
    [ph, phase_names] = model.getPhaseNames();
    
    [CNV, MB] = deal(zeros(1, numel(ph)));
    for i = 1:numel(b)
        sub = problem.indexOfEquationName(phase_names{i});
        B = 1./b{i};
        R = double(problem.equations{sub});
        B_avg = mean(B);
        CNV(i) = B_avg*dt*max(abs(R)./pv);
        MB(i) = dt*abs(B_avg*sum(R))/pvsum;
        evaluated(sub) = true;
    end

    % Check if material balance for each phase fullfills residual
    % convergence criterion
    converged_MB  = MB <= tol_mb;

    % Check maximum normalized residuals (maximum mass error)
    converged_CNV = CNV <= tol_cnv;

    converged = [converged_CNV, converged_MB];
    values = [CNV, MB];
    cnv_names = arrayfun(@(x) ['CNV_', x], ph, 'UniformOutput', false);
    mb_names = arrayfun(@(x) ['MB_', x], ph, 'UniformOutput', false);
    names = horzcat(cnv_names, mb_names);
end
