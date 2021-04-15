function [values, names, tolerances, evaluated] = getConvergenceValuesCNV(model, problem)
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
    % Grab tolerances
    tol_mb = model.toleranceMB;
    tol_cnv = model.toleranceCNV;
    v = model.verbose;

    evaluated = false(1, numeq(problem));
    pvt = model.PVTPropertyFunctions;
    b = pvt.get(model, state, 'ShrinkageFactors');
    pv = pvt.get(model, state, 'PoreVolume');
    % Value of pore-volume and total pore-volume
    pv = value(pv);
    pvtot = sum(pv);
    if isfield(state, 'cnvscale')
        pv = pv.*state.cnvscale;
    end
    
    nph = numel(b);
    CNV = zeros(1, nph);
    MB = zeros(1, nph);
    [shortPhase, phaseNames] = model.getPhaseNames();
    dt = problem.dt;
    isMass = isa(model, 'GenericReservoirModel');
    if isMass
        rhoS = pvt.get(model, state, 'SurfaceDensity');
    end
    active = true(nph, 1);
    for ph = 1:nph
        eq_ix = problem.indexOfEquationName(phaseNames{ph});
        active(ph) = any(eq_ix);
        if active(ph)
            eq = value(problem.equations{eq_ix});
            if isMass
                eq = eq./value(rhoS{ph});
            end
            if isempty(b{ph})
                B = 1;
            else
                B = 1./value(b{ph});
            end
            B_avg = mean(B);
            % Volume error: Maximum point-wise saturation error, scaled to
            % surface volume via average b-factor for phase
            [mv, m_ix] = max(abs(eq)./pv);
            CNV(ph) = B_avg*dt*mv;
            % Total mass balance error
            MB(ph) = abs(B_avg*sum(eq))/pvtot;
            evaluated(eq_ix) = true;
            if v > 1
                fprintf('CNV_%s: %d is the worst cell at %g\n',...
                                                shortPhase(ph), m_ix, CNV(ph));
            end
        end
    end
    mb_names = arrayfun(@(x) ['MB_', x], shortPhase, 'UniformOutput', false);
    cnv_names = arrayfun(@(x) ['CNV_', x], shortPhase, 'UniformOutput', false);

    % Combine
    tolerances = [repmat(tol_cnv, 1, sum(active)), repmat(tol_mb, 1, sum(active))];
    values = [CNV(active), MB(active)];
    names = [cnv_names(active), mb_names(active)];
end
