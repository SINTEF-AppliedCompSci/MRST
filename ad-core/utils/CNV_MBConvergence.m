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
    Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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

    fluid = model.fluid;
    state = problem.state;
    nc = model.G.cells.num;
    pv = model.operators.pv;
    pvsum = sum(pv);

    % Grab tolerances
    tol_mb = model.toleranceMB;
    tol_cnv = model.toleranceCNV;

    evaluated = false(1, numel(problem));
    
    subW = problem.indexOfEquationName('water');
    subO = problem.indexOfEquationName('oil');
    subG = problem.indexOfEquationName('gas');
    
    active = [false, false, false];
    if any(subW),
        assert(model.water);
        BW = 1./fluid.bW(state.pressure);
        RW = double(problem.equations{subW});
        BW_avg = sum(BW)/nc;
        CNVW = BW_avg*problem.dt*max(abs(RW)./pv);
        
        evaluated(subW) = true;
        active(1) = true;
    else
        BW_avg = 0;
        CNVW   = 0;
        RW     = 0;
    end

    % OIL
    if any(subO),
        assert(model.oil);
        if isprop(model, 'disgas') && model.disgas
            % If we have liveoil, BO is defined not only by pressure, but
            % also by oil solved in gas which must be calculated in cells
            % where gas saturation is > 0.
            BO = 1./fluid.bO(state.pressure, state.rs, state.s(:, 3)>0);
        else
            BO = 1./fluid.bO(state.pressure);
        end
        
        RO = double(problem.equations{subO});
        BO_avg = sum(BO)/nc;
        CNVO = BO_avg*problem.dt*max(abs(RO)./pv);
        
        evaluated(subO) = true;
        active(2) = true;
    else
        BO_avg = 0;
        CNVO   = 0;
        RO     = 0;
    end

    % GAS
    if any(subG),
        assert(model.gas);
        if isprop(model, 'vapoil') && model.vapoil
            BG = 1./fluid.bG(state.pressure, state.rv, state.s(:,2)>0); % need to fix index...
        else
            BG = 1./fluid.bG(state.pressure);
        end
        RG = double(problem.equations{subG});
        BG_avg = sum(BG)/nc;
        CNVG = BG_avg*problem.dt*max(abs(RG)./pv);
        
        evaluated(subG) = true;
        active(3) = true;
    else
        BG_avg = 0;
        CNVG   = 0;
        RG     = 0;
    end

    % Check if material balance for each phase fullfills residual
    % convergence criterion
    MB = problem.dt*abs([BW_avg*sum(RW), BO_avg*sum(RO), BG_avg*sum(RG)])/pvsum;
    converged_MB  = MB < tol_mb;

    % Check maximum normalized residuals (maximum mass error)
    CNV = [CNVW CNVO CNVG] ;
    converged_CNV = CNV <= tol_cnv;

    converged = converged_MB(active) & converged_CNV(active);
    values = [CNV(active), MB(active)];
    cnv_names = {'CNV_W', 'CNV_O', 'CNV_G'};
    mb_names = {'MB_W', 'MB_O', 'MB_G'};
    names = horzcat(cnv_names(active), mb_names(active));
end
