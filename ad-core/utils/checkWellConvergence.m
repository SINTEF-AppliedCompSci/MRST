function [convergence, values, evaluated] = checkWellConvergence(model, problem)
    % Compute convergence for wells.
    %
    % SYNOPSIS:
    %   [converged, values] = CNV_MBConvergence(model, problem)
    %
    % DESCRIPTION:
    %   Compute convergence for well equations. Uses the properties
    %       model.toleranceWellRate
    %       moedl.toleranceWellBHP
    %
    % REQUIRED PARAMETERS:
    %   model      - Subclass of PhysicalModel that contain equations that
    %                are of type well and perforations.
    %
    %   problem    - LinearizedProblem class instance we want to test for
    %                convergence.
    %
    %
    % RETURNS:
    %   convergence - Boolean indicating if the state used to produce the 
    %                 LinearizedProblem has converged.
    %                  
    %   values      - Residual inf of wells.
    %
    %   evaluated   - Logical array into problem.equations indicating which
    %                 residual equations we have actually checked 
    %                 convergence for.

    %{
    Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
    isperf = problem.indexOfType('perf');
    iswell = problem.indexOfType('well');
    
    evaluated = isperf | iswell;
    
    ratevals = problem.equations(isperf);
    bhpvals = problem.equations(iswell);

    values = [cellfun(@(x) norm(double(x), inf), ratevals), ...
              cellfun(@(x) norm(double(x), inf), bhpvals)];
    
    tmp = find(isperf | iswell);
    isperf = isperf(tmp);
    iswell = iswell(tmp);
    
    convergence = false(size(tmp));
    convergence(isperf) = values(isperf) < model.toleranceWellRate;
    convergence(iswell) = values(iswell) < model.toleranceWellBHP;
end
