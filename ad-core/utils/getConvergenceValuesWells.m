function [values, names, tolerances, evaluated] = getConvergenceValuesWells(model, problem)
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

    isperf = problem.indexOfType('perf');
    iswell = problem.indexOfType('well');
    isseg  = problem.indexOfType('seg');
    isnode = problem.indexOfType('node');
    isalpha  = problem.indexOfType('alpha');
    
    evaluated = (isperf | iswell | isseg | isnode | isalpha);

    values = cellfun(@(x) norm(value(x), inf), problem.equations(evaluated));
    
    tmp = find(evaluated);
    isperf = isperf(tmp);
    iswell = iswell(tmp);
    isseg  = isseg(tmp);
    isnode = isnode(tmp);
    isalpha  = isalpha(tmp);
    
    tolerances = zeros(size(tmp));
    tolerances(isperf)  = model.toleranceWellRate;
    tolerances(iswell)  = model.toleranceWellRate;
    tolerances(isseg)   = model.toleranceWellMS;
    tolerances(isnode)  = model.toleranceWellMS;
    tolerances(isalpha) =  model.toleranceWellMS;
    
    names = strcat(problem.equationNames(evaluated), ' (', problem.types(evaluated), ')');
end
