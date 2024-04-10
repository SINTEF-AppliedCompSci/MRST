function solver = setUpGeothermalSolverAD(model)
%Experimental function for setting up linear solver for geothermal systems

    solver = AMGCLSolverAD('block_size', 2);
    solver.maxIterations = 100;
    solver.tolerance = 1e-4;
    solver.setPreconditioner('relaxation');
    solver.setRelaxation('ilu0');
    solver.setSolver('gmres');
    
%     solver.amgcl_setup.block_size = 2;
    nc = model.G.cells.num;
    solver.keepNumber = 2*nc;
    solver.reduceToCell = false;
    order = [1:nc; (1:nc) + nc];
    order = order(:);
%     
%     o = getCellMajorReordering(nc, 2);
%     
    [solver.equationOrdering, solver.variableOrdering] = deal(order);
    
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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