function solver = setUpGeothermalWellboreLinearSolver(model, varargin)
%Set up a linear solver for composite geothermal models with WellboreModel

    % Optional input arguments
    opt = struct('useCPR', true);
    opt = merge_options(opt, varargin{:});

    if ~opt.useCPR
        solver = AMGCLSolverAD();
        solver.setPreconditioner('relaxation');
        solver.setRelaxation('ilu0');
        solver.setSolver('gmres');
    else
        solver = AMGCL_CPRSolverAD();
        solver.decoupling = 'quasiimpes';
    end
    % Set maximum iterations and tolerance
    solver.maxIterations = 200;
    solver.tolerance = 1e-4;
    
    % Set block size
    bz = 1 + model.submodels.Reservoir.thermal;
    solver.amgcl_setup.block_size = bz;
    % Set equation/variable order for reservor dofs
    ncr = model.submodels.Reservoir.G.cells.num;
    rorder = repmat(1:ncr, bz, 1) + (0:bz-1)'.*ncr;
    rorder = rorder(:);
    % Set equation/variable order for well dofs
    ncw = model.submodels.Wellbore.G.cells.num;
    worder = repmat(1:ncw, bz, 1) + (0:bz-1)'.*ncw + bz*ncr;
    worder = worder(:);
    % Set equation/variable order
    order = [rorder; worder];
    [solver.equationOrdering, solver.variableOrdering] = deal(order);
    % Set cells to keep in Schur complement
    cells = model.submodels.Wellbore.G.cells.type == 0;
    if isfield(model.submodels.Wellbore.G.cells, 'hybrid')
        cells = cells & model.submodels.Wellbore.G.cells.hybrid == 0;
    end
    ncw = nnz(cells);
    solver.keepNumber = bz*(ncr + ncw);
    solver.amgcl_setup.cell_size = ncr + ncw;
    solver.reduceToCell = false;
    
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