function lsolve = selectLinearSolverAD(model, varargin)
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

    opt = struct('useAMGCL',            true, ...
                 'useAGMG',             true,...
                 'useILU',              true,...
                 'useSYMRCMOrdering',   false,...
                 'useCPR',              true, ...
                 'useAMGCLCPR',         true, ...
                 'BackslashThreshold',  10000, ...
                 'tolerance',           1e-4);
    [opt, solver_arg] = merge_options(opt, varargin{:});
    solver_arg = ['tolerance', opt.tolerance, solver_arg];
    lsolve = BackslashSolverAD();
    ncomp = getComponentCount(model);
    ndof = ncomp*model.G.cells.num;
    if ndof <= opt.BackslashThreshold
        % We do not need a custom linear solver
        return
    end
    if opt.useAMGCL
        mrstModule add linearsolvers
        opt.useAMGCL = checkAMGCL();
    end
    isDiagonal = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
    
    if opt.useAMGCLCPR && opt.useAMGCL && opt.useCPR
        % AMGCL CPR
        lsolve = AMGCL_CPRSolverAD('maxIterations', 50,...
                                   'block_size', ncomp,...
                                   'relaxation', 'ilu0', ...
                                   's_relaxation', 'ilu0', ...
                                   solver_arg{:});
        setSolverOrderingReduction(model, lsolve, ncomp, opt);
    elseif opt.useCPR && ~isDiagonal
        % MATLAB CPR + AMGCL, AGMG or backslash for elliptic part
        parg = {'tolerance', 1e-3, 'maxIterations', 25};
        if opt.useAGMG && checkAGMG()
            psolve = AGMGSolverAD(parg{:});
        elseif opt.useAMGCL && checkAMGCL()
            psolve = AMGCLSolverAD(parg{:});
        else
            psolve = BackslashSolverAD();
        end
        lsolve = CPRSolverAD('ellipticSolver', psolve, solver_arg{:});
    elseif opt.useAMGCL
        % AMGCL as a preconditioner Krylov solver
        lsolve = AMGCLSolverAD('preconditioner', 'relaxation',...
                               'relaxation', 'ilu0', ...
                               'block_size', ncomp, ...
                               'maxIterations', 50, solver_arg{:});
        setSolverOrderingReduction(model, lsolve, ncomp, opt);

    elseif opt.useILU
        % Matlab GMRES + ILU(0)
        lsolve = GMRES_ILUSolverAD('maxIterations', 100,...
                                   'reorderEquations', false, ...
                                    solver_arg{:});
        setSolverOrderingReduction(model, lsolve, ncomp, opt);
    end
end

function setSolverOrderingReduction(model, lsolve, ncomp, opt)
    G = model.G;
    ndof = ncomp*G.cells.num;
    if opt.useSYMRCMOrdering
        sym_ordering = getGridSYMRCMOrdering(G);
    else
        sym_ordering = [];
    end
    
    ordering = getCellMajorReordering(G.cells.num, ncomp, ...
        'ndof', ndof, ...
        'cell_ordering', sym_ordering);
    lsolve.variableOrdering = ordering;
    lsolve.equationOrdering = ordering;
    if isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend') || ...
       isa(lsolve, 'AMGCLSolverAD')
        lsolve.reduceToCell = false;
        lsolve.keepNumber = ndof;
    else
        lsolve.reduceToCell = true;
    end
end

function ok = checkAMGCL()
    if mrstSettings('get', 'useMEX')
        [A, b] = getTestSystem();
        try
            callAMGCL(A, b);
            ok = true;
        catch
            disp('AMGCL test failed. May not be compiled.');
            ok = false;
        end
    else
        ok = false;
    end
end

function ncomp = getComponentCount(model)
    model = model.validateModel();
    names = model.getComponentNames();
    if isa(model, 'GenericReservoirModel')
        ncomp = numel(names);
    elseif isa(model, 'ThreePhaseBlackOilModel')
        ncomp = numel(names) + sum(model.getActivePhases);
    else
        % Probably compositional
        ncomp = numel(names) + model.water;
    end
    if isprop(model, 'thermal')
        ncomp = ncomp + model.thermal;
    end

end

function ok = checkAGMG()
    mods = mrstPath();
    if any(strcmpi(mods, 'agmg'))
        mrstModule add agmg
        [A, b] = getTestSystem();
        try
            agmg(A, b);
            ok = true;
        catch
            disp('AGMG test failed. May not be installed.');
            ok = false;
        end
    else
        ok = false;
    end
end

function [A, b, result] = getTestSystem()
    A = speye(5);
    b = ones(5, 1);
    result = A\b;
end
