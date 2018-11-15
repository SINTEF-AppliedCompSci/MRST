function lsolve = selectLinearSolverAD(model, varargin)
    opt = struct('useAMGCL',            true, ...
                 'useAGMG',             true,...
                 'useILU',              true,...
                 'useSYMRCMOrdering',   true,...
                 'useCPR',              true, ...
                 'useAMGCLCPR',         true, ...
                 'tolerance',           1e-4);
    [opt, solver_arg] = merge_options(opt, varargin{:});
    solver_arg = ['tolerance', opt.tolerance, solver_arg];
    lsolve = BackslashSolverAD();
    if opt.useAMGCL
        opt.useAMGCL = checkAMGCL();
    end
    ncomp = getComponentCount(model);
    ndof = ncomp*model.G.cells.num;
    if ndof  < 10000
        % We do not need a custom linear solver
        return
    end
    isDiagonal = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
    
    if opt.useAMGCLCPR && opt.useAMGCL
        % AMGCL CPR
        mrstModule add linearsolvers
        lsolve = AMGCL_CPRSolverAD('maxIterations', 50,...
                                   'block_size', ncomp,...
                                   solver_arg{:});
        setSolverOrderingReduction(model, lsolve, ncomp, opt);
    elseif opt.useAMGCL
        % AMGCL as a preconditioner Krylov solver
        mrstModule add linearsolvers
        lsolve = AMGCLSolverAD('preconditioner', 'relaxation',...
                               'relaxation', 'ilu0', ...
                               'maxIterations', 50, solver_arg{:});
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
    [A, b] = getTestSystem();
    try
        callAMGCL(A, b);
        ok = true;
    catch
        disp('AMGCL test failed. May not be compiled.');
        ok = false;
    end
end

function ncomp = getComponentCount(model)
    names = model.getComponentNames();
    if isa(model, 'ThreePhaseBlackOilModel')
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