function lsolve = selectLinearSolverAD(model, varargin)
    opt = struct('useAMGCL', true, 'useAGMG', true,...
                 'useCPR', true, 'tolerance', 1e-4);
    [opt, solver_arg] = merge_options(opt, varargin{:});
    solver_arg = ['tolerance', opt.tolerance, solver_arg];
    lsolve = BackslashSolverAD();
    if opt.useAMGCL
        opt.useAMGCL = checkAMGCL();
    end
    ncomp = getComponentCount(model);
    if ncomp * model.G.cells.num < 10000
        % We do not need a linear solver
        return
    end
    
    if opt.useAMGCL
        mrstModule add linearsolvers
        if opt.useCPR
            lsolve = AMGCL_CPRSolverAD('maxIterations', 50, solver_arg{:});
        else
            lsolve = AMGCLSolverAD('preconditioner', 'relaxation',...
                                   'relaxation', 'ilu0', ...
                                   'maxIterations', 50, solver_arg{:});
            ndof = ncomp*model.G.cells.num;
            ordering = getCellMajorReordering(G.cells.num, ncomp, ndof);
            lsolve.variableOrdering = ordering;
            lsolve.equationOrdering = ordering;

        end
    elseif opt.useCPR
        if opt.useAGMG && checkAGMG()
            psolve = AGMGSolverAD('tolerance', 1e-3, 'iterations', 25);
        end
        lsolve = CPRSolverAD('ellipticSolver', psolve, solver_arg{:});
    end
    
    if ~isa(lsolve, 'BackslashSolverAD')
        if isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend')
            lsolve.reduceToCell = false;
            lsolve.amgcl_setup.active_rows = model.G.cells.num*ncomp;
        end
    end
end

function ok = checkAMGCL()
    A = speye(5);
    b = ones(5, 1);
    try
        callAMGCL(A, b);
        ok = true;
    catch
        disp('Unable to call AMGCL - may not be compiled.');
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

end

function ok = checkAGMG()
    mods = mrstPath();
    if any(strcmpi(mods, 'agmg'))
        mrstModule add agmg
        A = speye(5);
        b = ones(5, 1);
        try
            agmg(A, b);
            ok = true;
        catch
            disp('Unable to call AGMG - not installed.');
            ok = false;
        end
    else
        ok = false;
    end
end