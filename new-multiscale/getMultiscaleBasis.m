function [basis, CG] = getMultiscaleBasis(CG, A, varargin)
    G = CG.parent;
    Nc = CG.cells.num;
    Nf = G.cells.num;
    
    opt = struct('type',             'rsb', ...
                 'useMex',           false, ...
                 'regularizeSys',    true, ...
                 'mexGrid',          [], ...
                 'iterations',       ceil(50*(Nf/Nc).^(1/G.griddim)), ...
                 'tolerance',        5e-3, ...
                 'implicitDual',    false, ...
                 'useControlVolume', true);

    opt = merge_options(opt, varargin{:});

    if opt.regularizeSys
        A = (A + A')/2;
        A = A - diag(sum(A, 2));
    end


    switch lower(opt.type)
        case {'msrsb', 'rsb', 'jacobi', 'smoothed', 'jacobi-mex'}
            assert(exist('mbasisSmoothed', 'file') > 0, 'MsRSB basis functions not available');
            if ~isfield(CG.cells, 'interaction')
                CG = storeInteractionRegion(CG);
            end
            if opt.useMex || strcmpi(opt.type, 'jacobi-mex')
                if isempty(opt.mexGrid)
                    opt.mexGrid = setupGridsForMex(A, CG);
                end
                B = mbasisSmoothed(opt.mexGrid, 'maxiter', opt.iterations, ...
                    'tolerance', opt.tolerance, 'omega', .66);
            else
                B = iteratedJacobiBasis(A, CG, 'iterations', opt.iterations,...
                    'incrementTol', opt.tolerance);
            end
        case {'mstpfa', 'tpfa'}
            assert(exist('createFaceBasis', 'file') > 0, 'MsTPFA basis functions not available');
            if ~isfield(CG.faces,'region')
                CG = partitionMSTPFA(CG);
            end

            faceb = createFaceBasis(CG, A);
            B = assembleCoarseOperatorsPartition(CG, faceb);
        case {'msfvm', 'msfv'}
            require msfvm
            assert(exist('createMSFVBasis', 'file') > 0, 'MsTPFA basis functions not available');
            assert(isfield(CG, 'dual'), 'The MsFV method requires a dual grid!');
            DG = CG.dual;
            if opt.implicitDual
                if ~isfield(DG, 'll')
                    DG.ll = DG.lineedge;
                end
                DG = getDualConstants(CG, DG);
                B = createMSFVBasis(A, DG, false);
            else
                if ~isfield(DG, 'explicit')
                    DG = makeExplicitDual(CG, DG);
                end
                b = constructLocalMSFVBasis(CG, DG, A);
                B = b.B;
            end
        otherwise
            error('Unknown basis function type')
    end

    if opt.useControlVolume
        R = controlVolumeRestriction(CG.partition);
    else
        R = B';
    end

    basis = struct('R', R, 'B', B, 'type', lower(opt.type));
end
