function [state, report] = incompMultiscale(state, CG, T, fluid, basis, varargin)
%Solve incompressible TPFA flow problem using approximate multiscale method
%
% SYNOPSIS:
%   [state, report] = incompMultiscale(state, CG, T, fluid, basis)
%
% DESCRIPTION:
%   This function supports the same types of problems as incompTPFA.
%   Requires precomputed basis and coarsegrids. Supports flux
%   reconstruction and iterative improvement of solution.
%
% REQUIRED PARAMETERS:
%   state  - See incompTPFA.
%
%   CG     - Coarsegrid used to generate basis functions
%
%   T      - See incompTPFA.
%
%   fluid  - See incompTPFA.
%
%   basis  - Basis functions from getMultiscaleBasis.
%
% OPTIONAL PARAMETERS:
%  any         - Additional options are passed onto `incompTPFA`.
%
%  getSmoother - Smoother function from getSmootherFunction. Required if
%                iterations is set to anything larger than zero.
%
%  iterations  - Number of multiscale iterations. Requires getSmoother
%                option. Default option is zero, which means that the
%                multiscale system will be solved a single time without any
%                application of smoothers.
%
%  tolerance   - Tolerance for iterative solver. Interpreted as 
%                |A*x - b|_2 / |b|_2 <= tolerance for convergence.
%  
%  useGMRES    - If enabled, GMRES will be used to accelerate the
%                iterations.
% 
%  LinSolve    - Linear solver function handle for coarse scale system.
%                Typically only required if the coarse grid has several
%                thousand blocks.
%
%  reconstruct - If enabled, the solver will reconstruct a divergence-free
%                velocity field suitable for transport. This is not
%                required if the number of iterations is large enough that
%                the pressure is completely resolved. Note that this
%                requires that the restriction operator uses a
%                control-volume form.
% RETURNS:
%
%  state       - Solved reservoir state.
%
%  report      - Report from linear iterative solver.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
    opt = struct('getSmoother',    [], ...
                 'iterations',  0,...
                 'tolerance',   1e-6, ...
                 'useGMRES',    false, ...
                 'eliminateWells', false, ...
                 'p0',          [], ...
                 'LinSolve',    @mldivide, ...
                 'reconstruct', true);
    [opt, incompOpt] = merge_options(opt, varargin{:});
    
    G = CG.parent;
    nc = G.cells.num;
    
    [A, q] = getSystemIncompTPFA(state, G, T, fluid, incompOpt{:});    
    if numel(q) > nc
        if opt.eliminateWells
            % Solve reduced system for wells
            [A, q, A_ww, A_wp, q_w] = eliminateWellEquations(A, q, nc);

            recover = @(p) recoverWellSolution(A_ww, A_wp, q_w, p);
        else
            % Add degrees of freedom corresponding to wells on coarse scale
            nw = numel(q) - nc;
            Bw = speye(nw, nw);
            rm = sparse(nc, nw);
            rw = sparse(nw, CG.cells.num);
            basis.B = [basis.B, rm; rw, Bw];
            basis.R = [basis.R, rw'; rm', Bw'];
            CG.partition = [CG.partition; (1:nw)' + max(CG.partition)];
            recover = @(p) p;
        end
    else
        recover = @(p) p;
    end
        
    [p_ms, report] = solveMultiscaleIteratively(A, q, opt.p0, basis,...
                                                             opt.getSmoother, ...
                                                             opt.tolerance,...
                                                             opt.iterations, ...
                                                             opt.LinSolve, ...
                                                             opt.useGMRES);
    
    state = setFluxes(state, CG, T, fluid, A, q, p_ms, recover, opt, incompOpt);
    
end

function state = setFluxes(state, CG, T, fluid, A, rhs, pressure, recover, opt, incompOpt)
    G = CG.parent;
    
    [N, isNNC] = getNeighbourship(G, 'topological', true);
    
    setFlux = @(p) ...
        incompTPFA(state, G, T, fluid, 'LinSolve', @(varargin) p, 'use_trans', numel(T) == size(N, 1), incompOpt{:});
    
    p_primal = recover(pressure);
    state = setFlux(p_primal);
    
    if opt.reconstruct
        sp = reconstructPressure(CG.partition, pressure, A, rhs);
        sp = recover(sp);
        
        state_o = setFlux(sp);
        flux = state_o.flux;

        flux(CG.faces.fconn) = state.flux(CG.faces.fconn);
        if any(isNNC)
            isCrossBlock = CG.partition(G.nnc.cells(:, 1)) ~= CG.partition(G.nnc.cells(:, 2));
            bndFaceNNC = isNNC;
            bndFaceNNC(isNNC) = bndFaceNNC(isNNC) & isCrossBlock;
            flux(bndFaceNNC) = state.flux(bndFaceNNC);
        end
        state.flux = flux;
        state.reconstructedPressure = sp(1:G.cells.num);
        
        tmp = struct('Wells', [], 'bc', []);
        [tmp, ~] = merge_options(tmp, incompOpt{:});
        if ~isempty(tmp.Wells)
            isBHP = arrayfun(@(x) strcmpi(x.type, 'bhp'), tmp.Wells);
            state.wellSol(isBHP) = state_o.wellSol(isBHP);
        end
        
        if ~isempty(tmp.bc)
            isP = strcmpi(tmp.bc.type, 'pressure');
            if any(isP)
                f = tmp.bc.face(isP);
                state.flux(f) = state_o.flux(f);
            end
        end
    end
end
