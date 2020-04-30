function state = incompMPFA(G, mpfastruct, W, varargin)
% Solve incompressible flow problem (fluxes/pressures) using MPFA-O method.
% Only Neumann bc and well input (src could be implemented too).
%
%
% SYNOPSIS:
%   function state = incompMPFA(G, mpfastruct, W, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G          - Grid
%   mpfastruct - Assembly structure as computed by `computeNeumannMultiPointTrans` or `blockComputeNeumannMultiPointTrans`
%   W          - Well structure as defined by functions 'addWell' and 'assembleWellSystem'
%   varargin   - see below
%
% KEYWORD ARGUMENTS:
%
%   LinSolve     - Handle to linear system solver software to which the
%                  fully assembled system of linear equations will be
%                  passed.  Assumed to support the syntax
%
%                        x = LinSolve(A, b)
%
%                  in order to solve a system Ax=b of linear equations.
%                  Default value: LinSolve = @mldivide (backslash).
%
%   Verbose      - true if verbose
%   outputFlux   - If true, add fluxes in state output
%   MatrixOutput - If true, save system matrix A in the state variable state.A.
%   
% RETURNS:
%   state - state updated with pressure values and well solution
%
% EXAMPLE:
%
% SEE ALSO:
%
%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('LinSolve', @mldivide,...
                 'Verbose', mrstVerbose, ...
                 'outputFlux', false, ...
                 'MatrixOutput', false); 
    opt = merge_options(opt, varargin{:}); 

    is_well_posed = false; % changed to true if pressure is set through well
                           % or boundary conditions.
    nc = G.cells.num; 
    
    A = mpfastruct.A;
    
    nw  = length(W); 

    rhs = sparse(nc + nw, 1); 
    C   = cell(nw, 1); 
    B   = cell(nw, 1);
    D   = sparse(nc, 1); 
    d   = sparse(nc, 1);
    
    for k = 1 : nw
        wc = W(k).cells; 
        nwc = numel(wc); 
        w = k + nc; 
        wi = W(k).WI; 

        d(wc) = d(wc) + wi;
        if strcmpi(W(k).type, 'bhp')
            is_well_posed = true;
            bhp = W(k).val;
            % wimax = max(wi);
            wimax = 1;
            rhs(w)  = rhs(w) + wimax*sum(wi)*bhp; 
            rhs(wc) = rhs(wc) + wi.*bhp; 
            C{k}    = wimax*sparse(ones(nwc, 1), wc, wi, 1, nc);
            B{k}    = sparse(nc, 1);
            D(k)    = wimax; 
        elseif strcmpi(W(k).type, 'rate')
            rate   = W(k).val;
            rhs(w) = rhs(w) - rate;
            B{k}   = - sparse(wc, ones(nwc, 1), wi, nc, 1);
            C{k}   = sparse(ones(nwc, 1), wc, wi, 1, nc);
            D(k)   = - sum(wi); 
        else
            error('Unsupported well type.'); 
        end
    end
    
    
    C = vertcat(C{:}); 
    B = horzcat(B{:}); 
    D = spdiags(D, 0, nw, nw); 
    A = [A, B; C D]; 
    A = A + sparse(1:nc, 1:nc, d, size(A, 1), size(A, 2)); 


    if ~is_well_posed
        A(1) = 2*A(1); 
    end

    rhs = full(rhs);
    x = opt.LinSolve(A, rhs); 

    % --------------------------------------------------------------------- 
    dispif(opt.Verbose, 'Computing fluxes, face pressures etc...\t\t'); 
    pressure = x(1 : nc); 
    wellvars = x((nc + 1) : end);
    

    state.pressure = pressure;
    if opt.outputFlux 
        F    = mpfastruct.F;
        flux = F*pressure;
        state.flux = flux;
    end
    
    for k = 1 : nw
        if strcmpi(W(k).type, 'bhp')
            pw = W(k).val;
        elseif strcmpi(W(k).type, 'rate')
            pw = wellvars(k);
        else
            error('Unsupported well type.'); 
        end
        wc = W(k).cells;
        wi = W(k).WI;
        state.wellSol(k).flux = wi.*(pw - pressure(wc));
        state.wellSol(k).pressure = pw;
    end
    
    if opt.MatrixOutput
        state.A = A;
    end
end


