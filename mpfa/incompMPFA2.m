function state = incompMPFA2(G, mpfastruct, varargin)
%Solve incompressible flow problem (fluxes/pressures) using MPFA-O method.
%
% SYNOPSIS:
%   state = incompMPFA(state, G, T, fluid)
%   state = incompMPFA(state, G, T, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell pressures at the next time step in a
%   sequential splitting scheme for the reservoir simulation problem
%   defined by Darcy's law and a given set of external influences (wells,
%   sources, and boundary conditions).
%
%   This function uses a multi-point flux approximation (MPFA) method with
%   minimal memory consumption within the constraints of operating on a
%   fully unstructured polyhedral grid structure.
%
% REQUIRED PARAMETERS:
%
%   G,    - Grid and half-transmissibilities as computed by the function
%            'computeMultiPointTrans'.
%
%  mpfastruct - Computed by computeMultiPointTrans2
%
% OPTIONAL PARAMETERS:
%   wells  - Well structure as defined by functions 'addWell' and
%            'assembleWellSystem'.  May be empty (i.e., W = struct([]))
%            which is interpreted as a model without any wells.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions to
%            the reservoir flow.  May be empty (i.e., bc = struct([])) which
%            is interpreted as all external no-flow (homogeneous Neumann)
%            conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = struct([])) which is
%            interpreted as a reservoir model without explicit sources.
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

%
% RETURNS:

%
% SEE ALSO:
%   `computeMultiPointTrans2`, `addBC`, `addSource`, `addWell`

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('bc', [], 'src', [], 'wells', [],...
                 'LinSolve', @mldivide,...
                 'Verbose', mrstVerbose); 
    opt = merge_options(opt, varargin{:}); 


    if all([isempty(opt.bc), isempty(opt.src), isempty(opt.wells)])
        warning(id('DrivingForce:Missing'),...
                ['No external driving forces present in model -- ',...
                 'state remains unchanged.\n']); 
    end
    

    
    is_well_posed = false; % changed to true if pressure is set through well
                           % or boundary conditions.
    nc = G.cells.num; 
    
    div  = mpfastruct.div;
    A    = mpfastruct.A;
    F    = mpfastruct.F;
    tbls = mpfastruct.tbls;
    
    rhs = zeros(size(A, 1), 1);
    
    if ~isempty(opt.bc)
        
        bc = opt.bc;
        bctbl.faces = bc.face;
        bctbl.num   = numel(bctbl.faces);
        
        is_flux = strcmpi('flux', bc.type);
        if any(is_flux)
            facenodetbl = tbls.facenodetbl;
            bcfluxtbl.faces = bc.faces(is_flux);
            bcfluxtbl.num = numel(bcfluxtbl.faces);
            fluxvals = bc.value(is_flux);
            map = setupTableMapping(bcfluxtbl, facenodetbl, {'faces'});
            fluxterm = map*fluxval;
            
            rhs = rhs - fluxterm;
        end
    
        is_press = strcmpi('pressure', bc.type);
        if any(is_press)
            is_well_posed = true;
            factor = A(1, 1); 
            assert(factor > 0)
            facenodeexttbl = tbls.facenodeexttbl;
            bcpresstbl.faces = bc.face(is_press);
            bcpresstbl.num = numel(bcpresstbl.faces);
            pressvals = bc.value(is_press);
            map = setupTableMapping(bcpresstbl, facenodeexttbl, {'faces'});
            pressvals = map*pressvals;
            [ind, ~] = find(map);
            nc = G.cells.num;
            e_ind = nc + ind;
            A(e_ind, :) = 0;
            A(e_ind, e_ind) = factor*speye(numel(ind));
            next = facenodeexttbl.num;
            rhs(nc + (1 : next)) = rhs(nc + (1 : next)) + factor*pressvals;
        end
    end
    nnp = length(rhs); 

    % Add well equations
    if ~isempty(opt.wells)
        nw  = length(opt.wells); 
        rhs = [rhs; zeros(nw, 1)]; 
        
        C = cell(nnp, 1); 
        D = zeros(nnp, 1); 
        W = opt.wells; 
        d = zeros(G.cells.num, 1); 
        
        for k = 1 : nw
            wc = W(k).cells; 
            nwc = numel(wc); 
            w = k + nnp; 
            wi = W(k).WI; 

            d(wc) = d(wc) + wi;
            if strcmpi(W(k).type, 'bhp')
                is_well_posed = true;
                bhp = W(k).val;
                % wimax = max(wi);
                wimax = 1;
                rhs(w)  = rhs(w) + wimax*sum(wi)*bhp; 
                rhs(wc) = rhs(wc) + wi.*bhp; 
                C{k}    = wimax*sparse(ones(nwc, 1), wc, wi, 1, nnp);
                B{k}    = sparse(nnp, 1);
                D(k)    = wimax; 
            elseif strcmpi(W(k).type, 'rate')
                rate   = W(k).val;
                rhs(w) = rhs(w) - rate;
                B{k}   = - sparse(wc, ones(nwc, 1), wi, nnp, 1);
                C{k}   = sparse(ones(nwc, 1), wc, wi, 1, nnp);
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

    end
        

    if ~is_well_posed
        A(1) = 2*A(1); 
    end

    x = opt.LinSolve(A, rhs); 

    % --------------------------------------------------------------------- 
    dispif(opt.Verbose, 'Computing fluxes, face pressures etc...\t\t'); 
    pressure    = x(1 : nc); 
    e_pressure  = x(1 : nnp); % include pressure at boundary
    bc_pressure = x((nc + 1) : nnp); % pressure at boundary
    if ~isempty(opt.wells)
        wellvars = x((nnp + 1) : end);
    end
    
    flux = F*e_pressure;

    state.pressure = pressure;
    state.bc_pressure = bc_pressure;
    state.flux = flux;
    
    % return well values.
    if ~isempty(opt.wells)
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
    end
end

% -------------------------------------------------------------------------- 
% Helpers follow.
% -------------------------------------------------------------------------- 

function s = id(s)
    s = ['incompMPFA:', s]; 
end


