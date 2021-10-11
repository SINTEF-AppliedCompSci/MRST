function state = incompMPFAlegacy(state, g, T, fluid, varargin)
% Solve incompressible flow problem (fluxes/pressures) using MPFA-O method.
%
% SYNOPSIS:
%   state = incompMPFAlegacy(state, G, T, fluid)
%   state = incompMPFAlegacy(state, G, T, fluid, 'pn1', pv1, ...)
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
%   state  - Reservoir and well solution structure either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'incompMPFAlegacy' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%   G, T   - Grid and half-transmissibilities as computed by the function
%            'computeMultiPointTransLegacy'.
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
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
%   MatrixOutput - Whether or not to return the final system matrix 'A' to
%                  the caller of function 'incompMPFAlegacy'.
%                  Logical.  Default value: MatrixOutput = FALSE.
%
%   Verbose      - Whether or not to time portions of and emit informational
%                  messages throughout the computational process.
%                  Logical.  Default value dependent on global verbose
%                  setting in function 'mrstVerbose'.
%
% RETURNS:
%   xr - Reservoir solution structure with new values for the fields:
%          - pressure     -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%          - boundaryPressure --
%                            Pressure values for all boundary interfaces in
%                            the discretised reservoir model, 'G'.
%          - flux         -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%          - A            -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
%
%   xw - Well solution structure array, one element for each well in the
%        model, with new values for the fields:
%           - flux     -- Perforation fluxes through all perforations for
%                         corresponding well.  The fluxes are interpreted
%                         as injection fluxes, meaning positive values
%                         correspond to injection into reservoir while
%                         negative values mean production/extraction out of
%                         reservoir.
%           - pressure -- Well pressure.
%
% NOTE:
%   If there are no external influences, i.e., if all of the structures
%   'wells', 'bc', and 'src' are empty and there are no effects of gravity,
%   then the input values 'xr' and 'xw' are returned unchanged and a
%   warning is printed in the command window. This warning is printed with
%   message ID
%
%           'incompMPFAlegacy:DrivingForce:Missing'
%
% EXAMPLE:
%    G   = computeGeometry(cartGrid([3,3,5]));
%    f   = initSingleFluid('mu' ,    1*centi*poise     , ...
%                          'rho', 1014*kilogram/meter^3);
%    rock.perm = rand(G.cells.num, 1)*darcy()/100;
%    bc  = pside([], G, 'LEFT', 2);
%    src = addSource([], 1, 1);
%    W   = verticalWell([], G, rock, 1, G.cartDims(2), ...
%                       (1:G.cartDims(3)), 'Type', 'rate', 'Val', 1/day(), ...
%                       'InnerProduct', 'ip_tpf');
%    W   = verticalWell(W, G, rock, G.cartDims(1),   G.cartDims(2), ...
%                       (1:G.cartDims(3)), 'Type', 'bhp', ...
%                       'Val',  1*barsa(), 'InnerProduct', 'ip_tpf');
%    T   = computeMultiPointTransLegacy(G, rock);
%
%    state = initState(G, W, 100*barsa);
%    state = incompMPFAlegacy(state, G, T, f, 'bc', bc, 'src', src, ...
%                       'wells', W, 'MatrixOutput',true);
%
%    plotCellData(G, xr.pressure);
%
% SEE ALSO:
%    `incompMPFA`, `private/incompMPFATensorAsseembly` `computeMultiPointTransLegacy`
%
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
% Written by Jostein R. Natvig, SINTEF ICT, 2009.

    opt = struct('bc'          , []       , ...
                 'src'         , []       , ...
                 'wells'       , []       , ...
                 'LinSolve'    , @mldivide,...
                 'MatrixOutput', false    ,...
                 'Verbose'     , mrstVerbose); 
    opt = merge_options(opt, varargin{:});
    opt = treatLegacyForceOptions(opt);

    g_vec = gravity(); 
    no_grav =~ (norm(g_vec) > 0); % (1 : size(g.nodes.coords, 2))) > 0); 
    if all([isempty(opt.bc),...
            isempty(opt.src),...
            isempty(opt.wells), no_grav])
        warning(id('DrivingForce:Missing'),...
                ['No external driving forces present in model -- ',...
                 'state remains unchanged.\n']); 
    end
    
    nc = g.cells.num; 
    nw = length(opt.wells); 

    [mob, totmob, omega, rho] = dynamic_quantities(state, fluid); 

    % Needed after introduction of gravity
    TT = T; 
    totFace_mob = ...
        1 ./ accumarray(g.cells.faces(:, 1), 1 ./ totmob(rldecode([1:g.cells.num]', diff(g.cells.facePos)))); 
    b = any(g.faces.neighbors == 0, 2); 
    totFace_mob(~b) = totFace_mob(~b); 
    n = size(TT.d1, 1); 
    tothface_mob_mat = sparse(1:n, 1:n, TT.d1 * totFace_mob); 
    % Tg = Tg * totmob_mat; 

    % Boundary conditions and source terms.
    % Note: Function 'computeRHS' is a modified version of function
    % 'computePressureRHS' which was originally written to support the
    % hybrid mimetic method.
    [~, gg, hh, ~, dF, ~] = computePressureRHS(g, omega, opt.bc, opt.src); 
    
    if any(~dF)
        % Account for each face being subdivided
        hh(~dF) = hh(~dF) ./ TT.counts(~dF); 
    end
    % add gravity contribution for each mpfa half face
    grav = -omega(TT.cno) .* (TT.R * reshape(g_vec(1:g.griddim), [], 1)); 
    
    b = any(g.faces.neighbors == 0, 2); 
    cf_mtrans = TT.cf_trans; 
    % Define div operaor form mpfa sides to cell values in addtion to the fluxes
    % out of boundary.  e_div = [TT.C,- TT.D(:, sb)]' * TT.Do;
    e_div = TT.e_div; 
    % Multiply fluxes with harmonic mean of mobility in order to avoid for
    % re-asssembly. For coupled reservoir simulation the value of
    % the sum of upwind mobility should be used.
    A = e_div * tothface_mob_mat * cf_mtrans; 
    dghf = TT.cf_trans_g * grav; 
    rhs_g = e_div * tothface_mob_mat * dghf; 
    
    hh_tmp = TT.d1 * hh; 
    rhs = [gg;- hh_tmp(TT.sb)]; 
    rhs = rhs + rhs_g; 
    % Dirichlet condition. If there are Dirichlet conditions, move contribution to
    % rhs.  Replace equations for the unknowns by speye( * ) * x(dF) = dC.

    factor = A(1, 1); 
    assert (factor > 0)
    if any(dF)
        dF_tmp = TT.d1(TT.sb, :) * dF; 
        ind = [false(g.cells.num, 1); dF_tmp>0]; 
        is_press = strcmpi('pressure', opt.bc.type); 
        face = opt.bc.face (is_press); 
        bcval = opt.bc.value (is_press); 
        dC_tmp = TT.d1(TT.sb, face) * bcval; 
        rhs = rhs - A(:, g.cells.num + 1:end) * dC_tmp; 
        rhs(ind) = factor * dC_tmp(dF_tmp>0); 
        A(ind, :) = 0; 
        A(:, ind) = 0; 
        A(ind, ind) = factor * speye(sum(ind)); 
    end
    nnp = length(rhs); 
    rhs = [rhs; zeros(nw, 1)]; 

    %%%%%%%%%%%%%%%%%%% 
    % add well equations
    C = cell (nnp, 1); 
    D = zeros(nnp, 1); 
    W = opt.wells; 
    d = zeros(g.cells.num, 1); 
    for k = 1 : nw
        wc = W(k).cells; 
        nwc = numel(wc); 
        w = k + nnp; 

        wi = W(k).WI .* totmob(wc); 

        dp = computeIncompWellPressureDrop(W(k), mob, rho, norm(gravity)); 
        d   (wc) = d   (wc) + wi; 
        state.wellSol(k).cdp = dp; 
        if strcmpi(W(k).type, 'bhp')
            ww = max(wi); 
            rhs (w) = rhs (w) + ww * W(k).val; 
            rhs (wc) = rhs (wc) + wi .* (W(k).val + dp); 
            C{k} = -sparse(1, nnp); 
            D(k) = ww; 

        elseif strcmpi(W(k).type, 'rate')
            rhs (w) = rhs (w) + W(k).val; 
            rhs (wc) = rhs (wc) + wi .* dp; 

            C{k} = -sparse(ones(nwc, 1), wc, wi, 1, nnp); 
            D(k) = sum(wi); 

            rhs (w) = rhs (w) - wi.' * dp; 

        else
            error('Unsupported well type.'); 
        end
    end

    C = vertcat(C{:}); 
    D = spdiags(D, 0, nw, nw); 
    A = [A, C'; C D]; 
    A = A + sparse(1:nc, 1:nc, d, size(A, 1), size(A, 2)); 

    if ~any(dF) && (isempty(W) || ~any(strcmpi({W.type }, 'bhp')))
        A(1) = 2*A(1); 
    end
    ticif(opt.Verbose); 
    p = opt.LinSolve(A, rhs); 

    tocif(opt.Verbose); 

    % --------------------------------------------------------------------- 
    dispif(opt.Verbose, 'Computing fluxes, face pressures etc...\t\t'); 
    ticif (opt.Verbose); 
    state.pressure(1 : nc) = p(1 : nc); 
    % Reconstruct face pressures and fluxes.
    b = any(g.faces.neighbors == 0, 2); 
    state.flux = TT.d1' * (tothface_mob_mat * cf_mtrans * p(1:nnp) - tothface_mob_mat * dghf); % 
    state.flux(~b) = state.flux(~b) / 2; 
    state.boundaryPressure = p(nc + 1 : nnp); 

    for k = 1 : nw
        wc = W(k).cells; 
        dp = state.wellSol(k).cdp; 
        state.wellSol(k).flux = W(k).WI .* totmob(wc) .* (p(nnp + k) + dp - p(wc)); 
        state.wellSol(k).pressure = p(nnp + k); 
    end

    if opt.MatrixOutput
        state.A = A; 
        state.rhs = rhs; 
    end

    tocif(opt.Verbose); 
end

% -------------------------------------------------------------------------- 
% Helpers follow.
% -------------------------------------------------------------------------- 

function s = id(s)
    s = ['incompMPFAlegacy:', s]; 
end

% -------------------------------------------------------------------------- 

function [mob, totmob, omega, rho] = dynamic_quantities(state, fluid)
   [rho, kr, mu] = getIncompProps(state, fluid);
    mob = bsxfun(@rdivide, kr, mu); 
    totmob = sum(mob, 2); 
    omega = sum(bsxfun(@times, mob, rho), 2) ./ totmob; 
end
