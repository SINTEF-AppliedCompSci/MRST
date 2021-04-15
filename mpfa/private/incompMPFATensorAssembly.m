function state = incompMPFATensorAssembly(G, mpfastruct, varargin)
% Solve incompressible flow problem (fluxes/pressures) using MPFA-O method.
%
%
% SYNOPSIS:
%   function state = incompMPFA(G, mpfastruct, W, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G          - Grid
%   mpfastruct - Assembly structure as computed by `computeMultiPointTrans` using the `useTensorAssembly` option
%   varargin   - see below
%
% KEYWORD ARGUMENTS:
%
%   wells        - Well structure as defined by functions 'addWell'
%   bc           - Boundary condition structure as defined by function 'addBC'
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
% SEE ALSO:
% `incomMPFA`, `computeMultiPointTrans`, `private/computeMultiPointTransTensorAssembly`
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


    opt = struct('bc', [], ...
                 'wells', []); 
    [opt, extra] = merge_options(opt, varargin{:}); 
    
    if ~isempty(opt.bc)
        bc = opt.bc;
        state = incompMPFANaturalBcTA(G, mpfastruct, bc, extra{:});
    elseif ~isempty(opt.wells)
        W = opt.wells;
        state = incompMPFANeumannTA(G, mpfastruct, W, extra{:});        
    else
        error('No stimulation input, either bc or W, given');
    end
    
end

function state = incompMPFANeumannTA(G, mpfastruct, W, varargin)
% Solve incompressible flow problem (fluxes/pressures) using MPFA-O method for
% Neumann boundary conditions (no flows) and wells.
%
%
% SYNOPSIS:
%   function state = incompMPFANeumannTA(G, mpfastruct, W, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G          - Grid
%   mpfastruct - Assembly structure as computed by `computeMultiPointTrans` using the `useTensorAssembly` option
%   wells      - Well structure as defined by functions 'addWell'
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
% SEE ALSO:
%

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



function state = incompMPFANaturalBcTA(G, mpfastruct, bc, varargin)
% Solve incompressible flow problem (fluxes/pressures) using MPFA-O method for
% natural boundary conditions
%
%
% SYNOPSIS:
%   function state = incompMPFANeumannTA(G, mpfastruct, W, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G          - Grid
%   mpfastruct - Assembly structure as computed by `computeMultiPointTrans` using the `useTensorAssembly` option
%   bc         - Boundary condition structure as defined by function 'addBC'
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
% SEE ALSO:
%


    opt = struct('LinSolve'  , @mldivide  , ...
                 'Verbose'   , mrstVerbose, ...
                 'outputFlux', false, ...
                 'MatrixOutput', false);
    
    opt = merge_options(opt, varargin{:}); 

    
    is_well_posed = false; % changed to true if pressure is set through well or
                           % boundary conditions.
    nc = G.cells.num; 
    
    
    
    A    = mpfastruct.A;
    tbls = mpfastruct.tbls;
    
    rhs = zeros(size(A, 1), 1);
    
    extfacenodetbl = tbls.extfacenodetbl;
    
    is_press = strcmpi('pressure', bc.type);
    
    if any(is_press)
        
        is_well_posed = true;
        
        bcpresstbl.faces = bc.face(is_press);
        bcpresstbl = IndexArray(bcpresstbl);
        bcpresstbl0 = bcpresstbl;
        pressvals = bc.value(is_press);
        
        % bcpresstbl contains indices of external face node where pressure is given
        bcpresstbl = crossIndexArray(extfacenodetbl, bcpresstbl, {'faces'});
        map = TensorMap();
        map.fromTbl = bcpresstbl0;
        map.toTbl = bcpresstbl;
        map.mergefds = {'faces'};
        map = map.setup();
        
        pressvals = map.eval(pressvals);
        
        % setup bctbl given external face nodes where the pressure is not given. They
        % correspond to degrees of freedom in A that we have to keep
        knownbcpind = bcpresstbl.get('extfnind');
        unknownbcpind              = true(extfacenodetbl.num, 1);
        unknownbcpind(knownbcpind) = false;
        unknownbcpind              = find(unknownbcpind);
        
        extfnind = extfacenodetbl.get('extfnind');
        bctbl.extfnind = extfnind(unknownbcpind);
        bctbl = IndexArray(bctbl);
        % expand bctbl with all the indices in extfacenodetbl
        bctbl = crossIndexArray(bctbl, extfacenodetbl, {'extfnind'});
        bctbl = bctbl.addLocInd('bcind');
        
        % The degrees of freedom (unknown) that are kept in the system
        unknownpind = [(1 : nc)'; nc + unknownbcpind];
        % The external degrees of freedom that are known and removed from the system
        knownpind   = nc + knownbcpind;
        
        B   = A(unknownpind, knownpind);
        rhs = rhs(unknownpind, 1);
        
        nc = G.cells.num;
        
        A = A(unknownpind, unknownpind);            
        
        prhs = B*pressvals;
        rhs = rhs - prhs;
        
    end        
    
    is_flux = strcmpi('flux', bc.type);
    if any(is_flux)
        
        % The flux is given on faces and is divided equally on face nodes

        facenodetbl = tbls.facenodetbl;
        bcfluxtbl.faces = bc.face(is_flux);
        bcfluxtbl = IndexArray(bcfluxtbl);
        
        fluxvals = bc.value(is_flux);
        
        map = TensorMap();
        map.fromTbl = bctbl;
        map.toTbl = bcfluxtbl;
        map.mergefds = {'faces'};
        map = map.setup();
        
        nfluxperface = map.eval(ones(bcflux.num, 1));
        coef = 1./nfluxperface;
        
        prod = TensorMap();
        prod.tbl1 = bcfluxtbl;
        prod.tbl2 = bctbl;
        prod.tbl3 = bctbl;
        prod.mergefds = {'faces'};
        prod = prod.setup();
        
        fluxvals = prod.eval(coef, fluxvals);
        
        bcind = bcctbl.get('bcind');
        ind = nc + bcind;
        rhs(ind) = rhs(ind) + fluxvals;
        
    end

    nnp = length(rhs); 

    if ~is_well_posed
        A(1) = 2*A(1); 
    end

    x = opt.LinSolve(A, rhs); 

    % --------------------------------------------------------------------- 
    dispif(opt.Verbose, 'Computing fluxes, face pressures etc...\t\t'); 
    pressure = x(1 : nc); 
    state.pressure = pressure;

    bc_pressure = zeros(extfacenodetbl.num, 1);
    
    % We recover the computed values of the pressure and insert them in the vector
    % bc_pressure
    
    unknown_bc_pressure = x((nc + 1) : nnp); % pressure at boundary
    
    map = TensorMap();
    map.fromTbl = bctbl;
    map.toTbl = extfacenodetbl;
    map.mergefds = {'faces', 'nodes', 'extfnind'};
    map = map.setup();
    
    unknown_bc_pressure = map.eval(unknown_bc_pressure);
    
    bc_pressure = bc_pressure + unknown_bc_pressure;

    % We insert the given Dirichlet pressure value them in the vector bc_pressure
    
    map = TensorMap();
    map.fromTbl = bcpresstbl;
    map.toTbl = extfacenodetbl;
    map.mergefds = {'faces', 'nodes', 'extfnind'};
    map = map.setup();
    
    known_bc_pressure = map.eval(pressvals);
    bc_pressure = bc_pressure + known_bc_pressure;    
    
    if opt.outputFlux 
        F = mpfastruct.F;
        p = [pressure; bc_pressure];
        flux = F*p;
        state.flux = flux;
    end

    if opt.MatrixOutput
        state.A = A;
    end
    
end

