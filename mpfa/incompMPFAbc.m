function state = incompMPFAbc(G, mpfastruct, bc, varargin)
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
%   MatrixOutput - Save system matrix A as state.A.
%

%
% RETURNS:

%
% SEE ALSO:
%   `computeMultiPointTrans2`, `addBC`, `addSource`, `addWell`

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

% -------------------------------------------------------------------------- 
% Helpers follow.
% -------------------------------------------------------------------------- 

function s = id(s)
    s = ['incompMPFA:', s]; 
end


