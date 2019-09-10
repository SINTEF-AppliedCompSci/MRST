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
    extfacenodetbl = addLocInd(extfacenodetbl, 'efnind');
    
    is_press = strcmpi('pressure', bc.type);
    if any(is_press)
        is_well_posed = true;
        bcpresstbl.faces = bc.face(is_press);
        bcpresstbl.num   = numel(bcpresstbl.faces);
        bcpresstbl.num   = numel(bcpresstbl.faces);
        bcpresstbl0 = bcpresstbl;
        pressvals = bc.value(is_press);
        [~, bcpresstbl] = setupTableMapping(extfacenodetbl, bcpresstbl, ...
                                                           {'faces'});
        map = setupTableMapping(bcpresstbl0, bcpresstbl, {'faces'});
        pressvals = map*pressvals;
        
        knownbcpind = bcpresstbl.efnind;
        unknownbcpind              = true(extfacenodetbl.num, 1);
        unknownbcpind(knownbcpind) = false;
        unknownbcpind              = find(unknownbcpind);
        
        clear bctbl
        bctbl.efnind = extfacenodetbl.efnind(unknownbcpind);
        bctbl.num = numel(bctbl.efnind);
        
        [~, bctbl] = setupTableMapping(bctbl, extfacenodetbl, {'efnind'});
        bctbl = addLocInd(bctbl, 'bcind');
        
        unknownpind = [(1 : nc)'; nc + unknownbcpind];
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
        facenodetbl = tbls.facenodetbl;
        bcfluxtbl.faces = bc.face(is_flux);
        bcfluxtbl.num   = numel(bcfluxtbl.faces);
        fluxvals = bc.value(is_flux);
        map = setupTableMapping(bcfluxtbl, bctbl, {'faces'});
        nfluxperface = diag(map'*map);
        fluxvals = (1./nfluxperface).*fluxvals;
        fluxvals = map*fluxvals;
        ind = nc + bctbl.bcind;
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
    
    unknown_bc_pressure = x((nc + 1) : nnp); % pressure at boundary
    map = setupTableMapping(bctbl, extfacenodetbl, {'faces', 'nodes'});
    unknown_bc_pressure = map*unknown_bc_pressure;
    bc_pressure = bc_pressure + unknown_bc_pressure;
    
    map = setupTableMapping(bcpresstbl, extfacenodetbl, {'faces', 'nodes'});
    known_bc_pressure = map*pressvals;
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


