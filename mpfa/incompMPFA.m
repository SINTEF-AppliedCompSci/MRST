function state = incompMPFA(state, G, mpfaT, fluid, varargin)
% Solve incompressible problem with mpfa transmissibilities.
%
%
% SYNOPSIS:
%   function state = incompMPFA(state, G, mpfaT, fluid, varargin)
%
% DESCRIPTION:
%
% Two versions available : 'legacy' (default) and 'tensor Assembly'.
%
% The legacy version is faster. It is limited to a mesh with grid cells where
% corners have the same number of faces as the spatial dimension (this is always
% the case in 2D but not in 3D). The tensor assembly version
% (computeMultiPointTransTensorAssembly) can handle the other cases but is
% slower (the implementation will be optimized in the future to run faster).%
%
% PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from functions 'initResSol' and 'initWellSol'
%            respectively, or the results from a previous call to function
%            'incompMPFAlegacy' and, possibly, a transport solver such as
%            function 'implicitTransport'.
%
%            not used in `tensor Assembly` version
%   G        - Grid
%   mpfaT    - transmissibilities (or assembly structure) as computed by computeMultiPointTrans
%   fluid  - Fluid object as defined by function 'initSimpleFluid' (not used for tensor assembly).
%   varargin - see below
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
%            (not supported yet in case of tensor assembly)
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
%                  the caller of function 'incompMPFA'.
%                  Logical.  Default value: MatrixOutput = FALSE.
%
%   Verbose      - Whether or not to time portions of and emit informational
%                  messages throughout the computational process.
%                  Logical.  Default value dependent on global verbose
%                  setting in function 'mrstVerbose'.
%
% RETURNS:
%   state - contains following fields:
%
%             - pressure         -- Pressure values for all cells in the
%                                   discretised reservoir model, 'G'.
%             - boundaryPressure -- Pressure values for all boundary interfaces in
%                                   the discretised reservoir model, 'G'.
%                                   (not returned in `tensor assembly` version)
%             - flux             -- Flux across global interfaces corresponding to
%                                   the rows of 'G.faces.neighbors'.
%             - A                -- System matrix.  Only returned if specifically
%                                   requested by setting option 'MatrixOutput'.
%             - wellSol          -- Well solution structure array, one element for each well in the
%                                    model, with new values for the fields:
%    
%                 - flux     -- Perforation fluxes through all perforations for
%                               corresponding well.  The fluxes are interpreted
%                               as injection fluxes, meaning positive values
%                               correspond to injection into reservoir while
%                               negative values mean production/extraction out of
%                               reservoir.
%                 - pressure -- Well pressure.
%
% EXAMPLE: `linearPressureTestMPFA`, `mpfaExample1`, `mpfaExample2`, `mpfatest`
%
% SEE ALSO: 
% `private/incompMPFAlegacy`, `private/incompMPFATensorAssembly`
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

    opt = struct('useTensorAssembly', false); 
    [opt, extra] = merge_options(opt, varargin{:});
    
    if ~opt.useTensorAssembly
        
        % use legacy implementation
        state = incompMPFAlegacy(state, G, mpfaT, fluid, varargin{:});
    
    else
        
        % use tensor assembly based implementation
        state = incompMPFATensorAssembly(G, mpfaT, extra{:});        

    end
    


end

