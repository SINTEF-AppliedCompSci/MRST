function varargout = computeMultiPointTrans(G, rock, varargin)
% Compute multi-point transmissibilities for MPFA using local flux mimetic approach
%
% SYNOPSIS:
%   function varargout = computeMultiPointTrans(G, rock, varargin)
%
% DESCRIPTION:
% 
% We follow the local flux mimetic approach as described in this reference paper:
%
%      title     = {Local flux mimetic finite difference methods},
%      author    = {Lipnikov, Konstantin and Shashkov, Mikhail and Yotov, Ivan},
%      journal   = {Numerische Mathematik},
%      volume    = {112},
%      number    = {1},
%      pages     = {115--152},
%      year      = {2009},
%      publisher = {Springer}
%
% Two versions are available : 'legacy' (default) and 'tensor Assembly' (set key option `useTensorAssembly` to true).
%
% This legacy version is faster. It is limited to mesh where all the grid cells
% have corners that have the same number of faces as the spatial dimension (this
% is always the case in 2D but not in 3D). The tensor assembly version can
% handle the other cases but is slower (We plan to optimize this implementation
% in the future to run faster).
%
% PARAMETERS:
%   G       - Grid data structure as described by grid_structure.
%
%   rock    - Rock data structure with valid field 'perm'.  The
%             permeability is assumed to be in measured in units of
%             metres squared (m^2).  Use function 'darcy' to convert from
%             (milli)darcies to m^2, e.g.,
%
%                 perm = convertFrom(perm, milli*darcy)
%
%             if the permeability is provided in units of millidarcies.
%
%             The field rock.perm may have ONE column for a scalar
%             permeability in each cell, TWO/THREE columns for a diagonal
%             permeability in each cell (in 2/3 D) and THREE/SIX columns
%             for a symmetric full tensor permeability.  In the latter
%             case, each cell gets the permeability tensor
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%   varargin - see below
%
% KEYWORD ARGUMENTS:
%   verbose           - Whether or not to emit informational messages throughout the
%                       computational process.  Default value depending on the
%                       settings of function 'mrstVerbose'.
%
%   useTensorAssembly - If true, uses  tensor assembly 
%
%   blocksize         - If non-empty, divide the nodes in block with the given block size and proceed
%                       with assembly by iterating on the blocks.
%                       This is necessary in case of large models (get otherwise memory problems with MATLAB)
%
%                       Only used/available for tensor assembly version
%
%   neumann           - If true, set up the problem for Neumann boundary conditions (no flux). 
%                       In this case, a lighter implementation is used.
%
%                       Only used/available for tensor assembly version
%   
%
%   facetrans         - Accounts for face transmissibilities (see `computeMultiPointTransLegacy`)
%
%                       Only used/available for legacy version
%
%   invertBlocks      -
%               Method by which to invert a sequence of small matrices that
%               arise in the discretisation.  String.  Must be one of
%                  - MATLAB -- Use an function implemented purely in MATLAB
%                              (the default).
%
%                  - MEX    -- Use two C-accelerated MEX functions to
%                              extract and invert, respectively, the blocks
%                              along the diagonal of a sparse matrix.  This
%                              method is often faster by a significant
%                              margin, but relies on being able to build
%                              the required MEX functions.
% RETURNS:
%   vargout - Two cases : 
%  
%      - for legacy : [T, T_noflow]  - half-transmissibilities for each local face of each grid cell
%                                      in the grid. The number of half-transmissibilities equal the
%                                      number of rows in G.cells.faces. T_noflow gives the half-transmissibilities 
%                                      only for internal faces
% 
%      - for tensor assembly : [mpfastruct] - Assembly structure as computed by
%                                             computeMultiPointTransTensorAssembly (contains transmissibilities and
%                                             IndexArrays). If `neumann` is set to false, extra mappings are returned to 
%                                             handle the boundary.
%
% EXAMPLE:
%
% SEE ALSO:
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

   opt = struct('useTensorAssembly', false, ...
                'blocksize', [], ...
                'neumann', false);
   
   [opt, extra] = merge_options(opt, varargin{:});
   
   if (nargout > 1) || ~opt.useTensorAssembly
       % Uses legacy version
       [varargout{1:nargout}] = computeMultiPointTransLegacy(G, rock, extra{:});
       
   else
       % Uses tensor assembly version
       extra = [extra, {'blocksize', opt.blocksize, 'neumann', opt.neumann}];
       varargout{1} = computeMultiPointTransTensorAssembly(G, rock, extra{:});
       
   end
   


end
