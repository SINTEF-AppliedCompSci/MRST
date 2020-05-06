function [T, T_noflow] = computeMultiPointTrans(G, rock, varargin)
% Compute multi-point transmissibilities for MPFA using local flux mimetic approach
%
% SYNOPSIS:
%   function [T, T_noflow] = computeMultiPointTrans(G, rock, varargin)
%
% DESCRIPTION:
%
% Two versions available : 'legacy' (default) and 'tensor Assembly' (use key option `useTensorAssembly`).
%
% This legacy version is faster. It is limited to a mesh with grid cells where
% corners have the same number of faces as the spatial dimension (this is always
% the case in 2D but not in 3D). The tensor assembly version
% (computeMultiPointTransTensorAssembly) can handle the other cases but is
% slower (the implementation will be optimized in the future to run faster).
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
%   useTensorAssembly - Uses the tensor assembly implemention
%
%   blocksize         - If non-empty, divide the nodes in block with the given block size and proceed
%                       with assembly by iteration on the blocks.
%                       This is necessary in case for large models (get otherwise memory problems with MATLAB)
%                       Only used/available for tensor assembly version
%
%   neumann           - If true, set up the problem for Neumann boundary conditions (no flux). 
%                       In this case, a lighter implementation is used.
%                       Only used/available for tensor assembly version
%   
%
%   facetrans         - Accounts for face transmissibilities (see `computeMultiPointTransLegacy`)
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
%   T - half-transmissibilities for each local face of each grid cell
%       in the grid.  The number of half-transmissibilities equal the
%       number of rows in G.cells.faces.
%
% EXAMPLE:
%
% SEE ALSO:
%

   opt = struct('useTensorAssembly', false, ...
                'blocksize', [], ...
                'neumann', false);
   
   [opt, extra] = merge_options(opt, varargin{:});
   
   if ~opt.useTensorAssembly
       [T, T_noflow] = computeMultiPointTransLegacy(G, rock, extra{:});
   else
       extra = [extra, {'blocksize', opt.blocksize, 'neumann', opt.neumann}];
       mpfastruct = computeMultiPointTransTensorAssembly(G, rock, extra{:});
       T = mpfastruct;
   end       

end
