function W = addWellVE(W, Gt, rock2D, cell, varargin)
% Add well directly to a top surface grid
%
% SYNOPSIS:
%   function W = addWellVE(W, Gt, rock2D, cell, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   W        - Well structure or empty if no other well exist.  Updated upon
%              return. 
%   Gt       - top surface grid data structure
%   rock2D   - rock structure associated with the top surface grid
%   cell     - (unique) grid cell where the well should be added
%   varargin - 
%
% RETURNS:
%   W - Updated (or freshly created) well structure.
%
% SEE ALSO:
% `convertwellsVE`

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

   opt = struct('innerProduct', 'ip_tpf',                     ...
                'name'        , sprintf('W%d', numel(W) + 1), ...
                'radius'      , 0.1,                          ...
                'type'        , 'bhp',                        ...
                'val'         , 0,                            ...
                'comp_i'      , [0 1],                        ...
                'compi'       , [],                           ....
                'refDepth'    , [],                            ...
                'sign'        , 0);
   
   opt = merge_options(opt, varargin{:});
   
   Wnew = addWell([], Gt, rock2D, cell, ...
                  'type'         , opt.type          , ...
                  'val'          , opt.val           , ...
                  'radius'       , opt.radius        , ...
                  'comp_i'       , opt.comp_i        , ...
                  'compi'        , opt.compi         , ...
                  'innerProduct' , opt.innerProduct  , ...
                  'name'         , opt.name          , ...
                  'refDepth'     , opt.refDepth      , ...
                  'sign'         , opt.sign);              

   Wnew(end).WI = Wnew(end).WI .* Gt.cells.H(cell);
   Wnew(end).h = Gt.cells.H(cell);
   Wnew(end).dZ = Gt.cells.H(cell) * 0;
   
   W = [W Wnew];
end
