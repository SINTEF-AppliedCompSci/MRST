function [xyz, zcorn] = grdeclXYZ(grdecl, varargin)
% grdeclXYZ get the pillar and zcorn coordinates of the cornerpoint format
%           in a more useful format.
%
% SYNOPSIS:
%   [xyz,zcorn]= grdeclXYZ(grdecl,varargin)
%   [xyz]= grdeclXYZ(grdecl,varargin)
%
% PARAMETERS:
%   grdecl   - Grid in Eclipse format
%
%
% RETURN:
%  xyz   - the pillar coordinates such that xyz(1:3,i,j) is the top
%          coordinates of pilar i,j
%  zcorn - vertical cordinates (z-value) for each corner for each cell,
%           ordered in increasing x, y, z.
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
% 'lefthanded_numbering' - should be put true if numbering of cells for the y
% axis is in decresing order
%
% SEE ALSO:
%

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

    opt = struct('lefthanded_numbering', false);
    opt = merge_options(opt, varargin{:});
    xyz = reshape(grdecl.COORD, [6, grdecl.cartDims(1:2) + 1]);
    if(opt.lefthanded_numbering)
        xyz([2, 5], :, :) = -xyz([2, 5], :, :);
    end
    zcorn = reshape(grdecl.ZCORN, 2 * grdecl.cartDims);

end
