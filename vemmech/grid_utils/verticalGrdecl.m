function grdecl = verticalGrdecl(grdecl, varargin)
%Transform GRDECL pillars into vertical pillars
%
% SYNOPSIS:
%   function grdecl = verticalGrdecl(grdecl, varargin)
%
% DESCRIPTION: Straightened up the pillars in a corner point grid, given in
% Eclipse, and make them vertical
%
% PARAMETERS:
%   grdecl   - Grid structure in Eclipse format
%
% OPTIONAL PARAMETERS:
%   method - 'all' : all pillars are made vertical
%
%            'sides' : Only the pillare on the sides are made vertical
%
%
% RETURNS:
%   grdecl - Grid structure with vertical pillars in Eclipse format.
%
% EXAMPLE:
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

    opt = struct('method', 'all');
    opt = merge_options(opt, varargin{:});
    top = [1, 2];
    bot = [4, 5];
    [xyz, zcorn] = grdeclXYZ(grdecl);
    switch opt.method
      case 'sides'
        xyz(bot, 1, :) = xyz(top, 1, :); % x-sides
        xyz(bot, end, :) = xyz(top, end, :);
        xyz(bot, :, 1) = xyz(top, :, 1); % y-sides
        xyz(bot, :, end) = xyz(top, :, end);
      case 'all'
        xyz(top, :, :) = xyz(bot, :, :); 
      otherwise
        error('no such method');     
    end
    grdecl.COORD = xyz(:);
end