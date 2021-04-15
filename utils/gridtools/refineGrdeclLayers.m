function ngrdecl = refineGrdeclLayers(grdecl, layers, ref)
% Refine a GRDECL structure in the vertical direction
%
% SYNOPSIS:
%   function ngrdecl = refineGrdeclLayers(grdecl, layers, ref)
%
% DESCRIPTION:
%   Refine corner-point grid in vertical direction
%
% PARAMETERS:
%   grdecl   - Grid in Eclipse format
%   layers   - Refine vertical layers from layers(1) to layers(2)
%   ref      - Vertical refinement factor for each cell
%
% RETURNS:
%   ngrdecl - Refined grid in eclipse format.
%
%
% SEE ALSO:
%   `refineGrdecl`

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

    [xyz, zcorn] = grdeclXYZ(grdecl);
    zcorn1 = zcorn(:, :, 1:2*(layers(1)-1));
    zcorn2 = zcorn(:, :, 2*(layers(2)-1)+1:end);
    cgrdecl = cutGrdecl(grdecl, [1 grdecl.cartDims(1);1 grdecl.cartDims(2);layers]);
    cgrdecl = refineGrdecl(cgrdecl, [1 1 ref]);
    [cxyz, czcorn] = grdeclXYZ(cgrdecl);
    ol = (layers(2)-layers(1)+1);
    newl = ol*ref;
    newcartdim = [grdecl.cartDims(1:2), grdecl.cartDims(3)-ol+newl];
    zcorn = nan(newcartdim*2);
    zcorn(:, :, 1:2*(layers(1)-1)) = zcorn1;
    zcorn(:, :, end-size(zcorn2, 3)+1:end) = zcorn2;
    zcorn(:, :, 2*(layers(1)-1)+1:(end-size(zcorn2, 3)+2)) = czcorn;
    ngrdecl = struct('cartDims', newcartdim, 'COORD', grdecl.COORD, 'ZCORN', zcorn);

end
