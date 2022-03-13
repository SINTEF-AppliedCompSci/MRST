function sqform = convertAtlasToStruct(meta_thick, meta_top, data_thick, data_top)
%Create GRDECL struct from thickness/top data from the CO2 Storage Atlas 
%
% SYNOPSIS:
%   grdecl = convertAtlasTo3D(m_thick, m_top, d_thick, d_top, 3)
%
% DESCRIPTION:
%   Given two datasets with possibily non-matching nodes, interpolate and
%   combine to form a GRDECL struct suitable for 3D simulations after a
%   call to processGRDECL or for further manipulation.
%
% REQUIRED PARAMETERS:
%   meta_thick - Metainformation for the thickness data as produced from
%                readAAIGrid
%
%   meta_top   - Metainformation for top data.
%
%   data_thick - Data for the thickness data as produced from
%                readAAIGrid
%
%   top_thick  - Data for the top data as produced from
%                readAAIGrid
%
%
% RETURNS:
%   sqform - GRDECL struct suitable for processGRDECL
%
%
% NOTES:
%       It is likely easier to use getAtlasGrid instead of calling this
%       routine directly.
%
% SEE ALSO:
%   `getAlasGrid`

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

    ndims = [meta_top.dims];
    dims = ndims - 1;
    h = meta_top.cellsize;

    
    grdecl.cartDims = reshape(dims, 1, []);
    
    xl = meta_top.xllcorner;
    yl = meta_top.yllcorner;
    
    % Create grids
    [X, Y]  = ndgrid(linspace(xl, dims(1)*h + xl, ndims(1)), ...
                        linspace(yl, dims(2)*h + yl, ndims(2)));
    

    F_top   = interpolateData(meta_top,   data_top);
    F_thick = interpolateData(meta_thick, data_thick);
    
    x = reshape(X, [], 1);
    y = reshape(Y, [], 1);
    thick = reshape(F_thick(x, y), ndims(1), ndims(2));
    top   = reshape(F_top(x, y)  , ndims(1), ndims(2));
   
    sqform=struct('nodeDims',ndims,'thick',thick,'top',top, ...
       'dx',h,'dy',h,'orig',[xl,yl],'X',X,'Y',Y);
    
    
   
end

function F = interpolateData(meta, data)
    dims = meta.dims;
    % We have a cell centered grid, so subtract by one
    gdims = dims - 1;
    h = meta.cellsize;
    xl = meta.xllcorner;
    yl = meta.yllcorner;
    [X, Y]  = ndgrid(linspace(xl, gdims(1)*h + xl, dims(1)), ...
                     linspace(yl, gdims(2)*h + yl, dims(2)));    
    F =@(x,y) interp2(X',Y',data', x, y);              
    %F = TriScatteredInterp(X(:), Y(:), data(:));
end
