function info = refineInfo(info_old, ref); 

   info = info_old; 
   meta_old = info_old.meta; 

   % refine meta
   meta = meta_old; 
   meta.ncols = (meta_old.ncols - 1) * ref + 1; 
   meta.nrows = (meta_old.nrows - 1) * ref + 1; 
   meta.cellsize = meta_old.cellsize / ref; 
   meta.dims = [meta_old.dims - 1] * ref + 1; 

   data_old = info_old.data; 

   [Xn, Yn] = meshgrid(linspace(1, meta_old.ncols, meta.ncols), linspace(1, meta_old.nrows, meta.nrows)); 
   data = interp2(data_old, Xn, Yn); 
   % clear Xn; 
   % clear Yn; 
   info.meta = meta; 
   info.data = data; 
end

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