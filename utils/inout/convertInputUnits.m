function grdecl = convertInputUnits(grdecl, u)
%Convert input data to MRST's strict SI conventions
%
% SYNOPSIS:
%   grdecl = convertInputUnits(grdecl, u)
%
% PARAMETERS:
%   grdecl - Input data structure as defined by function 'readGRDECL'.
%
%   u      - Input unit convention data structure as defined by function
%            'getUnitSystem'.  It is the caller's responsibility to request
%            a unit system that is compatible with the actual data
%            represented by 'grdecl'.
%
% RETURNS:
%   grdecl - Input data structure whos fields have values in consistent
%            units of measurements (e.g., permeability measured in m²,
%            lengths/depths measured in meters &c).
%
% SEE ALSO:
%   `readGRDECL`, `getUnitSystem`, `convertFrom`.

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

   for kw = reshape(fieldnames(grdecl), 1, [])
      key = kw{1};

      switch key
         case {'PERMX' , 'PERMXY', 'PERMXZ', ...
               'PERMYX', 'PERMY' , 'PERMYZ', ...
               'PERMZX', 'PERMZY', 'PERMZ' , 'PERMH' }
            grdecl.(key) = convertFrom(grdecl.(key), u.perm);

         case {'DXV'   , 'DYV'   , 'DZV'   , 'DEPTHZ', ...
               'COORD' , 'ZCORN'                     }
            grdecl.(key) = convertFrom(grdecl.(key), u.length);

         case {'TRANX', 'TRANY', 'TRANZ'}
            grdecl.(key) = convertFrom(grdecl.(key), u.trans);

         case {'cartDims',                        ...  % MRST specific
               'UnhandledKeywords',               ...  % MRST specific
               'ROCKTYPE',                        ...  % MRST specific
               'ACTNUM', 'NTG', 'PORO', 'SATNUM', ...
               'FAULTS', 'MULTFLT',               ...
               'MULTX' , 'MULTX_' ,               ...
               'MULTY' , 'MULTY_' ,               ...
               'MULTZ' , 'MULTZ_' ,               ...
               'VSH'   , 'GDORIENT' }
            continue;  % No unit conversion needed.

         otherwise
            error('No known unit conversion for ''%s''.', key);
      end
   end
end
