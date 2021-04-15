function rock = readSimRock(casename)
%Read SAM simulator ROCK input files.
%
% SYNOPSIS:
%   rock = readSimRock(casename)
%
% PARAMETERS:
%   casename - Simulation case base name.  String assumed to contain the
%              (relative) base name of a simulation case.
%
%              The rock parameters of a simulation case is defined by two
%              input files:
%                 Permeability: [caseName, '-perm.dat']
%                 Porosity:     [caseName, '-poro.dat']
%
% RETURNS:
%   rock - MRST rock structure.  Permeabilities measured in units of m^2
%          (standard SI units).
%
% SEE ALSO:
%   `readSimData`.

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


   rock = struct('perm', read([casename, '-perm.dat']), ...
                 'poro', read([casename, '-poro.dat']));

   assert (size(rock.perm,1) == size(rock.poro,1), ...
          ['Size of permeability data does not match ', ...
           'size of porosity data']);

   % Permeability in standard SI units (mD in on-disk format).
   rock.perm = convertFrom(rock.perm, milli*darcy);
end

%--------------------------------------------------------------------------

function v = read(fn)
   [fid, msg] = fopen(fn);

   if fid < 0,
      error('Failed to open rock data file ''%s'': %s', fn, msg);
   end

   % Rock data format
   % N [ number of lines/records ]
   % rec1
   % rec2
   % ...
   % recN
   v = fscanf(fid, '%f');

   fclose(fid);

   v = reshape(v(2:end), [], v(1)) .';
end
