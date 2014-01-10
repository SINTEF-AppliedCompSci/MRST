function ok = make_spe10_data
%Create on-disk (MAT file) representation of SPE 10 'rock' data.
%
% SYNOPSIS:
%   ok = make_spe10_data
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   ok - Status flag indicating successful creation of on-disk 'rock' data.
%
% NOTE:
%   The on-disk representation is a 'rock' structure stored in the file
%   'spe10_rock.mat' in the directory of function 'make_spe10_data'.
%
% SEE ALSO:
%   SPE10_rock.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   d      = fileparts(mfilename('fullpath'));
   fname  = @(f) fullfile(d, f);
   fexist = @(f) exist(fname(f), 'file') == 2;

   ok = fexist('spe10_rock.mat');

   if ~ok,
      % Data does not already exist (in MAT form) on disk.  Download and/or
      % extract official ZIP file and/or processes already unpacked ZIP
      % file contents.

      ok = fexist('spe_perm.dat') && fexist('spe_phi.dat');

      if ~ok,
         % Data files not available in unpacked form.  Extract local
         % archive *or* the archive file downloaded from the SPE web site.

         zipfile = 'por_perm_case2a.zip';
         if fexist(zipfile),
            url = fname(zipfile);
         else
            url = ['http://www.spe.org/web/csp/datasets/', zipfile];
         end

         dispif(mrstVerbose, ...
               ['Please wait while the second SPE10 ', ...
                'dataset is downloaded...']);
         unzip(url, d);
         dispif(mrstVerbose, 'Done\n');

         ok = fexist('spe_perm.dat') && fexist('spe_phi.dat');
      end

      if ok,
         % We should *always* end up here.

         %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         % Permeability data.
         %
         [fid, msg] = fopen(fname('spe_perm.dat'), 'rt');
         if fid < 0, error(msg); end

         rock.perm = reshape(fscanf(fid, '%f'), [], 3);
         fclose(fid);

         %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         % Porosity data.
         %
         [fid, msg] = fopen(fname('spe_phi.dat'), 'rt');
         if fid < 0, error(msg); end

         rock.poro = reshape(fscanf(fid, '%f'), [], 1);
         fclose(fid);

         %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         % Verify size of input data.
         %
         ncell = 60 * 220 * 85;

         assert (numel(rock.perm) == 3 * ncell);
         assert (numel(rock.poro) == 1 * ncell);

         %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         % Save data in form more amenable to subsequent M processing.
         %
         save(fname('spe10_rock'), 'rock')
      end
   end
end
