function downloadMBGL
%Download MATLAB BGL package from The MathWorks FileExchange
%
% SYNOPSIS:
%   downloadMBGL
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   Nothing.

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

   % Note: These URL fragments are likely to be unstable and must be
   % checked and/or updated from time to time.
   dl  = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads';
   pth = 'submissions/10922/versions/2/download';
   zip = 'zip';
   url = [dl, '/', pth, '/', zip];

   % Ouput directory.
   dest = fullfile(ROOTDIR, 'utils', '3rdparty', 'matlab_bgl');

   try
      fprintf('Downloading MatlabBGL ZIP archive ... ')
      t = tic;
      unzip(url, dest);
      t = toc(t);
      fprintf('done (%.02f [s])\n\n', t);
   catch
      fex       = 'https://www.mathworks.com/matlabcentral/fileexchange';
      fallback  = [fex, '/', '10922-matlabbgl'];
      print_url = ['<a href="', fallback, '">', fallback, '</a>'];

      error(['Failed to retrieve MatlabBGL package from known URL\n', ...
             'Please see fall-back location\n  * %s'], print_url);
   end

   if exist('assert', 'builtin')
      % MatlabBGL contains an 'assert.m' file in its 'test' directory that
      % conflicts with the built-in function of the same name.  Just use
      % the built-in, because the semantics of 'assert.m' are the same as
      % for ASSERT.
      %
      fprintf(['Patching MatlabBGL for compatibility ', ...
               'with recent MATLAB ... ']);
      delete(fullfile(dest, 'matlab_bgl', 'test', 'assert.m'));
      fprintf('done\n');
   end

   % Rewrite the module path for matlab_bgl to point to the correct
   % directory
   [fid, msg] = fopen(fullfile(ROOTDIR, 'startup_user.m'), 'at');
   if fid < 0
      warning('Open:Fail', 'Failed to open ''startup_user.m'': %s', msg);
   else
      fprintf(['Patching ''startup_user'' function for ', ...
               'MatlabBGL installation ... ']);

      fprintf(fid, 'mrstPath reregister matlab_bgl ...\n    ''%s''\n', ...
              fullfile(dest, 'matlab_bgl'));

      fclose(fid);
      fprintf('done\n');
   end

   fprintf('Update module mapping to new location of MatlabBGL ... ');
   mrstPath('reregister', 'matlab_bgl', fullfile(dest, 'matlab_bgl'));
   fprintf('done\n');

   % Create a modload.m that's compatible with the new, third-party
   % component.
   copyfile(fullfile(dest, 'private', 'modload.m.in'), ...
            fullfile(dest, 'private', 'modload.m'));

   clear modload modload_fallback
end
