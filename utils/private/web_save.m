function file = web_save(file, url, varargin)
%Restricted fall-back option for WEBSAVE on older releases of MATLAB
%
% SYNOPSIS:
%   file = web_save(file, url)
%   file = web_save(file, url, 'p1', v1, ...)
%
% PARAMETERS:
%   file - Name of local file into which contents of remote resource will
%          be saved.
%
%   url  - Uniform resource locator of remote resource (file contents).
%
%   'pn'/pv -
%          List of key/value pairs that will be passed through to the
%          underlying GET request of the URL.  Numeric arguments converted
%          to strings using function NUM2STR.
%
% RETURNS:
%   file - Unmodified input file name if successful, empty in case of
%          download failure.
%
% NOTE:
%   This function is only meaningful in the restricted context of an
%   encompassing download manager such as function githubDownload.
%
% SEE ALSO:
%   githubDownload.

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

   args = stringify_values(varargin{:});
   
   if ~isempty(args),
      args = [ { 'get' }, { args } ];
   end

   s = urlread(url, args{:});

   [fid, msg] = fopen(file, 'wb');

   if fid < 0,
      error('Open:Fail', 'Failed to open file ''%s'': %s', file, msg)
   end

   fprintf(fid, '%s', s);

   fclose(fid);
end

%--------------------------------------------------------------------------

function args = stringify_values(varargin)
   vals = varargin(2 : 2 : end);
   
   i = cellfun(@(x) isnumeric(x) || islogical(x), vals);
   strng = cellfun(@num2str, vals(i), 'UniformOutput', false);
   
   [vals{i}] = strng{:};

   args = reshape([ varargin(1 : 2 : end) ; vals ], 1, []);
end
