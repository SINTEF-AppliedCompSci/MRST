function [rstrt, rsspec] = readEclipseRestartFmt(prefix, varargin)
%Read formatted (text/ASCII) ECLIPSE restart data
%
% SYNOPSIS:
%   [restart, rsspec] = readEclipseRestartFmt(prefix)
%
% PARAMETERS:
%   prefix - Path-name prefix from which to construct list of summary file
%            names.  Specifically, this function reads files which match
%            the regular expressions
%
%                [prefix, '\.FRSSPEC']  (formatted restart specification)
%                [prefix, '\.F\d{4}']   (formatted restart files)
%
%            Use function 'readEclipseRestartUnFmt' to read unformatted
%            (binary) restart data.
%
% RETURNS:
%   restart - Restart data structure.  One field for each restart data item
%             in the set of restart files.  Individual restart data items
%             are stored in separate columns of the corresponding cell
%             array.
%
%   rsspec  - Restart specifiction obtained from the '.FRSSPEC' file.
%
% SEE ALSO:
%   `readEclipseSummaryFmt`, `readEclipseRestartUnFmt`.

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


   opt = struct('RestartFields', {{}});
   opt = merge_options(opt, varargin{:});

   is_open_pre = fopen('all');

   [dname, fp] = fileparts(prefix);
   if isempty(dname),
      dname = '.';
   end

   rstReader = @readEclipseOutputFileFmt;
   rsspec    = rstReader([prefix, '.FRSSPEC']);

   rstfiles  = matchResultFiles(dname, [fp, '\.F\d{4}']);
   rstrt     = readEclipseRestart(rstfiles, rstReader, opt);

   is_open_post = fopen('all');

   assert (all(size(is_open_pre) == size(is_open_post)) && ...
           all(is_open_pre == is_open_post),               ...
           'Restart reader leaks file identifiers.');
end
