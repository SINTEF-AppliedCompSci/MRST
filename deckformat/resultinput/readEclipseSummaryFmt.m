function [smry, smspec] = readEclipseSummaryFmt(prefix)
%Read formatted (text/ASCII) ECLIPSE summary data
%
% SYNOPSIS:
%   [summary, smspec] = readEclipseSummaryFmt(prefix)
%
% PARAMETERS:
%   prefix - Path-name prefix from which to construct list of summary file
%            names.  Specifically, this function reads files which match
%            the regular expressions
%
%                [prefix, '\.FSMSPEC']  (formatted summary specification)
%                [prefix, '\.A\d{4}']   (formatted summary files)
%
%            Use function 'readEclipseSummaryUnFmt' to read unformatted
%            (binary) summary data.
%
% RETURNS:
%   summary - Summary data structure.  All MINISTEP results of an ECLIPSE
%             run concatenated together.  Also includes an additional
%             sub-structure, UNITS, whose fields contain a textual
%             representation of unit of measurement of the corresponding
%             summary field such as
%
%                 summary.UNITS.TIME  = 'DAYS'
%                 summary.UNITS.WBHP  = 'BARSA'
%                 summary.UNITS.WOPR  = 'SM3/DAY'
%                 summary.UNITS.YEARS = 'YEARS'
%                 summary.UNITS.TCPU  = 'SECONDS'
%
%   smspec  - Summary specifiction obtained from the '.FSMSPEC' file.
%             Contains MINISTEP times and all summary vectors declared in
%             the run deck.  Additionally contains a field, '.RptTime' that
%             specifies the times at which restart files (3D data) is
%             reported.
%
% SEE ALSO:
%   `readEclipseSummaryUnFmt`.

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


   is_open_pre = fopen('all');

   [dname, fp] = fileparts(prefix);
   if isempty(dname),
      dname = '.';
   end

   smspec = readEclipseOutputFileFmt([prefix, '.FSMSPEC']);

   summaryfiles = matchResultFiles(dname, [fp, '\.A\d{4}']);

   smry = readEclipseSummary(summaryfiles, smspec,  ...
                             @(fn) fopen(fn, 'rt'), ...
                             @readFieldFmt);

   is_open_post = fopen('all');

   assert (all(size(is_open_pre) == size(is_open_post)) && ...
           all(is_open_pre == is_open_post),               ...
           'Summary reader leaks file identifiers.');
end
