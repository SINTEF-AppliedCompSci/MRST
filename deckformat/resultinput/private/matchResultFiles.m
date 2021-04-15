function fn = matchResultFiles(dname, patt)
%Match a set of files in a given directory to specified pattern
%
% SYNOPSIS:
%   files = matchResultFiles(dname, pattern)
%
% PARAMETERS:
%   dname   - Directory in which to match file names.  String.
%
%   pattern - Regular expression againts which the files of directory
%             'dname' will be matched.
%
% RETURNS:
%   files - List, represented as a cell array of strings, of file names in
%           the directory 'dname', that match the regular expression
%           pattern.  The list is lexicographically sorted.
%
% NOTE:
%   This is an internal function in the ECLIPSE result reader.  Its
%   semantics or existence may change without prior notice.
%
% SEE ALSO:
%   `dir`, `regexp`.

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


   d  = dir(dname);
   fn = regexp({ d.name }, patt, 'match', 'once');
   i  = ~cellfun(@isempty, fn);
   fn = strcat([dname, filesep], sort(fn(i)) .');
end
