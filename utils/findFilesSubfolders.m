function f = findFilesSubfolders(root)
%Find all files in a directory hierarchy
%
% SYNOPSIS:
%   files = findFilesSubFolders(root)
%
% PARAMETERS:
%   root - Full or relative path to root of directory hierarchy
%
% RETURNS:
%   files - Cell array of strings naming all files beneath 'root'.  The
%           pathnames include the 'root'.

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

   tmp = dir(root);

   subdir = { tmp([tmp.isdir]).name };
   subdir = subdir(~(strncmp(subdir, '.' , 1) | ...
                     strncmp(subdir, '..', 2)));

   fs = cellfun(@(s) findFilesSubfolders(fullfile(root, s)), ...
                sort(subdir), 'UniformOutput', false);

   f = cellfun(@(n) fullfile(root, n)    , ...
               { tmp(~[tmp.isdir]).name }, ...
               'UniformOutput', false);

   f = [ fs{:} , f ];
end
