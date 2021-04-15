function rstrt = readEclipseRestart(files, rstReader, opt)
%Read set of ECLIPSE restart files.
%
% SYNOPSIS:
%   restart = readEclipseRestart(files, rstReader, opt)
%
% PARAMETERS:
%   files     - List, represented as a cell array of strings, of files from
%               which to extract restart data.
%
%   rstReader - Restart-file reader.  Assumed to support the syntax
%
%                   item = rstReader(filename)
%
%               whence the restart data contained in the resource (an
%               on-disk file) denoted by 'filename' will be read into a
%               structure ('item') whose fields correspond to the
%               individual restart data in the resource.
%
%               Typically corresponds to one of @readEclipseOutputFileFmt
%               or @readEclipseOutputFileUnFmt.
%
%   opt       - Subject to change.
%
% RETURNS:
%   restart - Restart data structure.  One field for each restart data item
%             in the set of restart files.  Individual restart data items
%             are stored in separate columns of the corresponding cell
%             array.
%
% SEE ALSO:
%   `private/readEclipseSummary`, `readEclipseOutputFileFmt`.

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

   rstrt = [];

   for f = reshape(files, 1, [])
      item = rstReader(f{1});

      rstrt = append_restart_data(rstrt, item, opt);
   end
end

%--------------------------------------------------------------------------

function rstrt = append_restart_data(rstrt, item, opt)
   if ~isempty(opt.RestartFields)
      fields = opt.RestartFields;
   else
      fields = fieldnames(item);
   end

   for f = reshape(fields, 1, [])
      fld = f{1};

      if ~isfield(rstrt, fld), rstrt.(fld) = {}; end

      if  isfield(item, fld)
         rstrt.(fld) = [rstrt.(fld), { item.(fld).values }];
      end
   end
end
