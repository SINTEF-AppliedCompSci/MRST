function [name, output] = readEclipseOutputField(fid, varargin)
%Read next field of a formatted ECLIPSE output file
%
% SYNOPSIS:
%   [name, output] = readEclipseOutputField(fid)
%   [name, output] = readEclipseOutputField(fid, 'pn1', pv1, ...)
%
% PARAMETERS:
%   fid - Valid file identifier as obtained from FOPEN.  The file pointer
%         (FTELL(fid)) is assumed to be positioned directly before the next
%         keyword/field data.
%
% RETURNS:
%   name    - Name of this field.
%
%   output  - Field data.  A structure containing the following fields:
%              - type -- Data type of this field.
%                        String.  Supported values are
%                          - INTE -- Data values are Fortran INTEGERs.
%                          - REAL -- Data values are Fortran REALs.
%                          - DOUB -- Data values are Fortran DOUBLE PRECs.
%                          - LOGI -- Data values are Fortran LOGICALs.
%                          - CHAR -- Data values are character strings.
%
%              - values --
%                        Data values of this field.  A DOUBLE array for
%                        field data of type INTE, REAL, or DOUB.  A LOGICAL
%                        array for field data of type LOGI and a cell
%                        array of strings for field data of type CHAR.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - Verbose -- Whether or not to emit an informational
%                            message if a field with no associated data is
%                            encountered.
%                            Logical.  Default value: Verbose = FALSE (emit
%                            no message).
%
% NOTE:
%   This is a fairly low-level function which should generally not be
%   invoked directly from user code.
%
% SEE ALSO:
%   `readEclipseOutput`.

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


   opt = struct('Verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   lin              = fgetl(fid);
   lin(lin == '''') = '';

   a      = regexp(strtrim(lin), '\s+', 'split');
   name   =        a{1}       ;
   number = sscanf(a{2}, '%f');
   ttype  =        a{3}       ;

   switch lower(ttype),
      case {'inte', 'doub', 'real'},
         values = textscan(fid, '%f', number);
         values = values{1};
      case 'logi',
         values = textscan(fid, '%s', number);
         values = cellfun(@(s) s == 'T', values{1});    % -> LOGICAL
      case 'char',
         values = fscanf(fid, ' %*c%8c%*c', number);
         values = { cellstr(reshape(values, 8, []) .') };
      otherwise,
         if number == 0,
            dispif(opt.Verbose, ...
                   'Keyword ''%'' with no associated data.', name);
            output = struct('type', ttype, 'values', {});
            return
         else
            error(msgid('Type:Unknown'), ...
                  'Variable type ''%s'' is unexpected at this time.', ...
                  ttype);
         end
   end

   lin = fgetl(fid);  %#ok
   output = struct('type', ttype, 'values', values);
end
