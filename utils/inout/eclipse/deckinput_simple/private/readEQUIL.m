function equil = readEQUIL(fid, metric, varargin)
%Read EQUIL keyword
%
% SYNOPSIS:
%   equil = readEQUIL(fid, metric)
%   equil = readEQUIL(fid, metric, ntequil)
%
% PARAMETERS:
%   fid     - Valid file identifier as obtained from FOPEN.
%
%   metric  - Logical value indicating whether or not the input deck values
%             are given in metric or field units.  Typically corresponds to
%             'grdecl.METRIC'.
%
%   ntequil - Number of equilibration regions.  Integer.  Typically
%             corresponds to the ECLIPSE input deck parameter 'NTEQUIL'
%             entered in the keyword 'EQLDIMS' in the RUNSPEC section.
%             OPTIONAL.  Default value: ntequil = 1.
%
% RETURN:
%   equil - An 'ntequil'-by-11 array of DOUBLE values, the columns of which
%           correspond to the individual data items of the 'EQUIL' keyword
%           while the rows account for the individual equilibration
%           regions.
%
% NOTE:
%   The EQLIPSE edition of the EQUIL keyword only defines nine (9) columns.
%   The remaining two columns are included to support the FrontSim edition
%   of the EQUIL keyword.
%
% SEE ALSO:
%   `readDefaultedKW`.

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


   u = get_unit_system(metric);

   ntequil = 1;
   if nargin > 1 && isnumeric(varargin{1}) && numel(varargin{1}) == 1,
      ntequil = varargin{1};
   end

   template(1 : 11) = { '0' };
   template{10}     =   '1'  ;

   assert (ntequil >= 1);
   equil = readDefaultedKW(fid, template, 'NRec', ntequil);
   equil = cellfun(@(s) sscanf(s, '%f'), equil);

   depth = 1 : 2 : 6;
   press = 2 : 2 : 6;
   equil(:, depth) = convertFrom(equil(:, depth), u.depth);
   equil(:, press) = convertFrom(equil(:, press), u.press);
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function u = get_unit_system(metric)
   if metric,
      u = struct('depth', meter, ...
                 'press', barsa);
   else
      u = struct('depth', ft, ...
                 'press', psia);
   end
end
