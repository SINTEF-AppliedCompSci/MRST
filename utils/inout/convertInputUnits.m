function grdecl = convertInputUnits(grdecl, u, varargin)
%Convert input data to MRST's strict SI conventions
%
% SYNOPSIS:
%   grdecl = convertInputUnits(grdecl, u)
%   grdecl = convertInputUnits(grdecl, u, 'pn1', pv1, ...)
%
% PARAMETERS:
%   grdecl - Input data structure as defined by function 'readGRDECL'.
%
%   u      - Input unit convention data structure as defined by function
%            'getUnitSystem'.  It is the caller's responsibility to request
%            a unit system that is compatible with the actual data
%            represented by 'grdecl'.
%
% OPTIONAL PARAMETERS:
%   exclude - A cell array of keywords to exclude.
%
% RETURNS:
%   grdecl - Input data structure whos fields have values in consistent
%            units of measurements (e.g., permeability measured in m²,
%            lengths/depths measured in meters &c).
%
% SEE ALSO:
%   `readGRDECL`, `getUnitSystem`, `convertFrom`.

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

   opt = struct('exclude', {{}});
   opt = merge_options(opt, varargin{:});

   keywords = reshape(fieldnames(grdecl), 1, []);
   [~, idx] = intersect(keywords, opt.exclude);
   if ~isempty(idx)
      dispif(mrstVerbose, 'Excluding keyword ''%s''.\n', keywords{idx});
   end
   keywords(idx) = [];

   for kw = keywords
      key = kw{1};

      switch key
         case {'PERMX' , 'PERMXY', 'PERMXZ', ...
               'PERMYX', 'PERMY' , 'PERMYZ', ...
               'PERMZX', 'PERMZY', 'PERMZ' , 'PERMH' }
            grdecl.(key) = convertFrom(grdecl.(key), u.perm);

         case {'DXV'   , 'DYV'   , 'DZV'   , 'DEPTHZ', ...
               'COORD' , 'ZCORN'                     }
            grdecl.(key) = convertFrom(grdecl.(key), ...
                                       length_unit(grdecl, u));

         case {'TRANX', 'TRANY', 'TRANZ'}
            grdecl.(key) = convertFrom(grdecl.(key), u.trans);

         case 'MAPAXES'
            grdecl.(key) = convertFrom(grdecl.(key), map_units(grdecl, u));

         case {'cartDims',                        ...  % MRST specific
               'UnhandledKeywords',               ...  % MRST specific
               'ROCKTYPE',                        ...  % MRST specific
               'ACTNUM', 'NTG', 'PORO', 'SATNUM', ...
               'FAULTS', 'MULTFLT',               ...
               'MULTX' , 'MULTX_' ,               ...
               'MULTY' , 'MULTY_' ,               ...
               'MULTZ' , 'MULTZ_' ,               ...
               'VSH'   , 'GDORIENT',              ...
               'MAPUNITS', 'GRIDUNIT' }
            continue;  % No unit conversion needed.

         otherwise
            error('Unknown:Unit', ...
                  'No known unit conversion for ''%s''.', key);
      end
   end
end

%--------------------------------------------------------------------------

function unit = length_unit(grdecl, u)
   unit = u.length;
   if isfield(grdecl, 'GRIDUNIT') && ~isempty(grdecl.GRIDUNIT) && ...
         ~isempty(grdecl.GRIDUNIT{1})
      unit = umap(grdecl.GRIDUNIT{1});
   end
end

%--------------------------------------------------------------------------

function unit = map_units(grdecl, u)
   unit = u.length;
   if isfield(grdecl, 'MAPUNITS')
      unit = umap(grdecl.MAPUNITS);
   end
end

%--------------------------------------------------------------------------

function u = umap(mapunits)
   if ischar(mapunits)
      switch lower(mapunits)
         case { 'metres', 'meters' }
            u = meter;
         case { 'foot', 'feet', 'ft' }
            u = ft;
         case { 'centimetre', 'centimeter', 'cm' }
            u = centi*meter;
         otherwise
            warning('MAPUNITS:Unsupported', ...
                    'Unsupported MAPUNITS ''%s'' - Using ''metre''', ...
                    mapunits);

            u = meter;
      end
   else
      u = 1;
   end
end
