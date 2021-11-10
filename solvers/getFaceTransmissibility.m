function T = getFaceTransmissibility(G, rock, varargin)
%Compute face transmissibilities, accounting for input-specific multipliers
%
% SYNOPSIS:
%    T = getFaceTransmissibility(G, rock)
%    T = getFaceTransmissibility(G, rock,       'pn1', pv1, ...)
%    T = getFaceTransmissibility(G, rock, mult)
%    T = getFaceTransmissibility(G, rock, mult, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Computes transmissibilities per interface.  Incorporates multipliers if
%   provided by caller.
%
% PARAMETERS:
%   G    - Valid grid structure.
%
%   rock - Valid rock structure.  May contain transmissibility multiplier
%          information as derived using function `grdecl2Rock` or similar.
%          Those multipliers will be applied directly to the resulting face
%          transmissibility values returned from `getFaceTransmissibility`.
%
%   mult - User-defined multiplier structure.  This is optional.  If you
%          don't need any multipliers in addition to those that exist in
%          the `rock` structure, you do not need to pass this parameter.
%          Function `getFaceTransmissibility` supports both one-sided
%          ("half face") multipliers that are applied before reducing the
%          transmissibility to per-face values and per-face multipliers
%          that are applied after the reduction to per-face values.
%
%          We support two primary multiplier source types
%            - mult is numeric, shape numel(T)-by-n.  In this case, the
%            number of rows is used to infer whether the array contains
%            one-sided or per-face multipliers.
%
%            - mult is structure
%                - Field `onesided`, if present, is one of the following
%                types
%                  * numeric, shape size(G.cells.faces,1)-by-n
%
%                  * structure with one or both of the field types
%                     - direct: numeric, shape size(G.cells.faces,1)-by-n
%                     - face and value: Sparse representation.  Elements of
%                     'face' are treated as face indices and the elements
%                     of 'value' are treated as numeric multiplier values.
%
%                - Field `face`, if present, is one of the following types
%                  * numeric, shape G.faces.num-by-n
%
%                  * structure with one or both of the field types
%                     - direct: numeric, shape G.faces.num-by-n
%                     - cellface and value: Sparse representation. Elements
%                     of 'cellface' are treated as cell-face indices and
%                     the elements of 'value' are treated as multiplier
%                     values.
%
% OPTIONAL PARAMETERS:
%   All additional 'key'/value parameters are passed directly on to
%   underlying function `computeTrans`.
%
% RETURNS:
%   T - Transmissibilities, one per interface, including boundary faces if
%       present in the input grid `G`.
%
% SEE ALSO:
%   `computeTrans`, `computeGeometry`.

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

   multargname = '';
   if nargin > 2, multargname = inputname(3); end
   [mult, varargin] = get_multipliers(multargname, G, rock, varargin{:});

   T = compute_onesided_trans(G, rock, varargin{:});

   T = apply_onesided_multipliers(T, mult);
   T = apply_harmonic_reduction(T, G);
   T = apply_face_multipliers(T, mult);
end

%--------------------------------------------------------------------------

function [mult, varargin] = get_multipliers(multargname, G, rock, varargin)
   mult = [];

   if mod(numel(varargin), 2) ~= 0
      if is_multiplier_type(varargin{1})
         mult     = varargin{1};
         varargin = varargin(2 : end);

         if is_deck_structure(mult)
            warning('Syntax:Deprecated', ...
            	    ['The syntax\n\n  T = %s(G, rock, deck)\n\nis ', ...
                    'deprecated and will be removed in a future ', ...
                    'version of MRST.\nPlease switch to putting ', ...
                    'transmissiblity multiplier data in the ''rock''', ...
                    'structure.'], mfilename());

            if ~any(isfield(rock, {'multipliers', 'faultdata'}))
               mult = extract_gridsection_multipliers(G, mult.GRID);
            else
               warning('Multipliers:AccountedFor', ...
                      ['Transmissibility multipliers already ', ...
                       'accounted for in ''rock'' structure']);
               mult = [];
            end
         end

      else
         if isempty(multargname) || all(isspace(multargname))
            multargname = '<expression>';
         end

         error('Multiplier:Isnt', ...
              ['Multiplier argument ''%s'' of type ''%s'' is not a ', ...
               'valid multiplier type'], multargname, class(varargin{1}));
      end
   end

   if ~ (isempty(varargin) || is_key_type(varargin(1 : 2 : end)))
      error('Unexpected:Input', ...
            'Optional parameters must be provided as ''key''/value pairs');
   end

   if any(isfield(rock, {'multipliers', 'faultdata'}))
      if size(G.cells.faces, 2) < 2
          warning('FaceDirection:Missing', ...
                 ['Unable to process faults - G.cells.faces does not ', ...
                  'contain direction data as second column.']);
      else
          mult = incorporate_rock_multipliers(mult, G, rock);
      end
   end
end

%--------------------------------------------------------------------------

function T = compute_onesided_trans(G, rock, varargin)
   if has_cpgeometry(G)
      cpgeom = G.cells.cpgeometry;

      T = computeTrans(G, rock, varargin{:}, ...
                       'K_system'       , 'loc_xyz', ...
                       'cellCenters'    , cpgeom.centroids, ...
                       'cellFaceCenters', cpgeom.facecentroids);
   else
      T = computeTrans(G, rock, varargin{:});
   end
end

%--------------------------------------------------------------------------

function T = apply_harmonic_reduction(T, G)
   T = 1 ./ accumarray(G.cells.faces(:,1), 1 ./ T, [G.faces.num, 1]);
end

%--------------------------------------------------------------------------

function T = apply_onesided_multipliers(T, mult)
   T = apply_multipliers(T, mult, 'onesided', { 'cellface', 'value' });
end

%--------------------------------------------------------------------------

function T = apply_face_multipliers(T, mult)
   T = apply_multipliers(T, mult, 'face', { 'face', 'value' });
end

%--------------------------------------------------------------------------

function tf = is_multiplier_type(mult)
   tf = isempty(mult) || isnumeric(mult) || ...
      is_deck_structure(mult) || ...
      (isstruct(mult) && any(isfield(mult, { 'onesided', 'face' })));
end

%--------------------------------------------------------------------------

function tf = is_deck_structure(mult)
   tf = ~isempty(mult) && isstruct(mult) && ...
      all(isfield(mult, {'RUNSPEC', 'GRID'}));
end

%--------------------------------------------------------------------------

function mult = extract_gridsection_multipliers(G, griddata)
   mult = [];
   if isempty(griddata), return, end

   rock = grdecl2Rock(griddata, G.cells.indexMap);
   mult = incorporate_rock_multipliers(mult, G, rock);
end

%--------------------------------------------------------------------------

function mult = incorporate_rock_multipliers(mult, G, rock)
   rmult = calculate_rock_multipliers(G, rock);

   assert (~isempty(rmult) && isnumeric(rmult) && ...
           (size(rmult, 1) == G.faces.num), 'Internal Logic Error');

   if isempty(mult)
      mult = rmult;

   elseif isnumeric(mult) && (size(mult, 1) == G.faces.num)
      mult = [ mult, rmult ];

   elseif isnumeric(mult) && (size(mult, 1) == size(G.cells.faces, 1))
      mult = struct('onesided', mult, 'face', rmult);

   elseif isstruct(mult)
      if ~isfield(mult, 'face')
         mult.face = rmult;

      elseif ~isfield(mult.face, 'direct')
         mult.face.direct = rmult;

      else
         mult.face.direct = [ mult.face.direct, rmult ];
      end

   else
      warning('Mult:UnknownType', ...
             ['Don''t know how to incorporate rock-based multipliers ', ...
              'into multiplier data of type ''%s''.'], class(mult));
   end
end

%--------------------------------------------------------------------------

function rmult = calculate_rock_multipliers(G, rock)
   rmult = [];

   if isfield(rock, 'multipliers')
      rmult = accumulateCartesianMultipliers(G, rock.multipliers);
   end

   if isfield(rock, 'faultdata')
      rmult = [ rmult, expand_fault_multipliers(G, rock.faultdata) ];
   end
end

%--------------------------------------------------------------------------

function multf = expand_fault_multipliers(G, faultdata)
	faults = processFaults(G, faultdata);

   numf  = vertcat(faults.numf);
   multf = accumarray(vertcat(faults.faces), ...
                      rldecode(vertcat(faults.mult), numf), ...
                      [G.faces.num, 1], @prod, 1);
end

%--------------------------------------------------------------------------

function T = apply_multipliers(T, mult, sname, map_fields)
% Supported types:
%
%   mult is numeric, shape numel(T)-by-n
%   mult is structure with field 'sname'
%      - mult.(sname) is numeric, shape numel(T)-by-n
%      - mult.(sname) is structure with one or both of the fields
%          - map_fields: Sparse representation => ACCUMARRAY
%          - direct: numeric, shape numel(T)-by-n

   if is_numeric_shape_match(T, mult)
      T = multiply_direct(T, mult);

   elseif isstruct(mult) && isfield(mult, sname)
      T = apply_substructure_multipliers(T, mult.(sname), map_fields);
   end
end

%--------------------------------------------------------------------------

function T = apply_substructure_multipliers(T, mult, map_fields)
   if is_numeric_shape_match(T, mult)
      T = multiply_direct(T, mult);

   elseif isstruct(mult)
      T = apply_structured_multipliers(T, mult, map_fields);

   else
      warning('Mult:UnknownType', ...
             ['Multipliers must be an appropriately shaped numeric ', ...
              'array or a structure. Got unusable ''%s''.'], class(mult));
   end
end

%--------------------------------------------------------------------------

function tf = is_numeric_shape_match(T, mult)
   tf = isnumeric(mult) && (size(mult, 1) == numel(T));
end

%--------------------------------------------------------------------------

function T = multiply_direct(T, mult)
   T = T .* prod(mult, 2);
end

%--------------------------------------------------------------------------

function T = apply_structured_multipliers(T, mult, map_fields)
   if all(isfield(mult, map_fields))
      T = T .* accumarray(mult.(map_fields{1})(:), ...
                          mult.(map_fields{2})(:), ...
                          [numel(T), 1], @prod, 1);
   end

   if isfield(mult, 'direct') && is_numeric_shape_match(T, mult.direct)
      T = multiply_direct(T, mult.direct);
   end
end

%--------------------------------------------------------------------------

function tf = has_cpgeometry(G)
   if isfield(G, 'nodes')
       physdim = size(G.nodes.coords, 2);
       numc    = G.cells.num;
       numcf   = size(G.cells.faces, 1);

       tf = isfield(G.cells, 'cpgeometry') ...
          && all(isfield(G.cells.cpgeometry, {'centroids', 'facecentroids'})) ...
          && (ndims(G.cells.cpgeometry.centroids) == 2) ...
          && (ndims(G.cells.cpgeometry.facecentroids) == 2) ...
          && all(size(G.cells.cpgeometry.centroids) == [numc , physdim]) ...
          && all(size(G.cells.cpgeometry.facecentroids) == [numcf, physdim]); %#ok<ISMAT>
   else
       tf = false;
   end
end

%--------------------------------------------------------------------------

function tf = is_key_type(keys)
   tf = iscellstr(keys);
end
