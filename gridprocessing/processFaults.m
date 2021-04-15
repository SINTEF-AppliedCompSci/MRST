function [flt, flt_id] = processFaults(G, geomspec)
%Construct fault structure from input specification (keyword `FAULTS`)
%
% SYNOPSIS:
%   [faults, id] = processFaults(G, geom)
%
% PARAMETERS:
%   G    - Valid grid data structure.  Must contain the second (third)
%          column of G.cells.faces, the values of which must identify the
%          cardinal direction of each face within each cell.
%
%   geom - Geometry specification.  Typically corresponds to the `GRID`
%          section data structure defined by function `readEclipseDeck`
%          (i.e., `deck.GRID`).
%
%          If `geom` contains the fault-related keywords `FAULTS` and,
%          optionally, `MULTFLT` (matched case insensitively), then fault
%          structures will be generate for each named fault.  Otherwise,
%          empty return values will be produced by function
%          `processFaults`.
%
% RETURNS:
%   faults - An n-by-1 structure array, one element for each of the `n`
%            unique faults in the geometry specification.  Each array
%            element has the following fields:
%
%               name  - Fault name.  Copied from `geom`.
%
%               faces - Global grid faces (from the grid `G`) connected to
%               the given, named fault.
%
%               numf  - Number of global faces connected to given, named
%               fault (`== NUMEL(faces)`).
%
%               mult  - Fault transmissibility multiplier.
%               Numeric value 1 (one) unless redefined in `MULT`.
%
%   id     - Fault enumeration mapping.  Function handle supporting the
%            following syntax::
%
%                     i = id(name)
%
%            where `name` is one of the unique fault names defined in the
%            `FAULTS` field of `geom`.  The return value, `i`, is the index
%            into `faults` such that `faults(i)` contains the information
%            pertaining to the `name`d fault.
%
% SEE ALSO:
%   `readEclipseDeck`

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

   [faultspec, multflt] = fault_keywords(geomspec);

   if isempty(faultspec)
      [flt, flt_id] = deal([]);
   else
      if size(G.cells.faces, 2) < 2
         error(['Function ''%s'' is only supported in corner-point ', ...
                'grids which identify cardinal directions.'], mfilename);
      end

      [flt, flt_id] = build_faults(G, faultspec);

      if ~isempty(multflt)
         i             = cellfun(flt_id, multflt(:, 1));
         [flt(i).mult] = multflt{:, 2};
      end

      i = isfinite([ flt.mult ]) & ([ flt.numf ] > 0);
      if any(~i)
         % Remove any faults which have not been assigned multipliers or
         % which do not contain any faces.
         %
         remap    = zeros(size(i));
         remap(i) = 1 : sum(i);
         flt_id   = @(s) remap(flt_id(s));
         flt      = flt(i);
      end
   end
end

%--------------------------------------------------------------------------

function [spec, mult] = fault_keywords(geom)
   exists = @(p,a) ~cellfun(@isempty, regexpi(a, p, 'match', 'once'));
   fields = fieldnames(geom);

   i = exists('^faults$', fields);
   if any(i)
      assert (sum(double(i)) == 1);
      spec = geom.(fields{i});
   else
      spec = [];
   end

   j = exists('^multflt$', fields);
   if any(j)
      assert (sum(double(j)) == 1);
      mult = geom.(fields{j});
   else
      mult = [];
   end
end

%--------------------------------------------------------------------------

function [flt, flt_id] = build_faults(G, faultspec)
   [uflt, flt_id, p, rows] = enumerate_faults(faultspec);
   dir_id                  = identify_directions();

   act                   = zeros([prod(G.cartDims), 1]);
   act(G.cells.indexMap) = 1 : G.cells.num;

   flt(1 : numel(uflt))  = struct('name' , uflt, ...
                                  'faces', []  , ...
                                  'numf' , 0   , ...
                                  'mult' , NaN );

   on_fault = false([G.faces.num, 1]);

   for i = 1 : numel(uflt)
      % Fault spec rows defining this individual fault.
      %
      r = rows(p(i) : p(i+1) - 1);

      % Determine cells affected by this fault.
      %
      c = arrayfun(@(m) fault_cells(G, faultspec{m, 2:end-1}), r, ...
                   'UniformOutput', false);

      % Extract the global grid faces
      f = cellfun(@(c, d) stack_faces(G, act(c), dir_id(d)), ...
                  c, faultspec(r, end), 'UniformOutput', false);

      f = vertcat(f{:});
      on_fault(f) = true;

      % Assign faces to current fault.
      %
      flt(i).faces = find(on_fault);      on_fault(f) = false;
      flt(i).numf  = numel(flt(i).faces);
   end
end

%--------------------------------------------------------------------------

function [u, id, p, rows] = enumerate_faults(spec)
   u = unique(spec(:,1));

   try
      % Java's O(1) hash table string search support.
      ht = java.util.Hashtable;
      for n = 1 : numel(u)
         ht.put(u{n}, n);
      end
      id = @(s) ht.get(s);
   catch %#ok
      % Fall back to (probably) linear structure field name search if Java
      % is unavailable.
      ht = struct();
      for n = 1 : numel(u)
         ht.(regexprep(u{n}, '\W', '_')) = n;
      end
      id = @(s) ht.(regexprep(s, '\W', '_'));
   end

   a = sortrows([cellfun(id, spec(:,1)), (1 : size(spec, 1)) .']);
   p = cumsum([1; accumarray(a(:,1), 1, [numel(u), 1])]);
   rows = a(:,2);
end

%--------------------------------------------------------------------------

function id = identify_directions()
   s  = struct('x_', 1, 'i_', 1, 'x', 2, 'i', 2, ...
               'y_', 3, 'j_', 3, 'y', 4, 'j', 4, ...
               'z_', 5, 'k_', 5, 'z', 6, 'k', 6);
   id = @(d) s.(lower(regexprep(d, {'+', '-'}, {'', '_'})));
end

%--------------------------------------------------------------------------

function c = fault_cells(G, i1, i2, j1, j2, k1, k2)
   [i{1:3}] = ndgrid(i1 : i2, j1 : j2, k1 : k2);

   c = sub2ind(reshape(G.cartDims,  1, []), ...
               reshape(i{1}      , [],  1), ...
               reshape(i{2}      , [],  1), ...
               reshape(i{3}      , [],  1));
end

%--------------------------------------------------------------------------

function f = stack_faces(G, c, cftag)
   c = c(c > 0);

   if isempty(c)
      f = [];
   else
      ix = mcolon(G.cells.facePos(c), G.cells.facePos(c + 1) - 1) .';
      i  = G.cells.faces(ix   , 2) == cftag;
      f  = G.cells.faces(ix(i), 1);
   end
end
