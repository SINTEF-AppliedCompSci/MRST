function deck = readGRID(fid, dirname, deck)

% deck = readGRID(fid, dirname, deck)

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

   if deck.RUNSPEC.DUALPORO
      [dims, nc, np, nv, ncdp] = get_dimensions(deck);
   else
      [dims, nc, np, nv]       = get_dimensions(deck);
   end

   [grd, miss_kw] = get_state(deck);

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section
      switch kw
         case 'BOX'
            boxKeyword(fid);

         case 'ENDBOX'
            endboxKeyword;

         case 'SPECGRID'
            s = removeQuotes(readRecordString(fid));
            cartDims = reshape(sscanf(s, '%f', 3), 1, []);

            if deck.RUNSPEC.DUALPORO % if dual porosity then half the z cartdimes
                cartDims(3) = cartDims(3) / 2;

                defaultBox(cartDims);
                gridBox(defaultBox);

                dims = reshape(cartDims, 1, []);
                nc   = prod(dims);                 % Number of cells
                np   = prod(dims(1 : end-1) + 1);  % Number of pillars
                nv   = prod(dims + 1);             % Number of vertices
                grd.cartDims = cartDims;
                ncdp = 2*nc;                       % Number of fracture and matric cells for dp properties
            else
                defaultBox(cartDims);
                gridBox(defaultBox);

                dims = reshape(cartDims, 1, []);
                nc   = prod(dims);                 % Number of cells
                np   = prod(dims(1 : end-1) + 1);  % Number of pillars
                nv   = prod(dims + 1);             % Number of vertices
                grd.cartDims = cartDims;
            end

         case 'GDORIENT'
            tmpl = {'INC', 'INC', 'INC', 'DOWN', 'RIGHT'};

            grd.(kw) = readDefaultedRecord(fid, tmpl);           clear tmpl

         case {'DXV', 'DYV', 'DZV'}
            ix       = strcmp(kw, {'DXV', 'DYV', 'DZV'});
            grd.(kw) = readVector(fid, kw, dims(ix));

         case 'DEPTHZ'
            grd.(kw) = readVector(fid, kw, np);

         case {'DX', 'DY', 'DZ', 'DEPTH', 'TOPS'}
            grd = readGridBoxArray(grd, fid, kw, nc, NaN);

         case {'COORDX', 'COORDY', 'COORDZ'}
            grd.(kw) = readVector(fid, kw, nv);

         case 'COORD'
            grd.(kw) = readVector(fid, kw, 6 * np);

         case 'ZCORN'
            grd = read_zcorn(grd, fid, dims, nc);

         case 'FAULTS'
            tmpl(1:8) = { 'Default' };
            data = readDefaultedKW(fid, tmpl);  clear tmpl
            data(:, 2:end-1) = to_double(data(:, 2:end-1));

            if ~isfield(grd, kw), grd.(kw) = cell([0, 8]); end
            grd.(kw) = [grd.(kw); data];

         case 'MAPAXES'
            s    = removeQuotes(readRecordString(fid));
            data = reshape(sscanf(s, '%f', 6), 1, []);

            assert (numel(data) == 6, ...
                   ['The ''MAPAXES'' keyword must be followed by ', ...
                    'exactly six data items.']);

            grd.(kw) = data;  clear s data

         case 'MAPUNITS'
            data     = readDefaultedRecord(fid, { 'METRES' });
            grd.(kw) = data{1};  clear data

         case 'MULTFLT'
            tmpl = { 'FaultName', '1.0', '1.0' };
            data = readDefaultedKW(fid, tmpl);  clear tmpl
            data(:, 2:end) = to_double(data(:, 2:end));

            if ~isfield(grd, kw), grd.(kw) = cell([0, 3]); end
            grd.(kw) = [grd.(kw); data];

         case 'NNC'
            tmpl = repmat({ 'NaN' }, [1, 7]);
            data = readDefaultedKW(fid, tmpl);
            data = to_double(data);
            grd.(kw) = reshape([data{:}], size(data));           clear tmpl

         case 'PINCH'
            tmpl         = { '1.0e-3', 'GAP', 'Inf', 'TOPBOT', 'TOP' };
            data         = readDefaultedRecord(fid, tmpl);
            data([1, 3]) = to_double(data([1, 3]));  clear tmpl
            grd.(kw)     = data;

         case 'JFUNC'
            tmpl         = { 'BOTH', 'NaN', 'NaN', '0.5', '0.5', 'XY'};
            data         = readDefaultedRecord(fid, tmpl);
            data(2:5) = to_double(data(2:5));  clear tmpl
            grd.(kw)     = data;

         case 'PINCHREG'
            nrec = 0;

            if isfield(grd, 'PINCHNUM') && ...
                  isfield(deck.RUNSPEC, 'GRIDOPTS')

               nrec = deck.RUNSPEC.GRIDOPTS{3};

            elseif isfield(grd, 'FLUXNUM')
               if isfield(deck.RUNSPEC, 'REGDIMS')
                  nrec = deck.RUNSPEC.REGDIMS{4};
               end
               if isfield(deck.RUNSPEC, 'TABDIMS')
                  nrec = max(nrec, deck.RUNSPEC.TABDIMS(11));
               end
            end

            if nrec > 0
               tmpl = { '1.0e-3', 'GAP', 'Inf', 'TOPBOT', 'TOP' };
               data = readDefaultedKW(fid, tmpl, 'NRec', nrec);

               data(:, [1, 3]) = to_double(data(:, [1, 3])); clear tmpl;

               grd.(kw) = data;
            end

         case 'ACTNUM'
            if deck.RUNSPEC.DUALPORO
               grd = readGridBoxArray(grd, fid, kw, ncdp);
            else
               grd = readGridBoxArray(grd, fid, kw, nc, 1);
            end

         case {'SIGMAV', 'SIGMADV', 'DZMTRXV'}
            if deck.RUNSPEC.DUALPORO
               grd = readGridBoxArray(grd, fid, kw, nc, 0.0);
            end

         case 'SIGMA'
            if deck.RUNSPEC.DUALPORO
               %Not handled yet
            end

         case {'PORO'  , ...
               'PERMX' , 'PERMXY', 'PERMXZ', ...
               'PERMYX', 'PERMY' , 'PERMYZ', ...
               'PERMZX', 'PERMZY', 'PERMZ' , ...
               'THCONR', ...
               }
            if deck.RUNSPEC.DUALPORO
               grd = readGridBoxArrayDP(grd, fid, kw, nc);
            else
               grd = readGridBoxArray(grd, fid, kw, nc, NaN);
            end

         case {'NTG'    , 'MULTPV' ,  ...
               'MULTX'  , 'MULTX-' ,  ...
               'MULTY'  , 'MULTY-' ,  ...
               'MULTZ'  , 'MULTZ-' ,  ...
               'FLUXNUM', 'MULTNUM', 'OPERNUM', 'PINCHNUM'}
            if deck.RUNSPEC.DUALPORO
               grd = readGridBoxArrayDP(grd, fid, kw, nc);
            else
               grd = readGridBoxArray(grd, fid, kw, nc, 1.0);
            end

         case {'ADD', 'COPY', 'EQUALS', 'MAXVALUE', ...
               'MINVALUE', 'MULTIPLY'}
            grd = applyOperator(grd, fid, kw);

         case 'ADDZCORN'
            op  = read_zcorn_operator(fid, 'ALL');
            grd = add_zcorn(grd, op, dims);

         case 'EQLZCORN'
            op  = read_zcorn_operator(fid, 'TOP');
            grd = assign_zcorn(grd, op, dims);

         case 'MINPV'
            tmpl = { '1.0e-6' };  % In all unit systems
            data = readDefaultedRecord(fid, tmpl);
            grd.(kw) = to_double(data{1});

         case 'MINPVV'
            if deck.RUNSPEC.DUALPORO
               grd = readGridBoxArrayDP(grd, fid, kw, nc);
            else
               grd = readGridBoxArray(grd, fid, kw, nc, 1.0e-6);
            end

         case {'ECHO', 'NOECHO'}
            kw = getEclipseKeyword(fid);
            continue;  % Ignore.  Not handled in MRST

         %-----------------------------------------------------------------
         % Sectioning keywords below.  Modifies flow of control.
         % Don't change unless absolutely required...
         %
         case 'END'
            % Logical end of input deck.
            % Quite unusual (but nevertheless legal) in GRID.
            %
            in_section = false;
            deck.GRID  = grd;

            % Restore default input box at end of section
            gridBox(defaultBox);

         case 'EDIT'
            % Read next section (i.e., 'GRID' -> 'EDIT' -> 'PROPS')
            in_section = false;

            deck = set_state(deck, grd, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readEDIT(fid, dirname, deck);
            grd = deck.GRID;

         case 'PROPS'
            % Read next section (i.e., 'GRID' -> 'PROPS', no 'EDIT')
            in_section = false;

            deck = set_state(deck, grd, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readPROPS(fid, dirname, deck);

         case 'INCLUDE'
            % Handle 'INCLUDE' (recursion).
            deck = set_state(deck, grd, miss_kw);

            deck = readEclipseIncludeFile(@readGRID, fid, dirname, ...
                                          deck.RUNSPEC, deck);

            [grd, miss_kw] = get_state(deck);

         otherwise
            if ischar(kw)
               miss_kw = [ miss_kw, { kw } ];  %#ok
            end
      end

      % Get next keyword.
      kw = getEclipseKeyword(fid);
      in_section = in_section && ischar(kw);
   end

   deck = set_state(deck, grd, miss_kw);
end

%--------------------------------------------------------------------------

function [dims, nc, np, nv,ncdp] = get_dimensions(deck)
   assert (isstruct(deck)         && isfield(deck, 'RUNSPEC') && ...
           isstruct(deck.RUNSPEC) && isfield(deck.RUNSPEC, 'cartDims'));

   dims = reshape(deck.RUNSPEC.cartDims, 1, []);
   nc   = prod(dims);                 % Number of cells
   np   = prod(dims(1 : end-1) + 1);  % Number of pillars
   nv   = prod(dims + 1);             % Number of vertices
   ncdp = 2*nc;                         %Number of entries for DP model 2*nc
end

%--------------------------------------------------------------------------

function v = to_double(v)
   convert = @(s) sscanf(regexprep(s, '[dD]', 'e'), '%f');

   if ischar(v)
      v = convert(v);
   else
      v = cellfun(convert, v, 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function [grd, miss_kw] = get_state(deck)
   grd     = deck.GRID;
   miss_kw = deck.UnhandledKeywords.GRID;
end

%--------------------------------------------------------------------------

function deck = set_state(deck, grd, miss_kw)
   deck.GRID                   = grd;
   deck.UnhandledKeywords.GRID = unique(miss_kw);
end

%--------------------------------------------------------------------------

function grd = read_zcorn(grd, fid, dims, nc)
   i  = gridBox;
   kw = 'ZCORN';

   if all(i == defaultBox)
      % Input box covers entire grid.  This is the common case, so use
      % simple (and efficient) implementation here.

      grd.(kw) = readVector(fid, kw, 8 * nc);
   else
      % Input box covers a subset of the grid.  Need to take current box
      % limit into account when deciding how much data to read and where to
      % put it in the overall 'ZCORN' array.
      %
      % Note: This code is very similar to the implementation of function
      % 'readGridBoxArray', but *that* function assumes that the array has
      % a single scalar value for each (global) cell.  We may want to
      % refactor this code.

      if ~isfield(grd, kw), grd.(kw) = zeros([8 * nc, 1]); end

      [j{1 : 3}] = ndgrid((2*i(1) - 1) : (2 * i(2)), ...
                          (2*i(3) - 1) : (2 * i(4)), ...
                          (2*i(5) - 1) : (2 * i(6)));

      ix = sub2ind(2 .* dims           , ...
                   reshape(j{1}, [], 1), ...
                   reshape(j{2}, [], 1), ...
                   reshape(j{3}, [], 1));

      grd.(kw)(ix) = readVector(fid, kw, numel(ix));
   end
end

%--------------------------------------------------------------------------

function op = read_zcorn_operator(fid, dflt_action)
   tmpl = [{ 'NaN' }                                          , ... % const
           arrayfun(@int2str, gridBox, 'UniformOutput', false), ... % box
           repmat({ 'NaN' }, [1, 4])                          , ... % cont
           { dflt_action } ];                                       % action

   op  = {};
   rec = readDefaultedRecord(fid, tmpl);

   while ~isequal(rec, tmpl)
      op = [ op; rec ];                       %#ok % Willfully ignore MLINT

      tmpl(2 : 7) = rec(2 : 7);

      rec = readDefaultedRecord(fid, tmpl);
   end
end

%--------------------------------------------------------------------------

function grd = add_zcorn(grd, addzcorn, dims)
   assert (isfield(grd, 'ZCORN'), ...
           'Attempt to modify non-existent ''ZCORN'' array.');

   assert (size(addzcorn, 2) == 12, ...
           '''ZCORN'' operator must have twelve items per row.');

   for r = 1 : size(addzcorn, 1)
      grd.ZCORN = add_zcorn_action(grd.ZCORN, addzcorn(r, :), dims);
   end
end

%--------------------------------------------------------------------------

function grd = assign_zcorn(grd, eqlzcorn, dims)
   assert (isfield(grd, 'ZCORN'), ...
           'Attempt to modify non-existent ''ZCORN'' array.');

   assert (size(eqlzcorn, 2) == 12, ...
           '''ZCORN'' operator must have twelve items per row.');

   for r = 1 : size(addzcorn, 1)
      grd.ZCORN = assign_zcorn_action(grd.ZCORN, eqlzcorn(r, :), dims);
   end
end

%--------------------------------------------------------------------------

function z = add_zcorn_action(z, a, dims)
   [c, b, cont, affects] = zcorn_operator_components(a);

   [i, j, k] = zcorn_box(b, affects);

   if ~any(cellfun('isempty', {i, j, k}))
      [I, J, K] = ndgrid(i, j, k);
      ix_cont   = zcorn_continuity(dims, cont, b, i, j, k, z);

      ix    = [sub2ind(2 .* dims, I(:), J(:), K(:)); ix_cont];
      z(ix) = z(ix) + c;
   end
end

%--------------------------------------------------------------------------

function z = assign_zcorn_action(z, e, dims)
   [c, b, cont, affects] = zcorn_operator_components(e);

   nlayers = b(6) - b(5);
   if nlayers > 1
      warning('EQLZCORN:MultiLayers', ...
             ['''EQLZCORN'' operator box affects %d layers. Should ', ...
              'normally be a single layer only.'], nlayers);
   end

   [i, j, k] = zcorn_box(b, affects);

   if ~any(cellfun('isempty', {i, j, k}))
      [I, J, K] = ndgrid(i, j, k);
      ix_cont   = zcorn_continuity(dims, cont, b, i, j, k, z);

      ix    = [sub2ind(2 .* dims, I(:), J(:), K(:)); ix_cont];
      z(ix) = c;
   end
end

%--------------------------------------------------------------------------

function [c, b, cont, affects] = zcorn_operator_components(op)
   d       = to_double(op(1 : 11));
   affects = op{ end };

   c       = d{ 1 };          % Constant
   b       = [ d{2 :  7} ];   % Grid box
   cont    = [ d{8 : 11} ];   % Continuity specifiers
end

%--------------------------------------------------------------------------

function [i, j, k] = zcorn_box(b, affects)
   assert (~any(b < 0), ...
           'ZCORN bounding box must have non-negative limits.');

   lb = b([1, 3]);
   ub = b([2, 4]);

   assert ((all(lb > 0) && all(ub > 0) && all(lb <= ub)) || ...
           all(xor(lb > 0, ub > 0)), ...
          ['ZCORN lateral bounding box must identify at least ', ...
           'one corner or a non-empty cell range.']);

   assert ((b(5) > 0) && (b(6) > 0) && (b(5) <= b(6)), ...
          ['ZCORN vertical bounding box must identify a non-empty ', ...
           'cell range.']);

   i = index_range(lb(1), ub(1));
   j = index_range(lb(2), ub(2));
   k = index_range(b (5), b (6));

   assert (numel(k) >= 2, 'Internal error.');

   switch affects(1)
      case 'A'  % Operator affects ALL corners
         % Trivial.  Nothing to do.

      case 'B'  % Operator affects BOTTOM corners (TOP in b(5) excluded)
         k = k(2 : end);

      case 'T'  % Operator affects TOP corners (BOTTOM in b(6) excluded)
         k = k(1 : end - 1);

      otherwise
         warning('readGRID:ZcornBox:Unsupported', ...
                 'ZCORN action ''%s'' unsupported. Ignored.', affects);

         [i, j, k] = deal([]);
   end
end

%--------------------------------------------------------------------------

function ix = zcorn_continuity(dims, cont, b, i, j, k, z)
   ix = [];

   p      = isnan(cont);     % Defaulted continuity specifiers
   adjust = [-(b(1) > 1),       ...  % Left  if box starts at > 1
              (b(2) < dims(1)), ...  % Right if box ends   at < DIMENS(1)
             -(b(3) > 1),       ...  % Back  if box starts at > 1
              (b(4) < dims(2))];     % Front if box ends   at < DIMENS(2)

   cont(p) = b(p) + adjust(p);

   % Continuous surface modifaction to the left of box, exclude faults
   if cont(1) == b(1) - 1
      args = {{ 2*cont(1), j, k }, { i(1), j, k }};
      ix   = [ix ; box_boundary(dims, z, 1, args)];
   end

   % Continuous surface modifaction to the right of box, exclude faults
   if cont(2) == b(2) + 1
      args = {{ i(end), j, k }, { 2*cont(2) - 1, j, k }};
      ix   = [ix ; box_boundary(dims, z, 2, args)];
   end

   % Continuous surface modifaction to the rear of box, exclude faults
   if cont(3) == b(3) - 1
      args = {{ i, 2*cont(3), k }, { i, j(1), k }};
      ix   = [ix ; box_boundary(dims, z, 1, args)];
   end

   % Continuous surface modifaction to the front of box, exclude faults
   if cont(3) == b(4) + 1
      args = {{ i, j(end), k }, { i, 2*cont(4) - 1, k }};
      ix   = [ix ; box_boundary(dims, z, 2, args)];
   end
end

%--------------------------------------------------------------------------

function i = index_range(lb, ub)
   if lb == 0

      i = 2*ub - 0;

   elseif ub == 0

      i = 2*lb - 1;

   else

      assert (lb <= ub, 'Index limits must define non-empty range.');

      i = (2*lb - 1) : 2*ub;

   end
end

%--------------------------------------------------------------------------

function i = box_boundary(dims, z, col, args)
   [lb{1:3}] = ndgrid(args{1}{:});
   [ub{1:3}] = ndgrid(args{2}{:});

   ix = @(b) sub2ind(2 .* dims, b{1}(:), b{2}(:), b{3}(:));

   i  = [ ix(lb), ix(ub) ];
   p  = z(i(:,1)) == z(i(:,2));   % Exclude points on faults
   i  = i(p, col);
end
