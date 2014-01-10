function writeDeck(deck, dirname, varargin)
%Dump parts of an Eclipse deck to individual files.
%
% SYNOPSIS:
%   dumpdeck(filename, dirname)
%
% DESCRIPTION:
%   This function dumps the contents of selected keywords to individual
%   fi_les in a directory specified by the user. This is used as a
%   preliminary interface to the C/C++ simulator code developed for URC.
%
%   Currently, the following keywords are written to file
%
%   In the RUNSPEC section:
%        DIMENS
%
%   In the GRID section:
%        -- Eclipse cornerpoint
%        COORD
%        ZCORN
%
%        ACTNUM
%        PERMX
%        PERMY
%        PERMZ
%        PORO
%
%   In the PROPS section
%        SGOF
%        PVDO
%        PVDG
%        DENSITY
%
%   In the SOLUTION section
%        PRESSURE
%        SGAS
%        SOIL
%
%   In the SCHEDULE section
%        WELSPECS
%        COMPDAT
%        WCONINJE
%        WCONPROD
%        TSTEP
%
% REQUIRED PARAMETERS:
%   filename - name of Eclipse input file.
%
%   dirname  - name of directory in which to dump data.  If the directory
%              does not exist, it will be created.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   opt = struct('grid', [], 'rock', []);
   opt = merge_options(opt, varargin{:});

   if isempty(opt.grid) || isempty(opt.rock)
      require deckformat
   end

   if ~isdir(dirname),
      mkdir(dirname);
   end

   [fn, fn] = fileparts(dirname);                                      %#ok
   fn       = fullfile(dirname, [fn, '.DATA']);

   [fid, msg] = fopen(fn, 'w+');
   if fid < 0,
      error(msgid('Open:Failure'), 'Failed to open ''%s'': %s\n', fn, msg);
   end

   if isempty(opt.grid),
      G = initEclipseGrid(deck);
   else
      G = opt.grid;
      if ~isfield(G.cells, 'indexMap'),
         G.cells.indexMap = (1 : G.cells.num) .';
      end
   end

   if ~isfield(G, 'griddim'), G.griddim = 3; end

   if isempty(opt.rock),
      rock = initEclipseRock(deck);
   else
      rock = opt.rock;
   end

   fprintf(fid, '%s\n', '-- Generate deck from MRST writeDeck');
   fprintf(fid,'%s\n','RUNSPEC');
   if(isfield(deck.RUNSPEC,'TITLE'))
      fprintf(fid, '%s\n', 'TITLE');
      fprintf(fid, '%s\n', deck.RUNSPEC.TITLE);
      fprintf(fid, '\n');
   end
   fprintf(fid, '\n');
   fprintf(fid, '------------------------------------------------------\n');
   dump_runspec(fid,dirname, deck);
   fprintf(fid, '------------------------------------------------------\n');
   fprintf(fid, '%s\n', 'GRID');
   dump_grid(fid,dirname, deck, rock);
   fprintf(fid, '------------------------------------------------------\n');
   fprintf(fid, '%s\n', 'PROPS');
   dump_props(fid, dirname, deck);
   fprintf(fid, '------------------------------------------------------\n');
   fprintf(fid, '%s\n', 'SOLUTION');
   dump_solution(fid,dirname, deck);
   fprintf(fid, '------------------------------------------------------\n');
   fprintf(fid, '%s\n', 'SUMMARY');
   % MRST do not handle SUM
   %%{
   %myfields=fieldnames(deck.SUMMARY)
   if(isfield(deck,'UnhandledKeywords'))
      if(isfield(deck.UnhandledKeywords,'SUMMARY'))
         myfields=deck.UnhandledKeywords.SUMMARY;
         fprintf(fid,'\n');
         for i=1:numel(myfields)
            fprintf(fid,'%s\n/\n\n',myfields{i});
         end
      end
   end
   %%}
   fprintf(fid, '------------------------------------------------------\n');
   % to do
   rock = compressRock(rock, G.cells.indexMap);
   fprintf(fid, '%s\n', 'SCHEDULE');
   dump_schedule(fid,dirname, deck, G, rock);
   fclose(fid);
end

function dump_runspec(fid,dirname, deck)
   myfields={'DIMENS','EQLDIMS','TABDIMS','WELLDIMS'};%,'NUPCOL'}
   for i=1:numel(myfields)
      myfield=myfields{i};
      if(isfield(deck.RUNSPEC,myfield))
         dump_vector(fid,dirname, lower(myfield), '%9i\n', deck.RUNSPEC.(myfield));
      end
   end
   myfields={'OIL','WATER','GAS','METRIC','NOGRAV','FIELD'};
   for i=1:numel(myfields)
      myfield=myfields{i};
      if(isfield(deck.RUNSPEC,myfield))
         if(deck.RUNSPEC.(myfield)==1)
            fprintf(fid,'%s \n',myfield);
         end
      end
   end
   fprintf(fid,'\n');
   fprintf(fid,'START \n');
   dato=datestr(deck.RUNSPEC.START);
   dato=regexp(dato,'-','split');
   fprintf(fid, ' %s ''%s'' %s\n/\n', dato{1}, upper(dato{2}), dato{3});
end

%--------------------------------------------------------------------------

function dump_grid(fid,dirname, deck, rock)
   grid_support = true;

   if isfield(deck.GRID, 'COORD'),

      assert (isfield(deck.GRID, 'ZCORN'));
      coord = deck.GRID.COORD;
      zcorn = deck.GRID.ZCORN;

   elseif all(isfield(deck.GRID, {'DXV', 'DYV', 'DZV'}))

      [coord, zcorn] = block_centred_to_cpg(deck);

   else

      grid_support = false;

   end

   if grid_support,
      fprintf(fid, 'SPECGRID\n');
      fprintf(fid, '%d ', deck.GRID.cartDims);
      fprintf(fid, '1 F\n/\n\n');

      dump_vector(fid, dirname, 'coord', '%18.16e\n', coord);
      dump_vector(fid, dirname, 'zcorn', '%18.16e\n', zcorn);

      if isfield(deck.GRID, 'ACTNUM'),
         dump_vector(fid, dirname, 'actnum', '%d\n', deck.GRID.ACTNUM);
      end
   end

   K = permTensor(rock, 3);
   dump_vector(fid, dirname, 'permx', '%18.16e\n', K(:,1));
   dump_vector(fid, dirname, 'permy', '%18.16e\n', K(:,5));
   dump_vector(fid, dirname, 'permz', '%18.16e\n', K(:,9));
   dump_vector(fid, dirname, 'poro' , '%18.16e\n', rock.poro);
end

%--------------------------------------------------------------------------

function [coord, zcorn] = block_centred_to_cpg(deck)
   dims = deck.RUNSPEC.DIMENS;
   n    = prod(dims(1:2) + 1);

   x = cumsum([0; deck.GRID.DXV]);
   y = cumsum([0; deck.GRID.DYV]);
   z = cumsum([0; deck.GRID.DZV]);

   [X, Y, Z] = ndgrid(x, y, z);

   if isfield(deck.GRID, 'DEPTHZ'),
      Z = bsxfun(@plus, Z, reshape(deck.GRID.DEPTHZ, dims(1:2) + 1));
   end

   lines = zeros([n, 6]);
   lines(:, [1, 4]) = reshape(X(:,:,[1, end]), [n, 2]);
   lines(:, [2, 5]) = reshape(Y(:,:,[1, end]), [n, 2]);
   lines(:, [3, 6]) = reshape(Z(:,:,[1, end]), [n, 2]);
   coord = reshape(lines.', [], 1);

   % Assign z-coordinates
   % ind(d) == [1, 2, 2, 3, 3, ..., dims(d), dims(d), dims(d)+1]
   ind = @(d) 1 + fix((1 : 2*dims(d)) ./ 2);
   zcorn = reshape(Z(ind(1), ind(2), ind(3)), [], 1);
end

%--------------------------------------------------------------------------

function dump_props(fid, dirname, deck)
   for fld = reshape(fieldnames(deck.PROPS), 1, []),
      values = deck.PROPS.(fld{1});

      if strcmp(fld{1}, 'PVTO'),
         % Custom output routine
         dump_pvto(fid, dirname, values{1});
         continue
      end

      if iscell(values),
         values = values{1};
      end

      values(isnan(values)) = 0;
      dump_vector(fid, dirname, lower(fld{1}), '%18.16e\n', values');
   end
end

%--------------------------------------------------------------------------

function dump_pvto(fid, dirname, pvto)
   fprintf(fid, 'INCLUDE\npvto.txt /\n\n');

   [fid_pvto, msg] = fopen(fullfile(dirname, 'pvto.txt'), 'wt');
   if fid_pvto < 0,
      error('Failed to open ''pvto'' output file: %s', msg);
   end

   assert (numel(pvto.key) + 1 == numel(pvto.pos));

   fprintf(fid_pvto, 'PVTO\n');

   for r = 1 : numel(pvto.pos) - 1,
      fprintf(fid_pvto, '%.10e\n', pvto.key(r));

      i = pvto.pos(r) : pvto.pos(r + 1) - 1;
      fprintf(fid_pvto, '%.10e %.10e %.10e\n', pvto.data(i,:).');
      fprintf(fid_pvto, '/\n');
   end

   fprintf(fid_pvto, '/\n');
   fclose(fid_pvto);
end

%--------------------------------------------------------------------------

function dump_solution(fid,dirname, deck)
   myfields=fieldnames(deck.SOLUTION);
   for i=1:numel(myfields)
      myfield=myfields{i};

      v = deck.SOLUTION.(myfield);
      if iscell(v), v = v{1}; end

      dump_vector(fid,dirname, lower(myfield), '%18.16e\n', v);
   end
end

%--------------------------------------------------------------------------

function dump_schedule(fid,dirname, deck, G, rock)
   org_fid=fid;

   if isfield(deck.SCHEDULE, 'control'),
   fprintf(org_fid,'%s\n','INCLUDE');
   fprintf(org_fid,'%s\n','welspecs.txt');
   fprintf(org_fid,'/\n\n');
   fid = fopen(fullfile(dirname, 'welspecs.txt'), 'wt');
   wspecs = replace_default(deck.SCHEDULE.control(1).WELSPECS).';
   fprintf(fid,'%s\n',upper('welspecs'));
   fprintf(fid, '%s %s %d %d %18.16e %s %d %s %s %s %d %s %d/\n', wspecs{:});
   fprintf(fid, '/\n\n');
   fclose(fid);

   fprintf(org_fid,'%s\n','INCLUDE');
   fprintf(org_fid,'%s\n','compdat.txt');
   fprintf(org_fid,'/\n\n');
   fid = fopen(fullfile(dirname, 'compdat.txt'), 'wt');
   cdat = replace_default(set_WI(G, rock, deck.SCHEDULE.control(1).COMPDAT)) .';
   %cdat = replace_default(deck.SCHEDULE.control(1).COMPDAT);
   fprintf(fid,'%s\n',upper('compdat'));
   fprintf(fid, '%s %d %d %d %d %s %d %18.16e %18.16e %18.16e %d %s %s %18.16e/\n', cdat{:});
   fprintf(fid, '/\n\n');
   fclose(fid);

   fprintf(org_fid,'%s\n','INCLUDE');
   fprintf(org_fid,'%s\n','wconinje.txt');
   fprintf(org_fid,'/\n\n');
   fid = fopen(fullfile(dirname, 'wconinje.txt'), 'wt');
   wconinje = replace_default(deck.SCHEDULE.control(1).WCONINJE).';
   fprintf(fid,'%s\n',upper('wconinje'));
   s = sprintf(['%s %s %s %s %18.16e %18.16e ', ...
                '%18.16e %18.16e %d %18.16e ' , ...
                '%18.16e %18.16e %18.16e %18.16e /\n'], wconinje{:});
   fprintf(fid, '%s', regexprep(s, 'Inf|NaN', '1*', 'ignorecase'));
   fprintf(fid, '/\n\n');
   fclose(fid);

   fprintf(org_fid,'%s\n','INCLUDE');
   fprintf(org_fid,'%s\n','wconprod.txt');
   fprintf(org_fid,'/\n\n');
   fid = fopen(fullfile(dirname, 'wconprod.txt'), 'wt');
   wconprod = replace_default(deck.SCHEDULE.control(1).WCONPROD).';
   fprintf(fid,'%s\n',upper('wconprod'));
   s = sprintf('%s %s %s %18.16e %18.16e %18.16e %18.16e %18.16e %18.16e %d %d %d/\n', wconprod{:});
   fprintf(fid, '%s', regexprep(s, 'Inf|NaN', '1*', 'ignorecase'));
   fprintf(fid, '/\n\n');
   fclose(fid);
   end

   dump_vector(org_fid,dirname, 'tstep', '%18.16e\n', deck.SCHEDULE.step.val);
end

%--------------------------------------------------------------------------

function dump_vector(fid, dirname, field, fmt, v)
   org_fid = fid;

   fn         = fullfile(dirname, [field, '.txt']);
   [fid, msg] = fopen(fn, 'wt');
   if fid < 0,
      error('Unable to open ''%s'' for writing: %s', fn, msg);
   end

   fprintf(fid, '%s\n', upper(field));
   fprintf(fid, fmt, v);
   fprintf(fid, '/\n');

   fclose(fid);

   fprintf(org_fid, '%s\n', 'INCLUDE');
   fprintf(org_fid, '%s\n', [field, '.txt']);
   fprintf(org_fid, '/\n\n');
end

%--------------------------------------------------------------------------

function compdat = set_WI(G, rock, compdat)
   ijk = reshape([compdat{:,2:5}], [], 4);
   n   = 1 + (ijk(:,4) - ijk(:,3));

   %compdat = rldecode(compdat, n);

   WI   = vertcat(compdat{:,8});
   diam = vertcat(compdat{:, 9});
   diam(~(diam > 0)) = 1*ft;

   diam = abs(diam);

   if any(WI < 0),
      k   = mcolon(ijk(:,3), ijk(:,4)).';
      ijk = [rldecode(ijk(:,1:2), n), k];

      c = cart2active(G, sub2ind(G.cartDims, ijk(:,1), ijk(:,2), ijk(:,3)));

      K = permTensor(rock, G.griddim);
      K = K(c,:);

      k1 = K(:,1);
      k2 = K(:, 1 + G.griddim + 1);
      [d1, d2, ell] = cellDims(G, c);

      wc = 0.14;

      re1 = 2 * wc .* sqrt((d1.^2).*sqrt(k2 ./ k1) + ...
                        (d2.^2).*sqrt(k1 ./ k2));
      re2 = (k2 ./ k1).^(1/4) + (k1 ./ k2).^(1/4);

      re  = reshape(re1 ./ re2, [], 1);
      ke  = sqrt(k1 .* k2);

      Skin = vertcat(compdat{:,11});
      Skin(Skin < 0) = 0.0;

      Kh    = vertcat(compdat{:,10}); i = Kh < 0;
      Kh(i) = ell(i) .* ke(i);

      radius = diam ./ 2;
      WI2    = 2 * pi * Kh ./ (log(re ./ radius) + Skin);

      WI(WI < 0) = WI2(WI < 0);
   end

   compdat(:,8) = num2cell(WI);
   compdat(:,9) = num2cell(diam);
end

%--------------------------------------------------------------------------

function t = replace_default(t)
   i = cellfun(@ischar, t(1,:));
   t(:,i) = regexprep(t(:,i), 'Default', '1*');
end

%--------------------------------------------------------------------------

function [dx, dy, dz] = cellDims(G, ix)
% cellDims -- Compute physical dimensions of all cells in single well
%
% SYNOPSIS:
%   [dx, dy, dz] = cellDims(G, ix)
%
% PARAMETERS:
%   G  - Grid data structure.
%   ix - Cells for which to compute the physical dimensions (bounding
%        boxes).
%
% RETURNS:
%   dx, dy, dz -- Size of bounding box for each cell.  In particular,
%                 [dx(k),dy(k),dz(k)] is Cartesian BB for cell ix(k).

n = numel(ix);
[dx, dy, dz] = deal(zeros([n, 1]));

ixc = G.cells.facePos;
ixf = G.faces.nodePos;

for k = 1 : n,
   c = ix(k);                                     % Current cell
   f = G.cells.faces(ixc(c) : ixc(c + 1) - 1, 1); % Faces on cell
   e = mcolon(ixf(f), ixf(f + 1) - 1);            % Edges on cell

   nodes  = unique(G.faces.nodes(e, 1));          % Unique nodes...
   coords = G.nodes.coords(nodes,:);              % ... and coordinates

   % Compute bounding box
   m = min(coords);
   M = max(coords);

   % Size of bounding box
   dx(k) = M(1) - m(1);
   if size(G.nodes.coords, 2) > 1,
      dy(k) = M(2) - m(2);
   else
      dy(k) = 1;
   end

   if size(G.nodes.coords, 2) > 2,
      dz(k) = M(3) - m(3);
   else
      dz(k) = 0;
   end
end
end
