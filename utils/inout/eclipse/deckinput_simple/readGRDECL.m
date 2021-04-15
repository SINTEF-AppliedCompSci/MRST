function [grdecl, unrec] = readGRDECL(fn, varargin)
%Read subset of ECLIPSE GRID file
%
% SYNOPSIS:
%   grdecl         = readGRDECL(fn)
%   grdecl         = readGRDECL(fn, 'pn1', pv1, ...)
%   [grdecl,unrec] = readGRDECL(fn, 'pn1', pv1, ...)
%
% PARAMETERS:
%   fn      - String holding name of (readable) GRDECL specification.
%             Note: The GRDECL specification is assumed to be a physical
%             (seekable) file on disk, not to be read in through (e.g.) a
%             POSIX pipe.
%
% OPTIONAL PARAMETERS:
%   'verbose'        - Emit messages to screen while processing.
%                      Default value: FALSE.
%
%   'keywords'       - Cell array of strings that should be recognized
%                      outside the predefined set of keywords. NB! Are
%                      interpreted as one value per cell. Default value: {}
%
%   'missing_include - Callback function through which to handle missing
%                      `INCLUDE` files.  Must support the syntax ::
%
%                          missing_include(id, '%s', msg)
%
%                      with 'id' being a message ID and 'msg' being a
%                      string (diagnostic message).
%
%                      Default value: missing_include = @error (end input
%                      reading with an error/failure if we encounter a
%                      missing `INCLUDE` file).
%
%                      Other possible values are @warning (report warning),
%                      @(varargin) [] (ignore everything) and similar.
%
% RETURNS:
%   grdecl - Output structure containing the known, though mostly
%            unprocessed, fields of the GRDECL specification, i.e., the
%            `GRID` section of an ECLIPSE input deck.
%
%            With the exception of `SPECGRID` whose first three arguments
%            (grid cell dimensions NX, NY, NZ) are stored in `grdecl`
%            structure field `cartDims`, all currently recognized keywords
%            are stored in equally named `grdecl` structure fields.
%            Specifically, GRDECL keyword `COORD` is stored `grdecl`
%            structure field `COORD` and so forth.
%
%            The pillar description `COORD` is stored in a 6*nPillar
%            array (number of pillars, nPillar == (NX+1)*(NY+1)) of
%            bottom/top coordinate pairs.  Specifically, ::
%
%              grdecl.COORD((i-1)*6+(1:3)) -- (x,y,z)-coordinates of
%                                             pillar 'i' top point.
%              grdecl.COORD((i-1)*6+(4:6)) -- (x,y,z)-coordinates of
%                                             pillar 'i' bottom point.
%
%            The currently recognized keywords are::
%
%            'ACTNUM', 'COORD', 'DXV', 'DYV', 'DZV', 'DEPTHZ',
%            'DIMENS', 'INCLUDE', 'MULTX', 'MULTX-', 'MULTY', 'MULTY-',
%            'MULTZ', 'MULTZ-', 'NOGRAV', 'NTG', 'PERMX', 'PERMXY',
%            'PERMXZ', 'PERMY', 'PERMYX', 'PERMYZ', 'PERMZ', 'PERMZX',
%            'PERMZY', 'PORO', 'ROCKTYPE', 'SATNUM', 'ZCORN'
%
%   unrec -  Cell array with unrecognized keywords.  For convenience and
%            greater transparency, this data is also provided as the field
%            `UnhandledKeywords` of `grdecl` itself.
%
% SEE ALSO:
%   `processGRDECL`.

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

opt = struct('verbose'        , mrstVerbose, ...
             'cartDims'       , []         , ...
             'grdecl'         , []         , ...
             'keywords'       , {{}}       , ...
             'missing_include', @error);
opt = merge_options(opt, varargin{:});

if isempty(opt.grdecl)
   grdecl = struct;
else
   grdecl = opt.grdecl;
end

[fid, msg] = fopen(fn, 'rt');
if fid < 0, error([fn, ': ', msg]), end

cartDims = opt.cartDims;
numCell  = [-1, prod(cartDims)];
numCell  = numCell(1 + double(~isempty(cartDims)));

unrec = [];

while ~feof(fid)
   lin = fgetl(fid);
   if lin == -1
      msg = ferror(fid, 'clear');
      fclose(fid);
      error('readGRDECL:Input:Empty', ...
            'GRDECL file ''%s'' is unreadable.\nSystem reports: %s\n', ...
            fn, msg);
   end

   % Loop until next keyword
   kw = regexp(lin, '^[A-Z][A-Z0-9]{0,7}(|/)', 'match', 'once');
   while isempty(kw) && ~feof(fid)
      lin = fgetl(fid);
      if lin ~= -1
         kw = regexp(lin, '^[A-Z][A-Z0-9]{0,7}(|/)', 'match', 'once');
      end
   end

   if ~feof(fid)
      dispif(opt.verbose, 'Reading keyword ''%s''\n', kw);
      t0 = tic();

      switch kw
         case {'SPECGRID', 'DIMENS'}
            t = fscanf(fid, '%f', 3) .';
            if isempty(cartDims)
               cartDims = t;
            elseif ~all(cartDims == t)
               error(msgid('Grid:Initialized'), ...
                     'Attempting to change already defined grid size.');
            end
            numCell  = prod(cartDims);
            trash    = fgetl(fid); %#ok
            grdecl.cartDims = reshape(cartDims, 1, []);

         case 'INCLUDE'
            % checkDim(cartDims, numCell, kw, fid);
            inc_fn_tmp = fscanf(fid, '%s', 1);
            inc_fn = regexp(inc_fn_tmp, '[''"]?([-./\w]+)[''"]?', ...
                            'tokens', 'once');
            inc_fn = inc_fn{1};

            terminated = inc_fn_tmp(end) == '/';
            if inc_fn(end) == '/', inc_fn = inc_fn(1:end-1); end

            % Gobble up keyword-closing '/' character if not already read.
            if ~terminated
               slash = fscanf(fid, '%s', 1);
               if ~strcmp(slash, '/')
                  error(msgid('Include:WrongfulTermination'), ...
                        'INCLUDE keyword not correctly terminated.');
               end
            end

            inc_fn(inc_fn == '/') = filesep;
            if inc_fn(1) ~= filesep
               % Translate relative pathname to absolute pathname.
               inc_fn = fullfile(fileparts(fn), inc_fn);
            end

            if ~ exist(inc_fn, 'file')
               msg = sprintf('Include file ''%s'' missing', inc_fn);
               opt.missing_include('INCLUDE:Missing', '%s', msg);
               continue
            end

            dispif(opt.verbose, ' -> ''%s''\n', inc_fn);
            [grdecl, u] = readGRDECL(inc_fn, 'cartDims', cartDims, ...
                                     'verbose', opt.verbose,       ...
                                     'grdecl', grdecl,             ...
                                     'missing_include', ...
                                     opt.missing_include);
            dispif(opt.verbose, ' <- ''%s''\n', inc_fn);

            unrec = [ unrec, u ];                               %#ok<AGROW>

            if isempty(cartDims) && isfield(grdecl, 'cartDims') && ...
                  numel(grdecl.cartDims) == 3
               cartDims        = reshape(grdecl.cartDims, 1, []);
               grdecl.cartDims = cartDims;
            end

         case 'COORD'
            checkDim(cartDims, numCell, kw, fid);
            grdecl.COORD = readVector(fid, kw, ...
                                      6 * prod(cartDims(1:2) + 1));

         case 'ZCORN'
            checkDim(cartDims, numCell, kw, fid);
            grdecl.ZCORN = readVector(fid, kw, 8 * numCell);

         case {'PORO',                         ...
               'PERMX' , 'PERMXY', 'PERMXZ',   ...
               'PERMYX', 'PERMY' , 'PERMYZ',   ...
               'PERMZX', 'PERMZY', 'PERMZ' ,   ...
               'PERMH',                        ...
               'ACTNUM', 'SATNUM', 'ROCKTYPE', ...
               'MULTX' , 'MULTX-',             ...
               'MULTY' , 'MULTY-',             ...
               'MULTZ' , 'MULTZ-',             ...
               'NTG'   , 'VSH'   ,             ...
               }
            checkDim(cartDims, numCell, kw, fid);
            grdecl.(regexprep(kw, '\W', '_')) = ...
               readVector(fid, kw, numCell);

         case {'DXV', 'DYV', 'DZV'}
            ix          = strcmp(kw, {'DXV', 'DYV', 'DZV'});
            grdecl.(kw) = readVector(fid, kw, cartDims(ix));

         case 'DEPTHZ'
            checkDim(cartDims, numCell, kw, fid);
            grdecl.(kw) = readVector(fid, kw, prod(cartDims(1:2) + 1));

         case 'MAPAXES'
            grdecl.(kw) = readVector(fid, kw, 6);

         case 'MAPUNITS'
            data = readDefaultedRecord(fid, { 'METRES' });
            grdecl.(kw) = data{1};  clear data

         case {'ADD', 'COPY', 'EQUALS', 'MAXVALUE', ...
               'MINVALUE', 'MULTIPLY'}
            grdecl = ...
                applyOperatorSimple(grdecl, fid, cartDims, cartDims, kw);

         case 'FAULTS'
            tmpl(1:8) = { 'Default' };
            data = readDefaultedKW(fid, tmpl);  clear tmpl
            data(:, 2:end-1) = to_double(data(:, 2:end-1));

            if ~isfield(grdecl, kw), grdecl.(kw) = cell([0, 8]); end
            grdecl.(kw) = [grdecl.(kw); data];

         case 'MULTFLT'
            tmpl = { 'FaultName', '1.0', '1.0' };
            data = readDefaultedKW(fid, tmpl);  clear tmpl
            data(:, 2:end) = to_double(data(:, 2:end));

            if ~isfield(grdecl, kw), grdecl.(kw) = cell([0, 3]); end
            grdecl.(kw) = [grdecl.(kw); data];

         case 'GDORIENT'
            tmpl = {'INC', 'INC', 'INC', 'DOWN', 'RIGHT'};

            grdecl.(kw) = readDefaultedRecord(fid, tmpl);        clear tmpl

         case 'GRIDUNIT'
            grdecl.(kw) = readDefaultedRecord(fid, { 'METRES', '' });

         otherwise
            if ~isempty(opt.keywords) && any(strcmp(kw, opt.keywords))
               checkDim(cartDims, numCell, kw, fid);
               grdecl.(genvarname(kw)) = readVector(fid, kw, numCell);
            else
               unrec = [unrec, {kw}]; %#ok
            end
      end

      t0 = toc(t0);
      dispif(opt.verbose, '  -> Done [%.2f s]\n', t0);
   end
end

unrec = unique(unrec);
grdecl.UnhandledKeywords = unrec;

fclose(fid);
if isfield(grdecl, 'ACTNUM')
   grdecl.ACTNUM = int32(grdecl.ACTNUM);
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function checkDim(cartDims, numCell, kw, fid)
if isempty(cartDims) || numCell < 1 || any(cartDims < 1)
   % Don't leave open fd's in MATLAB's workspace when erroring out.
   fclose(fid);
   error('readGRDECL:Input:NoDim', ...
         'GRDECL keyword ''%s'' found before dimension specification', kw);
end

%--------------------------------------------------------------------------

function v = to_double(v)
convert = @(s) sscanf(s, '%f');
if ischar(v)
   v = convert(v);
else
   v = cellfun(convert, v, 'UniformOutput', false);
end
