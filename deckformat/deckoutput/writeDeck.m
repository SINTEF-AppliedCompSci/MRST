function writeDeck(deck, directory, varargin)
% Write deck-structure to files.
%
% SYNOPSIS:
%   writeDeck(deck, dirname)
%   writeDeck(deck, dirname, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function writes the contents of selected keywords to individual
%   files in a directory specified by the user.
%
% REQUIRED PARAMETERS:
%   deck      -  structure that (may) contain fields RUNSPEC, GRID, EDIT, PROPS, 
%                REGIONS, SOLUTION and SCHEDULE  
%
%   directory - name of directory in which to write data.  If the directory
%              does not exist, it will be created.
%
% OPTIONAL PARAMETERS:
%   'unit'   - Output choice of unit. Default is 'metric'.
%
%   Remaining options are formats passed to fprintf (see file for details) 
% 
% SEE ALSO:
%   model2Deck, writeSchedule

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

if isempty(varargin)
    dispif(mrstVerbose, 'Output unit system defaulted to METRIC\n');
end

opt = struct('unit', 'metric');
[opt, format_opts] = merge_options(opt, varargin{:});
deck = convertDeckUnits(deck, 'outputUnit', opt.unit);

formats = struct('case',         'upper', ...
                 'int',            '%8d', ...
                 'string',        '%10s', ...
                 'double',      '%12.6f', ...
                 'sci',         '%15.7e', ...
                 'doubleRange', [.1 100]);

formats = merge_options(formats, format_opts{:});
if ~isa(formats.case, 'function_handle')
    formats.case = str2func(formats.case);
end

% organize directory
if ~(exist('isfolder', 'file') && isfolder(directory)) || ~isdir(directory)
    mkdir(directory);
end

[~, fn] = fileparts(directory);                            
fn      = fullfile(directory, [upper(fn), '.DATA']);

[fid, msg] = fopen(fn, 'w+');
if fid < 0
    error(msgid('Open:Failure'), 'Failed to open ''%s'': %s\n', fn, msg);
end

% 
dashedline = [repmat('-', [1, 57]), '\n'] ;

% DATA-file header
printHeader(fid);
fprintf('\n')
%--- RUNSPEC -------------------------------------------------------------- 
fprintf(fid, dashedline);
fprintf(fid,'%s\n','RUNSPEC');
fprintf(fid, dashedline);
if(isfield(deck.RUNSPEC,'TITLE'))
    fprintf(fid, '%s\n', 'TITLE');
    fprintf(fid, '%s\n', deck.RUNSPEC.TITLE);
end
fprintf(fid, '\n');
dump_runspec(fid, directory, deck, formats);

%--- GRID ----------------------------------------------------------------- 
fprintf(fid, dashedline);
fprintf(fid, '%s\n', 'GRID');
fprintf(fid, dashedline);
fprintf(fid, '%s\n', 'INIT');
dump_grid(fid,directory, deck, formats);

%--- EDIT -----------------------------------------------------------------
if(isfield(deck,'EDIT'))
    fprintf(fid, dashedline);
    fprintf(fid, '%s\n', 'EDIT');
    fprintf(fid, dashedline);
    dump_edit(fid,directory, deck, formats);
end

%--- PROPS ----------------------------------------------------------------
% fixing rock for now for opm runs E300 (?)
if isfield(deck.PROPS, 'ROCK')
    deck.PROPS.ROCK = deck.PROPS.ROCK(:,1:2);
    fprintf(fid, dashedline);
end
fprintf(fid, '%s\n', 'PROPS');
fprintf(fid, dashedline);

dump_props(fid, directory, deck, formats);

%--- REGIONS --------------------------------------------------------------
if(isfield(deck,'REGIONS'))
    fprintf(fid, dashedline);
    fprintf(fid, '%s\n', 'REGIONS');
    fprintf(fid, dashedline);
    dump_regions(fid,directory, deck, formats)
end

%--- SOLUTION -------------------------------------------------------------
fprintf(fid, dashedline);
fprintf(fid, '%s\n', 'SOLUTION');
fprintf(fid, dashedline);
dump_solution(fid,directory, deck, formats);

%--- SUMMARY --------------------------------------------------------------
fprintf(fid, dashedline);
fprintf(fid, '%s\n', 'SUMMARY');
fprintf(fid, dashedline);
if(isfield(deck,'UnhandledKeywords'))
    if(isfield(deck.UnhandledKeywords,'SUMMARY'))
        myfields=deck.UnhandledKeywords.SUMMARY;
        fprintf(fid,'\n');
        for i=1:numel(myfields)
            fprintf(fid,'%s\n/\n\n',myfields{i});
        end
    end
end

%--- SCHEDULE -------------------------------------------------------------
fprintf(fid, dashedline);
fprintf(fid, '%s\n', 'SCHEDULE');
fprintf(fid, dashedline);
%dump_schedule(fid,directory, deck);
% external function
writeSchedule(fid, directory, deck.SCHEDULE, 'writeInclude', true, ...
                                             'includeName',  'SCHEDULE.TXT', ...
                                             'formats',       formats)
fclose(fid);
end

%--------------------------------------------------------------------------
%--- ----------------------------------------------------------------------

function dump_runspec(fid, dirname, deck, f)
flds = {'DIMENS','EQLDIMS'};%,'WELLDIMS'};%,'NUPCOL'}
for i=1:numel(flds)
    fld = flds{i};
    if(isfield(deck.RUNSPEC, fld))
        v   = deck.RUNSPEC.(fld);
        fmt = getFmtStr(f.int, numel(v)); 
        dump_vector(fid, dirname, lower(fld), fmt, v, false);
    end
end
if(isfield(deck.RUNSPEC,'TABDIMS'))
    v   = deck.RUNSPEC.TABDIMS(1:6);
    fmt = getFmtStr(f.int, numel(v)); 
    dump_vector(fid,dirname, 'tabdims', fmt, v, false);
end
if(isfield(deck.RUNSPEC,'WELLDIMS'))
    v   = deck.RUNSPEC.WELLDIMS(1:4);
    fmt = getFmtStr(f.int, numel(v)); 
    dump_vector(fid,dirname, 'welldims', fmt, v, false);
end
flds={'OIL','WATER','GAS','DISGAS','VAPOIL','METRIC','NOGRAV','FIELD','UNIFOUT','FMTOUT'};
for i=1:numel(flds)
    fld = flds{i};
    if(isfield(deck.RUNSPEC, fld))
        if(deck.RUNSPEC.(fld)==1)
            fprintf(fid, '%s\n',fld);
        end
    end
end
if isfield(deck.RUNSPEC, 'GRIDOPTS')
    fmt = getFmtStr(f.string, f.int, f.int);
    dump_vector(fid,dirname, 'gridopts', fmt, deck.RUNSPEC.GRIDOPTS(1:3), false);
end
if isfield(deck.RUNSPEC, 'ENDSCALE')
    fmt = getFmtStr(f.string, f.string, f.int, f.int);
    dump_vector(fid,dirname, 'endscale', fmt , deck.RUNSPEC.ENDSCALE(1:4), false);
end
if isfield(deck.RUNSPEC, 'START')
    fprintf(fid,'START \n');
    dat = datestr(deck.RUNSPEC.START);
    dat = regexp(dat,'-','split');
    fprintf(fid, ' %s ''%s'' %s\n/\n\n', dat{1}, upper(dat{2}), dat{3});
end
end

%-------------------------------------
function dump_edit(fid,dirname, deck, f)
flds = {'DEPTH','PORV', 'TRANX', 'TRANY', 'TRANZ'};
for i=1:numel(flds)
    if(isfield(deck.EDIT,flds{i}))
        dump_vector(fid, dirname, lower(flds{i}), [f.sci, '\n'], deck.EDIT.(flds{i}));
    end
end

end

%-------------------------------------
function dump_regions(fid,dirname, deck, f)
flds = fieldnames(deck.REGIONS);
for i=1:numel(flds)
    if(isfield(deck.REGIONS,flds{i}))
        dump_vector(fid, dirname, lower(flds{i}), [f.int, '\n'], deck.REGIONS.(flds{i}));
    end
end
end


%--------------------------------------------------------------------------

function dump_grid(fid, dirname, deck, f)
grid_support = true;

if isfield(deck.GRID, 'COORD')
    assert (isfield(deck.GRID, 'ZCORN'));
    coord = deck.GRID.COORD;
    zcorn = deck.GRID.ZCORN;
elseif all(isfield(deck.GRID, {'DXV', 'DYV', 'DZV'}))
    [coord, zcorn] = block_centred_to_cpg(deck);
else
    grid_support = false;
    if all(isfield(deck.GRID, {'DX', 'DY', 'DZ'}))
        fprintf(fid, 'SPECGRID\n');
        fprintf(fid, '%d ', deck.GRID.cartDims);
        fprintf(fid, '1 F\n/\n\n');
        assert(isfield(deck.GRID, 'TOPS'));
        dump_vector(fid, dirname, 'dx', getFmtStr(f.sci), deck.GRID.DX);
        dump_vector(fid, dirname, 'dy', getFmtStr(f.sci), deck.GRID.DY);
        dump_vector(fid, dirname, 'dz', getFmtStr(f.sci), deck.GRID.DZ);
        dump_vector(fid, dirname, 'tops', getFmtStr(f.sci), deck.GRID.TOPS);
    end
end

if grid_support
    fprintf(fid, 'SPECGRID\n');
    fprintf(fid, '%d ', deck.GRID.cartDims);
    fprintf(fid, '1 F\n/\n\n');
    
    dump_vector(fid, dirname, 'coord', getFmtStr(f.sci) , coord);
    dump_vector(fid, dirname, 'zcorn', getFmtStr(f.sci), zcorn);
end
% write other data from deck
flds = {'PERMX', 'PERMY', 'PERMZ', 'PORO', 'ACTNUM'};
for k = 1:numel(flds)
    if isfield(deck.GRID, flds{k})
        if strcmp(flds{k}, 'ACTNUM')
            fmt = getFmtStr(f.int);
        else
            fmt = getFmtStr(f.sci);
        end
        dump_vector(fid, dirname, flds{k}, fmt, deck.GRID.(flds{k}));
    end
end

if isfield(deck.GRID, 'NNC')
    nnc = deck.GRID.NNC;
    fmt = getFmtStr(f.int, f.int, f.int, f.int, f.int, f.int, f.sci, '/');
    dump_vector(fid, dirname, 'nnc' , fmt, nnc');
end
end

%--------------------------------------------------------------------------

function [coord, zcorn] = block_centred_to_cpg(deck)
dims = deck.RUNSPEC.DIMENS;
n    = prod(dims(1:2) + 1);

x = cumsum([0; deck.GRID.DXV]);
y = cumsum([0; deck.GRID.DYV]);
z = cumsum([0; deck.GRID.DZV]);

[X, Y, Z] = ndgrid(x, y, z);

if isfield(deck.GRID, 'DEPTHZ')
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

function dump_props(fid, dirname, deck, f)
for fld = reshape(fieldnames(deck.PROPS), 1, [])
    values = deck.PROPS.(fld{1});
    
    if any(strcmp(fld{1}, {'PVTO','PVTG'}))
        % Custom output routine
        
            dump_pvt(fid,lower(fld{1}), dirname, values, f);
        
        %continue
%     elseif any(strcmp(fld{1}, {'SCALECRS'}))
%         fprintf(fid, '%s\n', fld{1});
%         fprintf(fid, '%s\n', values{1});
%         fprintf(fid, '\n'); 
    else
        values = convertValuesToCell(values);
    
        for regn = 1:numel(values)
            values{regn}(isnan(values{regn})) = 0;
        end
        if iscellstr(values)
            dump_vector(fid,dirname, lower(fld{1}), getFmtStr(f.string), values);
        else
            fmt = getFmtStr(f.sci, size(values{1}, 2));
            dump_multiple(fid, dirname, lower(fld{1}), fmt, values);
        end
    end
end
end

function values = convertValuesToCell(values)
if ~iscell(values)
    if(size(values,2)>1)
        tmp = cell(size(values,1),1);
        for i=1:numel(tmp)
            tmp{i} = values(i,:);
        end
        values = tmp;
    else
        values = {values};
    end
end
end

%--------------------------------------------------------------------------
% function dump_pvto(fid, dirname, pvto, f)
% fncase = @upper;
% fprintf(fid, ['INCLUDE\n', fncase('pvto.txt'),  '/\n\n']);
% 
% [fid_pvto, msg] = fopen(fullfile(dirname, fncase('pvto.txt')), 'wt');
% if fid_pvto < 0
%     error('Failed to open ''pvto'' output file: %s', msg);
% end
% 
% assert (numel(pvto.key) + 1 == numel(pvto.pos));
% 
% fprintf(fid_pvto, 'PVTO\n');
% 
% for r = 1 : numel(pvto.pos) - 1
%     fprintf(fid_pvto, [f.sci, '\n'], pvto.key(r));
%     
%     i = pvto.pos(r) : pvto.pos(r + 1) - 1;
%     fprintf(fid_pvto, getFmtStr(f.sci, 3), pvto.data(i,:).');
%     fprintf(fid_pvto, '/\n');
% end
% 
% fprintf(fid_pvto, '/\n');
% fclose(fid_pvto);
% end

%--------------------------------------------------------------------------

function dump_pvt(fid,name, dirname, values, f)
fncase = @upper;
filename = fncase([name,'.txt']);
fprintf(fid, ['INCLUDE\n''', filename, ''' /\n\n']);

[fid_pvt, msg] = fopen(fullfile(dirname, filename), 'wt');
if fid_pvt < 0
    error('Failed to open ''pvt'' output file: %s', msg);
end



fprintf(fid_pvt, [upper(name),'\n']);
for i=1:numel(values)
    pvt=values{i};
    assert (numel(pvt.key) + 1 == numel(pvt.pos));
    for r = 1 : numel(pvt.pos) - 1
        fprintf(fid_pvt, [f.sci, '\n'], pvt.key(r));
    
    
        i = pvt.pos(r) : pvt.pos(r + 1) - 1;
        fprintf(fid_pvt, getFmtStr(f.sci, 3), pvt.data(i,:).');
        fprintf(fid_pvt, '/\n');
    end
    fprintf(fid_pvt, '/\n');
end
fclose(fid_pvt);
end

%--------------------------------------------------------------------------

function dump_solution(fid,dirname, deck, f)
myfields=fieldnames(deck.SOLUTION);
for i=1:numel(myfields)
    myfield=myfields{i};
    
    v = deck.SOLUTION.(myfield);
    %if iscell(v), v = v{1}; end
    if strcmp(myfield,'EQUIL')
        v=v(:,1:9);
        fmt = getFmtStr(f.sci, f.sci, f.sci, f.sci, f.sci, f.sci, f.int, f.int, f.int);
        dump_multiple(fid,dirname, lower(myfield), fmt, v);
    elseif strcmp(myfield,'THPRES')
        v=v(:,1:3);
        fmt = getFmtStr(f.int, f.int, f.sci,'/');
        %dump_multiple(fid,dirname, lower(myfield), fmt, v, );
        dump_vector(fid, dirname, lower(myfield), fmt, v');
    elseif strcmp(myfield, 'AQUFETP')
        % might have nans
         v = v(:, 1:9);
         fmts = {f.int, f.sci, f.sci, f.sci, f.sci, f.sci, f.int, f.sci, f.sci};
%         fmts = repmat({f.sci}, [1,9]);
         [v, nix] = replace_nan(v);
         fmts(nix) = {f.string};
         dump_vector(fid, dirname, lower(myfield), getFmtStr(fmts{:}, ' /'), v');
    elseif strcmp(myfield, 'AQUANCON')
        % replace negative multiplyer
        v = v(:, 1:11);
        fmts = {f.int, f.int, f.int, f.int, f.int, f.int, f.int, f.string, f.sci, f.double, f.string};
        [v, nix]  = replace_negative(v);
        fmts(nix) = {f.string};
        dump_vector(fid, dirname, lower(myfield), getFmtStr(fmts{:}, ' /'), v');    
    else
        values = v;
        values = convertValuesToCell(values);
        for regn = 1:numel(values)
            values{regn}(isnan(values{regn})) = 0;
        end
        if iscellstr(values)
            dump_vector(fid,dirname, lower(myfield), getFmtStr(f.string), values);
        else
            fmt = getFmtStr(f.sci, size(values{1}, 2));
            dump_multiple(fid,dirname, lower(myfield), fmt, values);
        end
    end
end
end

%--------------------------------------------------------------------------

function dump_vector(fid, dirname, field, fmt, v, newFile)
if nargin < 6
    newFile = true;
end
fncase = @upper;
org_fid = fid;

if newFile
    fn         = fullfile(dirname, fncase([field, '.txt']));
    [fid, msg] = fopen(fn, 'wt');
    if fid < 0
        error('Unable to open ''%s'' for writing: %s', fn, msg);
    end
end

fprintf(fid, '%s\n', upper(field));
if ~iscell(v)
    fprintf(fid, fmt, v);
else
    fprintf(fid, fmt, v{:});
end
fprintf(fid, '/\n\n');

if newFile
    fclose(fid);
    fprintf(org_fid, '%s\n', 'INCLUDE');
    fprintf(org_fid, '''%s''\n', fncase([field, '.txt']));
    fprintf(org_fid, '/\n\n');
end
end

%--------------------------------------------------------------------------

function dump_multiple(fid, dirname, field, fmt, v)
fncase = @upper;
org_fid = fid;

fn         = fullfile(dirname, fncase([field, '.txt']));
[fid, msg] = fopen(fn, 'wt');
if fid < 0
    error('Unable to open ''%s'' for writing: %s', fn, msg);
end
fprintf(fid, '%s\n', upper(field));
if iscellstr(v)
   fprintf(fid, fmt, v{:});
   fprintf(fid, '/\n');
elseif iscell(v)
    for k = 1:numel(v)
        fprintf(fid, fmt, v{k}');
        fprintf(fid, '/\n');
    end
else
    for k = 1:size(v,1)
        fprintf(fid, fmt, v(k,:));
        fprintf(fid, '/\n');
    end
end
fclose(fid);

fprintf(org_fid, '%s\n', 'INCLUDE');
fprintf(org_fid, '''%s''\n', fncase([field, '.txt']));
fprintf(org_fid, '/\n\n');
end

%--------------------------------------------------------------------------

function [t, ix] = replace_negative(t)
if ~iscell(t)
    ix = t(1,:) < 0;
    if any(ix)
        t = num2cell(t);
    end
else
    ix = cellfun(@(x)isnumeric(x)&&x<0, t(1,:));
end
    if any(ix)
        t(:,ix) = repmat({'1*'}, [size(t,1), nnz(ix)]);
    end
end

%--------------------------------------------------------------------------

function [t, ix] = replace_nan(t)
if ~iscell(t)
    ix = isnan(t(1,:));
    if any(ix)
        t = num2cell(t);
    end
else
    ix = cellfun(@(x)isnumeric(x)&&isnan(x), t(1,:));
end
if any(ix)
    t(:,ix) = repmat({'1*'}, [size(t,1), nnz(ix)]);
end
end

%--------------------------------------------------------------------------

function str = getFmtStr(varargin)
n = 1;
str = varargin;
if isnumeric(str{end})
    n   = str{end};
    str = str(1:end-1);
end
str = cellfun(@(s)[s, ' '], str, 'UniformOutput', false);
str = horzcat(str{:});

if n > 1
    str = repmat(str, [1, n]);
end
str = [str, '\n'];
end

%--------------------------------------------------------------------------

function printHeader(fid)
str = ...
{'---------------------------------------------------------', ...
 '---                                                   ---', ...                     
 '---           M   M   RRRR     SSSS   TTTTT           ---', ... 
 '---           MM MM   R   R   S         T             ---', ...  
 '---           M M M   RRRR     SSS      T             ---', ...     
 '---           M   M   R R         S     T             ---', ...     
 '---           M   M   R  RR   SSSS      T             ---', ...
 '---                                                   ---', ...
 '--------------------------------- www.sintef.no/mrst  ---'};
fprintf(fid, '%s\n', str{:});
fprintf(fid, '\n%s\n\n', '-- Generated deck from MRST function writeDeck');
end
