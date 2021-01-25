function [smry, smspec] = readEclipseSummaryUnFmt(prefix, keyWords)
%Undocumented utility function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

smspec = readEclipseOutputFileUnFmt([prefix, '.SMSPEC']);

if isfield(smspec, 'WGNAMES')
   nameList = smspec.WGNAMES.values;
elseif isfield(smspec, 'NAMES')
   nameList = smspec.NAMES.values;
else
   error('SmryNames:Unknown', ...
        ['No known summary vector name list in SMSPEC fields\n\t', ...
         '{ %s}'], sprintf('%s ', sort(fieldnames(smspec))));
end
kwrdList = smspec.KEYWORDS.values;

% replace empty entries by 'empty'
[nameList{strcmp(nameList, '')}] = deal('empty');
[kwrdList{strcmp(kwrdList, '')}] = deal('empty');

numFull = numel(nameList);

if nargin > 1
    keyWords = keyWords(:);
    rowInx = ismember(kwrdList, keyWords);
    nameList = nameList(rowInx);
    kwrdList = kwrdList(rowInx);
else
    rowInx = (1:numel(nameList))';
end

[names, nInx, nInx] = unique(nameList);                                %#ok
[kwrds, kInx, kInx] = unique(kwrdList);                                %#ok
nlist  = numel(nameList);

smry_file = [prefix, '.UNSMRY'];
if exist(smry_file, 'file') % unified summary
    smry_files = {smry_file};
    %estimate number of ministeps
    d = dir(smry_file);
    estNum = floor(d.bytes/(4*numFull));
else % multiple summary-files
    [dname, fp] = fileparts(prefix);
    if isempty(dname)
        dname = '.';
    end
    smry_files = matchResultFiles(dname, [fp, '\.S\d{4}']);
    %estimate number of ministeps
    nbts = 0;
    for k = 1:numel(smry_files)
        d = dir(smry_files{k});
        nbts = nbts + d.bytes;
    end
    estNum = floor(nbts/(4*numFull));
end

data = zeros(nlist, estNum);
ministeps = nan(1, estNum);
curStep = 0;

dispif(mrstVerbose, ['Reading info from roughly ', num2str(estNum), ' ministeps:      '])
for f = reshape(smry_files, 1, [])
    [fid, msg] = fopen(f{1}, 'r', 'ieee-be');
    if fid < 0, error([f{1}, ': ', msg]); end
    % jump to start
    fseek(fid, 4, 'cof');
    while ~feof(fid)
        [name, field] = readFieldUnFmt(fid);
        if strcmp(name, 'MINISTEP')
            curStep = curStep +1;
            ministeps(curStep) = field.values;
            if mod(curStep, 100) == 0
                dispif(mrstVerbose, '\b\b\b\b%3d%%', round(100*curStep/estNum));
            end
        end
        if strcmp(name, 'PARAMS')
            data(:, curStep) = field.values(rowInx);
            %ministeps(ministep+1) = true;
        end
    end
    fclose(fid);
end

smry.WGNAMES  = names;
smry.KEYWORDS = kwrds;
smry.UNITS    = smspec.UNITS.values;
smry.STARTDAT = smspec.STARTDAT.values([3 2 1])';
smry.nInx     = nInx;
smry.kInx     = kInx;
smry.data     = data(:, 1:curStep);
smry.ministep = ministeps(1:curStep);

smry = addSmryFuncs(smry);

dispif(mrstVerbose, '\b\b\b\b%3d%%', 100);
dispif(mrstVerbose, ['\nActual number of ministeps: ', num2str(curStep), '\n']);
end


%--------------------------------------------------------------------------

function smry = addSmryFuncs(smry)
smry.get    = @(nm,kw,ms)getData(smry, nm, kw, ms);
smry.getInx = @(nm,kw)getRowInx(smry,nm,kw);
smry.getNms = @(kw)getNames(smry,kw);
smry.getKws = @(nm)getKeywords(smry, nm);
smry.getUnit= @(nm,kw)getUnit(smry, nm, kw);
end

function nms = getNames(smry,kw)
rInx   = getRowInx(smry, [], kw);
nmPos  = unique(smry.nInx(rInx));
nms    = smry.WGNAMES(nmPos);
end

function kws = getKeywords(smry,nm)
rInx   = getRowInx(smry, nm, []);
kwPos  = unique(smry.kInx(rInx));
kws    = smry.KEYWORDS(kwPos);
end

function s = getData(smry, nm, kw, ms)
[rInx, order] = getRowInx(smry, nm, kw);
if isempty(ms) || (ischar(ms) && strcmp(ms, ':'))
   ms = 1 : size(smry.data, 2);
end
s = smry.data(rInx, ms);
s = s(order,:);
end

function [rInx, order] = getRowInx(smry, nm, kw)
nm_multi = iscellstr(nm) && (numel(nm) > 1);
kw_multi = iscellstr(kw) && (numel(kw) > 1);
assert (~ (nm_multi && kw_multi), ...
       ['At most one of Keyword or Entity may be a cell array ', ...
        'of character vectors']);
nlist = numel(smry.kInx);
if ~isempty(nm)
    if ischar(nm), nm = { nm }; end
    [pick, nm_order] = string_lookup(smry.WGNAMES, nm);
    rInxN = pick(smry.nInx);
else
    rInxN = true(nlist, 1);
end
if ~isempty(kw)
    if ischar(kw), kw = { kw }; end
    [pick, kw_order] = string_lookup(smry.KEYWORDS, kw);
    rInxK = pick(smry.kInx);
else
    rInxK = true(nlist, 1);
end
rInx = and(rInxN, rInxK);
if nm_multi
    order = nm_order;
elseif kw_multi
    order = kw_order;
else
    order = (1 : sum(rInx)) .';
end
end

function u = getUnit(smry, nm, kw)
rInx = getRowInx(smry, nm, kw);
u = smry.UNITS{rInx};
end

%--------------------------------------------------------------------------

function [pick, order] = string_lookup(list, search)
ncol = numel(search);
[i, j] = blockDiagIndex(numel(list), ncol);
match = strcmp(reshape(list(i)  , [], ncol), ...
               reshape(search(j), [], ncol));
pick = any(match, 2);
[row, ~] = find(match);
[~, reverse_order] = sort(row);
order = zeros(size(reverse_order));
order(reverse_order) = 1 : numel(reverse_order);
end
