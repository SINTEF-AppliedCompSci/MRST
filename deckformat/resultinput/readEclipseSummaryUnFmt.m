function [smry, smspec] = readEclipseSummaryUnFmt(prefix, keyWords)
%Input Summary Vectors From Unformatted ECLIPSE Result Files
%
% SYNOPSIS:
%    smry          = readEclipseSummaryUnFmt(prefix)
%    smry          = readEclipseSummaryUnFmt(prefix, keywords)
%   [smry, smspec] = readEclipseSummaryUnFmt(...)
%
% PARAMETERS:
%   prefix   - Casename prefix (character vector) such that
%                 [ prefix, '.SMSPEC' ]
%              names the corresponding summary specification file.  This
%              function knows about both unified (*.UNSMRY) and separate
%              (*.S000n) summary result files.
%
%   keywords - List, represented as a cell-array of character vectors, of
%              summary vector names to load from the result set.  OPTIONAL.
%              Load all summary results if not present.
%
% RETURNS:
%   smry   - Summary result aggregation structure.  Contains at least the
%            following fields:
%              - WGNAMES  - List (cell-array of character vectors) of named
%                           entities (e.g., wells or groups) in the result
%                           set.  Typically contains at least the
%                           "sentinel"/no-name string ':+:+:+:+' associated
%                           to vectors that don't have a name (e.g., TIME).
%
%              - KEYWORDS - List (cell-array of character vectors) of named
%                           summary vector keywords (e.g., 'TIME' or 'WOPT').
%
%              - UNITS    - List (cell-array of character vectors) of named
%                           summary vector units.
%
%              - nInx     - Name index.  Index into WGNAMES.
%
%              - kInx     - Keyword index.  Index into KEYWORDS.
%
%              - data     - Summary results.  Data(i,j) is the numerical
%                           value of summary vector
%                              [KEYWORDS{kInx(i)}, ':', WGNAMES{nInx(i)}]
%                           at time (ministep) 'j'.
%
%              - get      - Result data accessor.  Function handle.  Must be
%                           called as
%                              x = smry.get(names, keywords, ministeps)
%                           in which 'names' is a named entity (e.g., a
%                           well), 'keywords' is a set of summary vectors
%                           for that named entity, and 'ministeps' is the
%                           set of ministeps (columns) for which to extract
%                           the result values.  If 'ministeps' is an empty
%                           array or the character vector ':' (colon), then
%                           returns all time values for the selected
%                           summary vector.
%
%                           Either of 'names' or 'keywords', but not both,
%                           may be a cell-array of character vectors in
%                           which case GET will return a single vector
%                           (e.g., oil production rate) for multiple named
%                           entities or multiple vectors for a single named
%                           entity.  The values are returned in the same
%                           order as the input cell-array.
%
%              - getUnit  - Unit string accessor.  Function handle.  Must
%                           be called as
%                              ustr = smry.getUnit(name, keyword)
%                           which will return the appropriate unit string
%                           for the named summary vector.
%
%   smspec - Raw summary specification structure as defined by function
%            `readEclipseOutputFileUnFmt`.
%
% EXAMPLES:
%   smry = readEclipseSummaryUnFmt('CASE');
%
%   % 1) Plot well-level oil production rate for all wells named PROD*.
%   time  = smry.get(repmat(':+', [1, 4]), 'TIME', []); % name = ':+:+:+:+'
%   wells = regexp(smry.WGNAMES, '^PROD');
%   wells = smry.WGNAMES(~ cellfun('isempty', wells));
%
%   wopr = smry.get(wells, 'WOPR', []);
%   plot(time, wopr, 'x-')
%   title(['Oil Production Rate [', smry.getUnit(wells{1}, 'WOPR'), ']'])
%   legend(wells, 'Location', 'Best')
%
%   % 2) Plot oil, water, and liquid production rate for well PROD1.
%   time = smry.get(repmat(':+', [1, 4]), 'TIME', []); % name = ':+:+:+:+'
%   kws  = {'WOPR', 'WWPR', 'WLPR'};
%
%   wxpr = smry.get('PROD1', kws, []);
%   plot(time, wxpr, 'x-')
%   title('Production Rates')
%   legend(kws, 'Location', 'Best')
%
% SEE ALSO:
%   readEclipseOutputFileUnFmt, readEclipseRestartUnFmt.

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

[names, nInxOrder, nInx] = unique(nameList);
[kwrds, kInxOrder, kInx] = unique(kwrdList);
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

smry.WGNAMES   = names;
smry.KEYWORDS  = kwrds;
smry.UNITS     = smspec.UNITS.values;
smry.STARTDAT  = smspec.STARTDAT.values([3 2 1])';
smry.nInx      = nInx;
smry.nInxOrder = nInxOrder;
smry.kInx      = kInx;
smry.kInxOrder = kInxOrder;
smry.data      = data(:, 1:curStep);
smry.ministep  = ministeps(1:curStep);

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
    [pick, nm_order] = string_lookup(smry.WGNAMES, smry.nInxOrder, nm);
    rInxN = pick(smry.nInx);
else
    rInxN = true(nlist, 1);
end
if ~isempty(kw)
    if ischar(kw), kw = { kw }; end
    [pick, kw_order] = string_lookup(smry.KEYWORDS, smry.kInxOrder, kw);
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

function [pick, result_order] = string_lookup(list, list_order, search)
ncol = numel(search);
[i, j] = blockDiagIndex(numel(list), ncol);
match = strcmp(reshape(list(i)  , [], ncol), ...
               reshape(search(j), [], ncol));
pick = any(match, 2);
[row, ~] = find(match);
[~, reverse] = sort(list_order(row));
result_order = zeros(size(reverse));
result_order(reverse) = 1 : numel(reverse);
end
