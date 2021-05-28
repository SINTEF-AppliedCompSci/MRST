function [smry, smspec] = readSummaryLocal(prefix, keyWords)
%Undocumented Utility Function

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

nameList = smspec.WGNAMES.values;
kwrdList = smspec.KEYWORDS.values;

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
%estimate number of ministeps
d = dir(smry_file);
estNum = floor(d.bytes/(4*numFull));

data = zeros(nlist, estNum);
ministeps = false(1, estNum);

[fid, msg] = fopen(smry_file, 'r', 'ieee-be');
if fid < 0, error([smry_file, ': ', msg]); end

fprintf(['Reading info from roughly ', num2str(estNum), ' ministeps...\n'])
nums = 0;
while ~feof(fid)
    [name, field] = readFieldUnFmt(fid);
    if strcmp(name, 'MINISTEP')
        ministep = field.values;
        if (50*ministep/estNum) > nums
            fprintf('*');nums=nums+1;
        end
    end
    if strcmp(name, 'PARAMS')
        data(:, ministep+1) = field.values(rowInx);
        ministeps(ministep+1) = true;
    end
end
fclose(fid);
fprintf(['\nActual number of ministeps: ', num2str(ministep+1), '\n']);

data = data(:, ministeps);

smry.WGNAMES  = names;
smry.KEYWORDS = kwrds;
smry.nInx     = nInx;
smry.kInx     = kInx;
smry.data     = data;

smry = addSmryFuncs(smry);

%--------------------------------------------------------------------------

function smry = addSmryFuncs(smry)
smry.get    = @(nm,kw,ms)getData(smry, nm, kw, ms);
smry.getInx = @(nm,kw)getRowInx(smry,nm,kw);
smry.getNms = @(kw)getNames(smry,kw);
smry.getKws = @(nm)getKeywords(smry, nm);


function nms = getNames(smry,kw)
rInx   = getRowInx(smry, [], kw);
nmPos  = unique(smry.nInx(rInx));
nms    = smry.WGNAMES(nmPos);

function kws = getKeywords(smry,nm)
rInx   = getRowInx(smry, nm, []);
kwPos  = unique(smry.kInx(rInx));
kws    = smry.KEYWORDS(kwPos);

function s = getData(smry, nm, kw, ms)
rInx = getRowInx(smry, nm, kw);
s = smry.data(rInx, ms);

function rInx = getRowInx(smry, nm, kw)
nlist = numel(smry.kInx);
if ~isempty(nm)
    i = find(strcmp(smry.WGNAMES, nm));
    rInxN = (smry.nInx==i);
else
    rInxN = true(nlist, 1);
end
if ~isempty(kw)
    j = find(strcmp(smry.KEYWORDS, kw));
    rInxK = (smry.kInx==j);
else
    rInxK = true(nlist, 1);
end
rInx = and(rInxN, rInxK);



