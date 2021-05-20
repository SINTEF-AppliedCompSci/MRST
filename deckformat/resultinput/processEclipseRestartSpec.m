function [spec, specLGR] = processEclipseRestartSpec(prefix, restartAmount)
% Read unformatted (binary) ECLIPSE restart specification file (*RSSPEC) and 
% produce index structure for subsequent flexible/efficient restart reading.  
%
% SYNOPSIS:
%   spec = readEclipseOutputFileFmt(filename, outputAmount)
%
% PARAMETERS:
%   rsspec   - Struct as read from readEclipseOutputFile* or filename 
%              (string) of file ECLIPSE specifications (with or without
%              file extension).
%
%   restartAmount - ('default' , 'all'). Value 'default' builds indexes to
%              only a subset of the restart fields (typically the most 
%              interesting ones, see helper-function 'defaultFields'). 
%
%
% RETURNS:
%   spec    - data structure containing keywords info/pointers for reading 
%            (parts of) binary restart file(s)
%   specLGR - array of structs containing same as spec but for each local 
%             grid 
%
% SEE ALSO:
%   `readEclipseOutputFileUnFmt`, `readRestartUnFmt`.

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

if nargin < 2
    restartAmount = 'default';
end

[pth, nm, ext] = fileparts(prefix);
if isempty(ext), ext = '.RSSPEC'; end

if ~strcmp(ext, '.RSSPEC')
   error(['Unexpected file format: ', ext])
end

filename = fullfile(pth, [nm, ext]);
rsspec = readEclipseOutputFileUnFmt(filename);

% get unit
ih = rsspec.INTEHEAD.values;
unit = {'metric', 'field', 'lab'};
unit = unit{ih(3)};

% time/date
[t, d] = deal(rsspec.TIME.values, rsspec.ITIME.values);
nSteps = numel(t);
date = reshape(d, [numel(d)/nSteps, nSteps])';
% report number (needed for multiple restart)
repNum = date(:,1);
date = date(:, 2:4);


% determine restart type (unified or multiple)
firstField = rsspec.NAME.values{1};
if strcmp(firstField, 'SEQNUM')
    type = 'unified';
elseif strcmp(firstField, 'INTEHEAD')
    type = 'multiple';
else
    type = 'unknown';
    warning(['Unexpected format of file: ', fn, '\n'])
end

% keywords/pointers/precision
nf     = find(strcmp(firstField, rsspec.NAME.values));
if ~(numel(nf)==nSteps)
    warning('Unable to output keyword-info');
    [keywords, pointers, prec, num] = deal({});
else
    [keywords, pointers, prec, num] = deal(cell(1, nSteps));

    ix = [nf; numel(rsspec.NAME.values)+1];
    for k = 1:nSteps
        subix = ix(k) : (ix(k+1)-1);
        nms   = rsspec.NAME.values( subix );
        % Reduce index-set if restartAmount is set to 'default'
        if strcmp(restartAmount, 'default')
            subix = subix( fieldIx(nms, {}) );
        end

        keywords{k} = rsspec.NAME.values( subix );

        % Pointer 'B' is high address, 'pointer' is low address.
        pointers{k} = pow2(rsspec.POINTERB.values(subix), 31) + ...
                      rsspec.POINTER.values(subix);

        prec{k}     = mapPrec(rsspec.TYPE.values( subix ));
        num{k}      = rsspec.NUMBER.values( subix );
    end
end 

spec = struct('time', t, 'date', date, 'unit', unit, ...
               'type', type, 'keywords', {keywords}, ...
               'pointers', {pointers}, 'prec', {prec}, ...
               'num', {num}, 'repNum', repNum);


% Add file-names in case multiple output
if strcmp(type, 'multiple')
    f    = dir([prefix, '.X*']);
    fnms = {f.name};
    pattern = '.X[0-9][0-9][0-9][0-9]';
    mtch = regexp(fnms, pattern, 'match');
    ix   =  ~cellfun(@isempty, mtch);
    if ~all(ix)
        [mtch, fnms] = deal(mtch(ix), fnms(ix));
    end
    num  = cellfun(@(c)str2double(c{1}(3:end)), mtch);
    if min(num) == 0
        num = num+1;
    end
    nspec = numel(spec.time);
    if ~issorted(num)
        [num, six] = sort(num);
        fnms = fnms(six);
    end
    fnms = cellfun(@(fn)fullfile(fileparts(prefix), fn), fnms(ix), ...
                   'UniformOutput', false);
    if numel(fnms)==nspec % we  have all steps
        spec.fnames = sort(fnms);
    else % we have some subset, match to repNum
        subix = ismember(spec.repNum, num);
        assert(nnz(subix) == numel(num), ...
             'Unable to match multiple restart-files to RSSPEC');
        spec.fnames = cell(1, nspec);
        spec.fnames(subix) = fnms;
    end
end


% Process possible multiple time appearing ZTRACER 
for k = 1:numel(spec.time)
    % Append tracer name to ZTRACER to uniquify
    ixz = find(strcmp('ZTRACER', spec.keywords{k}));
    for kz = reshape(ixz, 1, [])
        spec.keywords{k}{kz} = ['ZTRACER_', spec.keywords{k}{kz+1}];
    end
end
  
% Process possible multiple appearing ICAQ, SCAQ, ACAQ
flds = {'ICAQNUM', 'SCAQNUM', 'ACAQNUM'};
for k = 1:numel(spec.time)
    for fk = 1:numel(flds)
        ix = find(strcmp(flds{fk}, spec.keywords{k}));
        for fkix = 1:numel(ix)
            spec.keywords{k}{ix(fkix)} = [spec.keywords{k}{ix(fkix)}, '_', num2str(fkix)]; 
            spec.keywords{k}{ix(fkix)+1} = [spec.keywords{k}{ix(fkix)+1}, '_', num2str(fkix)]; 
        end
    end
end

% Process LGR sectiona and output
outputSpecLGR = nargout > 1;
if outputSpecLGR
    nlgr = max(cellfun(@(x)nnz(strcmp('LGR', x)), spec.keywords));
    specLGR = repmat(spec, 1, nlgr);
    lgrNames = {};
    ixLGR    = 0;
end

for k = 1:numel(spec.time) 
    % process LGR-info (remove from spec and move to specLGR if requested)
    ixl = find(strcmp('LGR', spec.keywords{k}));
    if ~isempty(ixl)
        ixr = find(strcmp('ENDLGR', spec.keywords{k}));
        keepix = true(numel(spec.keywords{k}),1);
        for lgrNum = 1:numel(ixl)
            ix = ixl(lgrNum):ixr(lgrNum);
            keepix(ix) = false;
            if outputSpecLGR
                ixLGR   = ixLGR+1;
                curName = rsspec.LGRNAME.values{ixLGR};
                strctIx = find(strcmp(curName, lgrNames));
                if isempty(strctIx)
                    lgrNames = [lgrNames, curName]; %#ok
                    strctIx  = numel(lgrNames);
                end
                for f = {'keywords','prec', 'pointers', 'num'}
                    specLGR(strctIx).(f{1}){k} = specLGR(strctIx).(f{1}){k}(ix);
                end
                specLGR(strctIx).name = curName;
            end
        end
        % only keep non-LGR fields    
        for f = {'keywords','prec', 'pointers', 'num'}
            spec.(f{1}){k} = spec.(f{1}){k}(keepix);
        end
    end
end
end

%--------------------------------------------------------------------------

function types = mapPrec(types)
% Set precision of records for various types (messages (type MESS) have empty records):
    tpstr = {'CHAR',            'DOUB',         'INTE',       'LOGI',       'MESS',    'REAL'};
    pcstr = {'840*uchar=>char', '1000*float64', '1000*int32', '1000*int32', 'blunder', '1000*float32'};
    for k = 1:numel(tpstr)
        mtch = strcmp(types, tpstr{k});
        types(mtch) = pcstr(k);
    end
end

%--------------------------------------------------------------------------

function fldnms = defaultFields()
   fldnms = { ...
      'ACAQ'    , ...
      'DEN'     , ...
      'DLYTIM'  , ...
      'DOUBHEAD', ...
      'DUDF'    , ...
      'DUDW'    , ...
      'FLR'     , ...
      'IAAQ'    , ...
      'IACN'    , ...
      'IACT'    , ...
      'IAQN'    , ...
      'ICAQ'    , ...
      'ICON'    , ...
      'IDWL'    , ...
      'IGPH'    , ...
      'IGRP'    , ...
      'INTEHEAD', ...
      'IUAD'    , ...
      'IUAP'    , ...
      'IUDQ'    , ...
      'IWEL'    , ...
      'IWLS'    , ...
      'LGR'     , ...
      'LOGIHEAD', ...
      'PADS'    , ...
      'PBUB'    , ...
      'POLYMAX' , ...
      'POLYMER' , ...
      'PRESSURE', ...
      'RAQN'    , ...
      'REGDIMS' , ...
      'REGRPT'  , ...
      'RS'      , ...
      'RV'      , ...
      'SAAQ'    , ...
      'SACN'    , ...
      'SACT'    , ...
      'SCAQ'    , ...
      'SCON'    , ...
      'SDWL'    , ...
      'SEQNUM'  , ...
      'SGAS'    , ...
      'SGRP'    , ...
      'SOIL'    , ...
      'SWAT'    , ...
      'SWEL'    , ...
      'VGAS'    , ...
      'VISC'    , ...
      'VOIL'    , ...
      'XAAQ'    , ...
      'XCON'    , ...
      'XGRP'    , ...
      'XMF'     , ...
      'XWEL'    , ...
      'YMF'     , ...
      'ZACN'    , ...
      'ZACT'    , ...
      'ZDWL'    , ...
      'ZGRP'    , ...
      'ZLACT'   , ...
      'ZMF'     , ...
      'ZUDL'    , ...
      'ZUDN'    , ...
      'ZWEL'    , ...
      'ZWLS'    , ...
   };
end

%--------------------------------------------------------------------------

function ix = fieldIx(nms, additionalFields)
    inc = [defaultFields, additionalFields];
    ix  = false(numel(nms), 1);
    
    for k = 1:numel(inc)
        ix = or(ix, strncmp(inc{k}, nms, numel(inc{k})));
    end
end
