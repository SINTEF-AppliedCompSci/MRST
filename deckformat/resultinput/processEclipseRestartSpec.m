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
    % check and repair pointers if # > 2^31 (only if unified)
    if strcmp(type, 'unified')
        rsspec = pointersToDouble(rsspec);
    end
    ix = [nf; numel(rsspec.NAME.values)+1];
    for k = 1:nSteps
        subix = ix(k) : (ix(k+1)-1);
        nms   = rsspec.NAME.values( subix );
        % Reduce index-set if restartAmount is set to 'default'
        if strcmp(restartAmount, 'default')
            subix = subix( fieldIx(nms, {}) );
        end
        %keywords{k} = cellfun(@fixVarName, rsspec.NAME.values( subix ), ...
        %                      'UniformOutput', false);
        keywords{k} = rsspec.NAME.values( subix );
        pointers{k} = rsspec.POINTER.values( subix );
        prec{k}     = mapPrec(rsspec.TYPE.values( subix ));
        num{k}      = rsspec.NUMBER.values( subix );
    end
end 

spec = struct('time', t, 'date', date, 'unit', unit, ...
               'type', type, 'keywords', {keywords}, ...
               'pointers', {pointers}, 'prec', {prec}, ...
               'num', {num});


% Add file-names in case multiple output
if strcmp(type, 'multiple')
    f    = dir([prefix, '.X*']);
    fnms = {f.name};
    pattern = '.X[0-9][0-9][0-9][0-9]';
    ix = ~cellfun(@isempty, regexp(fnms, pattern));
    fnms = sort(fnms(ix));
    fnms = cellfun(@(fn)fullfile(fileparts(prefix), fn), fnms, ...
            'UniformOutput', false);
    assert(numel(fnms)==numel(spec.time), ...
        'Unable to match multiple restart-files to RSSPEC');
    spec.fnames = fnms;
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
    fldnms = {'ACAQ',       ...
              'DOUBHEAD',   ...
              'ENDLGR',     ...
              'FLR',        ...
              'IAAQ',       ...
              'ICAQ',       ...
              'ICON',       ...
              'INTEHEAD',   ...
              'IWEL',       ...
              'LOGIHEAD',   ...
              'LGR',        ...
              'PADS',       ...
              'POLYMAX',    ...
              'POLYMER',    ...
              'PRESSURE',   ...
              'RS',         ...
              'RV',         ...
              'XMF',        ...
              'YMF',        ...
              'ZMF',        ...
              'VGAS',       ...
              'VOIL',       ...
              'VISC',       ...
              'DEN',        ...
              'SCON',       ...
              'SEQNUM',     ...
              'SGAS',       ...
              'SOIL',       ...
              'SWAT',       ...
              'SWEL',       ...
              'XAAQ',       ...
              'XCON',       ...
              'XWEL',       ...
              'ZWEL'};
end

%--------------------------------------------------------------------------

function ix = fieldIx(nms, additionalFields)
    inc = [defaultFields, additionalFields];
    ix  = false(numel(nms), 1);
    
    for k = 1:numel(inc)
        ix = or(ix, strncmp(inc{k}, nms, numel(inc{k})));
    end
end

%--------------------------------------------------------------------------

function rsspec = pointersToDouble(rsspec)
% repair pointer values exceeding single precision maximum
    jmp = diff(rsspec.POINTER.values)<0;
    if any(jmp)
        fx = [0 ; cumsum(jmp)];
        rsspec.POINTER.values = rsspec.POINTER.values + fx*2^31;
    end
end

    
              
              

        
    
