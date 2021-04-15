function writeSchedule(fn, dirname, SCHEDULE, varargin)
% Write schedule file/section to file. 
%
% SYNOPSIS:
%   writeSchedule(fn, dirname, SCHEDULE)
%   writeSchedule(fn, dirname, SCHEDULE, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function writes the contents of deck-type SCHEDULE-structure to file
%
% REQUIRED PARAMETERS:
%   fn       - Either file-name or file identifier. 
%   dirname  - name of directory in which to write data.  
%   SCHEDULE - deck-type schedule-structure
%
% OPTIONAL PARAMETERS:
%   'writeInclude' - If true and fn is a file identifier, schedule will be 
%                    written to a seperate file (other than fn). 
%                    Default value: false
%   'includeName'  - Name of include file. Default value: SCHEDULE.TXT
%   'onlyWells'    - If true, only write well-specific fields. Default
%                    value: false
%   'fields'       - Specify what fields to output. Default: all fields
%   'formats'      - Specify formats passed to fprintf (see below for details)
%
% SEE ALSO:
%   convertScheduleToDeck, model2Deck, writeDeck

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

opt = struct('writeInclude', false, ...
             'includeName',   '', ...
             'onlyWells',  false, ...
             'fields',      {{}}, ...
             'formats',        []);
opt = merge_options(opt, varargin{:});

f = opt.formats;
if isempty(f)
    f = struct('case',          @upper, ...
               'int',            '%8d', ...
               'string',        '%10s', ...
               'double',      '%12.6f', ...
               'sci',         '%15.7e', ...
               'doubleRange', [.1 100]);
end
fncase = f.case;

[writeCOMPDAT, writeWELSPECS, writeWCONPROD, writeWCONHIST, writeWCONINJE, writeTSTEP] = deal(true);
if ~isempty(opt.fields)
    opt.onlyWells = true;
    writeWELSPECS = ismember('WELSPECS', opt.fields);
    writeCOMPDAT  = ismember('COMPDAT',  opt.fields);
    writeWCONPROD = ismember('WCONPROD', opt.fields);
    writeWCONHIST = ismember('WCONHIST', opt.fields);
    writeWCONINJE = ismember('WCONINJE', opt.fields);
    writeTSTEP    = ismember('TSTEP',    opt.fields);
end
doCloseFile = false;
if isa(fn, 'char')  % assume file-name
    doCloseFile = true;
    [fn, msg] = fopen(fullfile(dirname, fn), 'w+');
    if fn < 0
        error(msgid('Open:Failure'), 'Failed to open ''%s'': %s\n', fn, msg);
    end
    opt.writeInclude = false;
end

if opt.writeInclude
    if isempty(opt.includeName)
        opt.includeName = fncase('schedule.txt');
    end
    fid_inc = fopen(fullfile(dirname, opt.includeName), 'wt');
else
    fid_inc = fn;
end

if isfield(SCHEDULE,'RPTSCHED') && ~opt.onlyWells
    % only needed to get restart output
    fprintf(fn,'RPTSCHED\n');
    fprintf(fn,'PRES SGAS RS WELLS\n');
    fprintf(fn,'/\n');
    fprintf(fn,'RPTRST\n');
    fprintf(fn,'BASIC=1\n');
    fprintf(fn,'/\n');
end

if isfield(SCHEDULE, 'control')
    if opt.writeInclude
        fprintf(fn,'%s\n','INCLUDE');
        fprintf(fn,'%s/', opt.includeName);
        fprintf(fn,'\n\n');
    end
    
    control = SCHEDULE.control;
    for cstep = 1:numel(control)
        wspecs = SCHEDULE.control(cstep).WELSPECS;
        if ~isempty(wspecs) && writeWELSPECS
            fprintf(fid_inc,'%s\n',upper('welspecs'));
            wspecs = replace_default(wspecs(:, 1:13)).';
            % WELSPECS        1         2        3      4      5       6         7      8          9       10       11       12      13
            %                nm        grp       I      J      dep    pph      drain  infl      shutin    cr-fl     pvtn    dens     fip    
            fmt = getFmtStr(f.string, f.string, f.int, f.int, f.sci, f.string, f.int, f.string, f.string, f.string, f.int, f.string, f.int, '/');
            fprintf(fid_inc, fmt, wspecs{:});
            fprintf(fid_inc, '/\n\n');
        end
        
        cdat = SCHEDULE.control(cstep).COMPDAT;
        if ~isempty(cdat) && writeCOMPDAT
            fprintf(fid_inc,'%s\n',upper('compdat'));
            cdat = replace_default(cdat(:, 1:14)).';
            fmts = {f.string, f.int, f.int, f.int, f.int, f.string, f.int, f.sci, f.sci, f.sci, f.sci, f.string, f.string, f.double};
            % all numeric values should be positive, otherwise they are
            % defaulted -> set to -1 to be clear
            [cdat, ix] = replace_negative(cdat);
            fmts(ix(:,1)) = {f.int};
            % COMPDAT           1      2       3      4      5       6        7       8      9    10     11     12        13       14
            %                   nm     I       J      K1     K2      flag    satnum   wi     d    Kh    skin   D-fac     dir       r0
            fmt = getFmtStr(fmts{:}, '/');
            fprintf(fid_inc, fmt, cdat{:});
            fprintf(fid_inc, '/\n\n');
        end
        
        wconinje = SCHEDULE.control(cstep).WCONINJE;
        if ~isempty(wconinje) && writeWCONINJE
            wconinje = replace_default(wconinje(:, 1:14)).';
            fprintf(fid_inc,'%s\n',upper('wconinje'));
            % WCONINJE      1          2        3           4        5      6      7      8     9      10     11      12    13     14    
            %flds =        {'name',  'type',   'flag',   'cntr',   'rate','resv','bhp', 'thp', 'vfp', 'rv',  'th',  'ro',  'rw',  'rg'}; 
            fmt = getFmtStr(f.string, f.string, f.string, f.string, f.sci, f.sci, f.sci, f.sci, f.int, f.sci, f.sci, f.sci, f.sci, f.sci, '/');
            %writeExplain(fid_inc, flds, fmt);
            s = sprintf(fmt, wconinje{:});
            fprintf(fid_inc, f.string, regexprep(s, 'Inf|NaN', '1*', 'ignorecase'));
            fprintf(fid_inc, '/\n\n');
        end
        
        wconprod = SCHEDULE.control(cstep).WCONPROD;
        if ~isempty(wconprod) && writeWCONPROD
            wconprod = replace_default(wconprod(:, 1:12)).';
            fprintf(fid_inc,'%s\n',upper('wconprod'));
            % WCONPROD       1        2          3           4    5      6      7      8      9      10     11     12          
            %                nm       flag      cntr       orat  wrat   grat   lrat   resv   bhp    thp    vfp    lift
            fmt = getFmtStr(f.string, f.string, f.string, f.sci, f.sci, f.sci, f.sci, f.sci, f.sci, f.sci, f.int, f.int, '/');
            s = sprintf(fmt, wconprod{:});
            fprintf(fid_inc, f.string, regexprep(s, 'Inf|NaN', '1*', 'ignorecase'));
            fprintf(fid_inc, '/\n\n');
        end
        wconhist = SCHEDULE.control(cstep).WCONHIST;
        if ~isempty(wconhist) && writeWCONHIST
            wconhist = replace_default(wconhist(:, 1:10)).';
            fprintf(fid_inc,'%s\n',upper('wconhist'));
            % WCONPROD       1        2          3           4    5      6      7      8        9      10              
            %                nm       flag      cntr       orat  wrat   grat   vfptab  alqwell  thp    bhp  
            fmt = getFmtStr(f.string, f.string, f.string, f.sci, f.sci, f.sci, f.int,   f.int,  f.sci, f.sci,'/');
            s = sprintf(fmt, wconhist{:});
            fprintf(fid_inc, f.string, regexprep(s, 'Inf|NaN', '1*', 'ignorecase'));
            fprintf(fid_inc, '/\n\n');
        end
        
        
        if writeTSTEP
            ind = SCHEDULE.step.control == cstep;
            fprintf(fid_inc, '%s\n', upper('tstep'));
            fprintf(fid_inc, [f.double, '\n'], SCHEDULE.step.val(ind));
            fprintf(fid_inc, '/\n');
        end
    end
    if opt.writeInclude
        fclose(fid_inc);
    end
end
if doCloseFile
    fclose(fn);
end
end

%--------------------------------------------------------------------------

function t = replace_default(t)
i = cellfun(@ischar, t(1,:));
t(:,i) = regexprep(t(:,i), 'Default', '1*');
end

%--------------------------------------------------------------------------

function [t, ix] = replace_negative(t)
if ~iscell(t)
    ix = t(1,:) < 0;
    if any(ix)
        t = num2cell(t);
    end
else
    ix = cellfun(@(x)isnumeric(x)&&x<0, t);
end
    if any(ix)
        t(ix) = {-1};
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
function writeExplain(fid, flds, fmt)
% can be used to add comment-line above 
tmp = regexp(fmt, '\%(?<len>\d+)', 'names');
l = arrayfun(@(x)str2double(x.len), tmp);
l(1) = l(1)-2;
n = numel(l);
nfmt = arrayfun(@(x)['%', num2str(x), 'd '], l, 'UniformOutput', false);
nfmt = horzcat(nfmt{:});
sfmt = arrayfun(@(x)['%', num2str(x), 's '], l, 'UniformOutput', false);
sfmt = horzcat(sfmt{:});
fprintf(fid, ['--', nfmt, '\n'], (1:n));
fprintf(fid, ['--', sfmt, '\n'], flds{:});
end
