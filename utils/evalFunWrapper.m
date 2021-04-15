function [status, str] = evalFunWrapper(fn, args, varargin)
% Utility to launch (mrst-)function evaluation in seperate matlab session 
%
% SYNOPSIS:
%   evalFunWrapper(fn, args, varargin)
%
% DESCRIPTION:
%   This function launches a seperate matlab session and perfoms the evaluation 
%       feval(fn, args{:})
%   by invoking evalFunStandalone
%
% PARAMETERS:
%   fn      - function handle of non-anonymous function. Corresponding
%             m-file must be on matlab/mrst-path (see option moduleList).  
%   args    - list of function arguments. Supported arguments are strings, 
%             function handles, cellstrings and doubles. Numeric arguments 
%             must be scalar or on row form. 
%
% OPTIONAL PARAMETERS
%  progressFileNm  - If non-empty a file with name progressFileNm will be 
%                    written before launcing new session. File is deleted
%                    when new session is done computing fn. File name
%                    should include full path. 
%  moduleList      - List of mrst-modules that will be added in new session
%                    Default is currently active modules
%
%  pathList        - List of (non-mrst) paths that will be added
%
%  background      - Option indicating if new session is run in background.
%                    (default true)
%
%  singleCompThread - Option indicating whether new session should be run
%                     using single computational thread (default true)
%
%  matlabBinary    - Command to start external Matlab-session. Default is  
%                    fullfile(matlabroot(), 'bin', 'matlab'), e.g., same
%                    Matlab-version as the one currently running.
%
%  matlabOpts      - String of matlab startup options. Default: see below
% 
%  exitWhenDone    - Exit/quit new session when computation is done 
%                    (default true)
% RETURNS:
%   status         - when run in background, status will be 0. 
%
% EXAMPLE:
%   evalFunWrapper(@arrayfun, {@why, 1:20}, 'exitWhenDone', false)
%
% SEE ALSO
% evalFunStandalone

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
opt = struct('progressFileNm',       '', ...
             'moduleList',         {{}}, ...
             'pathList',           {{}}, ...
             'background',         true, ...
             'singleCompThread',   true, ...
             'matlabBinary',         '', ...
             'matlabOpts',           '', ...
             'exitWhenDone',       true);         
opt = merge_options(opt, varargin{:});

if isempty(opt.moduleList)
    opt.moduleList = mrstModule();
end
       
assert(iscell(args), 'Function arguments must come in a list (cell)');

% warn if function is anonymous as this is not recommended
if isa(fn, 'function_handle')
    checkFunction(fn)
end

% write progress file 
fnm = opt.progressFileNm;
if ~isempty(fnm)
    flag = writeProgressFile(fnm);
    if ~flag
        return;
    end
end

% Set matlab binary
mbin = opt.matlabBinary;
if isempty(mbin)
    mbin = fullfile(matlabroot(), 'bin', 'matlab');
    % Support folder names with spaces (e.g., 'Program Files' on win)
    mbin = strcat('"', mbin, '"');
    % If 'matlabBinary' is defaulted on a Windows-system, force background-option 
    % to be false to avoid opening of additional cmd.exe-window(s). To enable 
    % backgrounds-calls on Windows-systems (e.g., such that main command-window 
    % will wait for external session to finnish) setting option 'matlabBinary'
    % to 'matlab' can be used, allthough this will launch default version 
    % (as of R2020b/Windows 10 Enterprise/Pro).
    if ispc
        opt.background = false;
    end
end

% handle matlab options
mopt = opt.matlabOpts;
if isempty(mopt)
    if ispc
        mopt = '-nodisplay -nosplash -nodesktop -minimize';
    else
        mopt = '-nojvm -nodisplay -nosplash';
    end
end

% startup folder option (location of evalFunStandalone)
sdopt = ['-sd "', fileparts(mfilename('fullpath')), '"'];

% single comp thread opt
stopt = '';
if opt.singleCompThread
    stopt = '-singleCompThread';
end

% main function call
expr = getExpression(fn, args, opt);

% run in background opt
bgopt = '&';
if ~opt.background
    bgopt = '';
end

% Add extra space to be sure
[mbin, mopt, sdopt, stopt, expr, bgopt] = addBlanks(mbin, mopt, sdopt, stopt, expr, bgopt);

% final string for pass to system
str = [mbin, mopt, sdopt, stopt, '-r', expr, bgopt];

status = system(str);
end

%% -------------------------------------------------------------------------
function expr = getExpression(fn, args, opt)
n  = numel(args);
cl = cellfun(@class, args, 'UniformOutput', false);

ix = find(strcmp('function_handle', cl));
for k = 1:numel(ix)
    checkFunction(args{ix(k)});
end
% accept fn char or function-handle
fn = convert2char(fn, class(fn));

argstr = cell(1, 2*n);
% include argument classes
argstr(1:2:end) = cl;
% convert arguments to strings
argstr(2:2:end) = cellfun(@(a,c)convert2char(a,c), args, cl, 'UniformOutput', false);

% opt-type input
optstr = {};
if ~isempty(opt.progressFileNm)
    optstr = {'progressFile', opt.progressFileNm};
end
if ~isempty(opt.moduleList)
    optstr = [optstr, {'moduleList', opt.moduleList}]; % module-list all strings
end
if ~isempty(opt.pathList)
    optstr = [optstr, {'pathList', opt.pathList}]; % path-list all strings
end

% function input to:
% evalFunStandalone(fn, args, startupPath, varargin)
% args in {}
inpstr = commacat( quote([{fn}, {argstr}, {ROOTDIR}, optstr]) );

%expr = ['"cd(', quote(curpath), '); evalFunStandalone(', inpstr, '); quit"'];
if opt.exitWhenDone
    expr = ['"evalFunStandalone(', inpstr, '); quit"'];
else
    expr = ['"evalFunStandalone(', inpstr, ');"'];
end
end

%%
function s = convert2char(s, cl)
if ~strcmp(cl, 'char')
    if any(strcmp(cl, {'double', 'logical'}))
        s = num2str(s);
        assert(size(s,1)<=1, 'Only row vector numerical arguments are supported')
    elseif strcmp(cl, 'function_handle')
        s = func2str(s);
    else
        assert(iscellstr(s), 'Arguments expected to be of class, ''char'', ''double'', ''function_handle'' or cell of strings'); %#ok
    end
end
end

%% 
function s = quote(c)
if ischar(c)
    s = ['''', c, ''''];
else
    s = cellfun(@quote, c, 'UniformOutput', false);
end
end

%%
function s = commacat(c)
if isempty(c)
    s = '';
else
    ix = cellfun(@iscell, c);
    if any(ix)
        c(ix) = cellfun(@(x)['{', commacat(x), '}'], c(ix), 'UniformOutput', false);
    end
    s = sprintf([repmat('%s,', 1, numel(c)-1), '%s'], c{:});
end
end

%%
function flag = writeProgressFile(fnm)
    if isfile(fnm)
        warning('Attempt to run simulation that is allready running, aborting');
        flag = false;
    else
        progress = clock;
        save(fnm, 'progress');
        flag = true;
    end
end

%%
function [varargout] = addBlanks(varargin)
varargout = varargin;
for k = 1:numel(varargout)
    if ~isempty(varargout{k})
        varargout{k} = [' ', varargout{k}, ' '];
    end
end
end

%% 
function checkFunction(fn)
fnh = functions(fn);
if strcmp(fnh.type, 'anonymous'), ...
    warning('Passing anonymous functions to external sessions is not recommended: \n   %s  ', fnh.function);
end
end


