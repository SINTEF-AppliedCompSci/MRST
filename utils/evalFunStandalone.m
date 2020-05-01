function warns = evalFunStandalone(fn, args, startupPath, varargin)
% Utility for performing stand-alone function evaluation in new matlab session
% Only intended to be used with system-calls from evalFunWrapper

% For debugging call evalFunWrapper with option exitWhenDone = false, and see 
% output warnings in 'ans'

opt = struct('moduleList',        {{}}, ...
             'pathList',          {{}}, ...
             'progressFile',        '', ...
             'saveWarnings',     false);
         
% All input should be strings / cellstrings          
% args come in pairs : {class(a), 'a', class(b), 'b', ...}

% setup mrst
run(fullfile(startupPath, 'startup.m'));
opt = merge_options(opt, varargin{:});

% add some default modules if list is empty
if isempty(opt.moduleList)
    opt.moduleList = {'optimization', 'ad-blackoil', 'ad-props', 'ad-core', 'deckformat', 'incomp', 'diagnostics'};
end
mrstModule('add', opt.moduleList{:});
if ~isempty(opt.pathList)
    addpath(opt.pathList{:});
end


% convert main function
fn = str2func(fn);

warns = {};
% convert args
[args, flags] = cellfun(@(s,c)convertFromChar(s,c), args(2:2:end), args(1:2:end), 'UniformOutput', false);
if ~all([flags{:}])
    warns = {'Argument conversion failed for one or more inputs'};
end

try
    feval(fn, args{:})
catch me
    warning('Function evaluation failed')
    warns = [warns, {me}];
end

if ~isempty(opt.progressFile)
    if isfile(opt.progressFile)
        delete(opt.progressFile);
    else
        warning('Progress file not found');
        warns = [warns, {'Progress file not found'}];
    end
end

if opt.saveWarnings && ~isempty(warns)
    wnm = ['warings', datestr(now,'yyyymmdd_HH_MM_SS_FFF')];
    pth = fileparts(mfilename('fullpath'));
    save(fullfile(pth, [wnm, '.mat']), 'warns');
end
    
end

%% ------------------------------------------------------------------------
function [s, flag] = convertFromChar(s, cl)
flag = true;
if ~strcmp(cl, 'char')
    if strcmp(cl, 'double')
        s = str2num(s);                                             %#ok
    elseif strcmp(cl, 'logical')
        s = logical(str2num(s));                                    %#ok
    elseif strcmp(cl, 'function_handle')
        s = str2func(s);
    else
        if ~iscellstr(s)                                            %#ok
            flag = false;
            warning(['Argument conversion failed, got class: ', cl])
        end
    end
end
end

