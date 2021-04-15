function varargout = mrstSettings(verb, varargin)
%Configure persistent settings for MRST for advanced functionality
%
% SYNOPSIS:
%   useMEX = mrstSettings('get', 'useMEX')
%   mrstSettings('set', 'useMEX', useMEX)
%
% REQUIRED PARAMETERS:
%   verb - Must be one of:
%            'get'   - Get the value of a setting. Second input should be a
%                      string, with the name of a valid setting.
%            'set'   - Set the value of a setting. The second input is the
%                      new value.
%            'setup' - Perform a full setup of all possible settings.
%
%
% RETURNS:
%   out - If verb is 'get', the value of the settings. If called with no
%         input arguments, the output will be a struct of all current
%         settings. Other input arguments will result in an error.
%
% EXAMPLE:
%   % Equal to list
%   mrstSettings()
%   % List current settings
%   mrstSettings('list')
%   % Reset and do set up
%   mrstSettings('setup', 'reset-no-check')
%   % Reset settings (after prompt) and do setup
%   mrstSettings('setup', 'reset')
%   % Perform setup wizard
%   mrstSettings('setup')
%   % Set option
%   mrstSettings('set', 'useMEX', true)
%   % Query a setting
%   useMEX = mrstSettings('get', 'useMEX');
%   % Query setting and output extra info
%   [useMEX, d] = mrstSettings('get', 'useMEX');
%   disp(d)
%

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

    persistent SETTINGS
    isDesktop = mrstPlatform('desktop'); % Check for GUI etc
    need_save = false;
    if isempty(SETTINGS)
        SETTINGS = loadSettings();
        if isempty(SETTINGS)
            % Didn't manage to load them - we are restarting
            SETTINGS = firstTimeSetup(SETTINGS, isDesktop);
            need_save = true;
        else
            % We have added some settings between MRST sessions. Won't
            % happen in a release, but might happen when migrating settings
            % files between MRST releases or for developers. We just copy
            % over the default value and print a notification when this
            % happens.
            default = getDefaultMRSTSettings(false);
            f = fieldnames(SETTINGS);
            df = fieldnames(default);
            missing = setdiff(df, f);
            for i = 1:numel(missing)
                nm = missing{i}; 
                v = default.(nm);
                if isstruct(v)
                    fprintf('New setting discovered: %s. Setting default value: %s\n', nm, v.value);
                end
                SETTINGS.(nm) = v; 
            end
            need_save = numel(missing) > 0;
        end
    end
    if nargin == 0 
        if nargout == 1
            verb = '';
        else
            verb = 'list';
        end
    end
    if nargin < 2
        sarg = '';
    else
        sarg = varargin{1};
    end
    % Initialize output
    varargout = cell(1, nargout);
    switch lower(verb)
        case 'setup'
            assert(nargout == 0);
            switch lower(sarg)
                case ''
                    % Redo (but keep current options)
                    SETTINGS = firstTimeSetup(SETTINGS, isDesktop, true);
                case 'reset'
                    % Reset (but prompt first!)
                    doReset = prompt('Reset settings', ...
                        'Reset all MRST settings and set up again? Warning: No way to undo!', isDesktop);
                    if doReset
                        SETTINGS = hardReset();
                        SETTINGS = firstTimeSetup(SETTINGS, isDesktop);
                    else
                        return;
                    end
                case 'reset-no-check'
                    % Reset (no prompt!)
                    SETTINGS = hardReset();
                    SETTINGS = firstTimeSetup(SETTINGS, isDesktop);
                otherwise
                    error('%s not supported', sarg)
            end
            need_save = true;
        case 'set'
            assert(nargin > 2);
            assert(nargout == 0);
            checkSetting(SETTINGS, sarg);
            new = varargin{2};
            if ischar(new) && isstruct(SETTINGS.(sarg))
                if islogical(SETTINGS.(sarg).value)
                    switch lower(new)
                        case {'true', 'on'}
                            new = true;
                        case {'false', 'off'}
                            new = false;
                        case 'toggle'
                            new = ~SETTINGS.(sarg).value;
                        otherwise
                            error('Operation %s not supported for ''set''', new);
                    end
                end
            end
            if ~isempty(regexpi(sarg, '\s*Directory'))                     %#ok
                assert(ischar(new), 'Directory must be a string.');
                if isdir(new)                                              %#ok
                    new = getCanonicalPath(new);
                    SETTINGS.(sarg) = new;
                    need_save = true;
                else
                    warning(['Supplied directory ''', new, ''' is not a directory. ', ...
                             sarg, ' has not been changed']);
                end
            else
                SETTINGS.(sarg).value = new;
                SETTINGS.(sarg).defaulted = false;
                need_save = true;
            end
        case 'get'
            if isempty(sarg)
                varargout{1} = SETTINGS;
            else
                checkSetting(SETTINGS, sarg);
                setting =  SETTINGS.(sarg);
                if ischar(setting)
                    varargout{1} = setting;
                else
                    varargout{1} = setting.value;
                end

                if nargout > 1
                    varargout{2} = setting;
                end
            end
        case 'list'
            % Default
            listSettings(SETTINGS, isDesktop)
        case ''
            varargout{1} = SETTINGS;
        otherwise
            error('Unknown verb ''%s''.', verb);
    end

    if need_save
        storeSettings(SETTINGS);
    end
end
% Setup / wizards / GUI
function settings = firstTimeSetup(settings, isDesktop, doWizard)
    if isempty(settings)
        settings = getDefaultMRSTSettings(true);
    end
    if nargin < 3
        if mrstPlatform('octave')
            doWizard = false;
        else
            doWizard = prompt('MRST settings', ...
                ['MRST has several advanced configuration settings. Would you like to', ...
                ' set these up now? Otherwise, these settings will use reasonable defaults', ...
                ' and can be configured by calling ''mrstSettings'' later.'], isDesktop);
        end
    end
    if doWizard
        if isDesktop
            [names, folders, folderTexts] = getNames(settings);
            for i = 1:numel(names)
                name = names{i};
                v = settings.(name);
                c = v.value;
                l = v.label;
                de = v.description;
                if isempty(de)
                    de = l;
                end
                helpstr = sprintf('%s:\n%s\n', l, de);

                if v.defaulted
                    s = 'default';
                else
                    s = 'current';
                end
                if c
                    tmp = 'Enabled';
                else
                    tmp = 'Disabled';
                end
                defstr = sprintf('Keep %s (%s)', s, tmp);
                result = triplePrompt(name, helpstr, {true, false, nan}, {'Enable', 'Disable', defstr}, defstr);
                if ~isnan(result)
                    settings.(name).value = result;
                    settings.(name).defaulted = false;
                end
            end
            for i = 1:numel(folders)
                name = folders{i};
                helpstr = sprintf('Choose current %s:\n%s', name, folderTexts{i});
                result = triplePrompt(name, helpstr, {true, false}, {'Select new', 'Keep current'}, 'Keep current');
                if ~isnan(result) && result
                    current = settings.(name);
                    selpath = uigetdir(current, sprintf('Set new MRST directory: %s', name));
                    settings.(name) = selpath;
                end
            end
        else
            % Use simple version
            mrstSettings list
        end
    end
end

function [names, folders, folderText] = getNames(settings)
    names = fieldnames(settings);
    folders = {'outputDirectory', 'dataDirectory'};
    folderText = {['Directory where persistent MRST output is stored, for', ...
                   ' example when running simulations that output results to disk'], ...
                   'Directory where MRST stores downloaded data sets (geological models, grids, etc) for simulations.'};
    names = setdiff(names, folders);
end

function result = triplePrompt(head, txt, choices, labels, default)
    answer = questdlg(txt, head, labels{:}, default);
    if isempty(answer)
        result = nan;
    else
        result = choices{strcmp(labels, answer)};
    end
end

function wasYes = prompt(head, txt, isDesktop)
    if isDesktop
        answer = questdlg(txt, head, 'Yes', 'No', 'No');
        % Handle response
        switch lower(answer)
            case 'yes'
                wasYes = true;
            otherwise
                wasYes = false;
        end
    else
        s = sprintf('%s: %s\nyes/no [no]: ', head, txt);
        v = input(s, 's');
        if isempty(v)
            wasYes = false;
        else
            switch lower(v)
                case {'y', 'yes'}
                    wasYes = true;
                otherwise
                    wasYes = false;
            end
        end
    end
end

% Listing
function listSettings(settings, isDesktop)
    [names, fldrs] = getNames(settings);
    n = 0;
    for i = 1:numel(names)
        name = names{i};
        v = settings.(name);
        n = max(n, numel(v.label));
    end
    fprintf('Current MRST settings (persistent between sessions):\n\n');
    nchar = fprintf('      Setting | Current value  | %-*s | Default\n', n, 'Description');
    fprintf('%s\n', repmat('-', nchar + 1, 1));
    for i = 1:numel(names)
        name = names{i};
        v = settings.(name);
        c = v.value;
        d = v.default;
        l = v.label;
        if isDesktop
            current = formatToggle(name, c);
        else
            current = formatOption(c);
        end
        fprintf(' %12s |  %12s  | %-*s | %5s\n', name, current, n, l, formatOption(d));
    end
    fprintf('\nConfigured directories:\n');
    for i = 1:numel(fldrs)
        f = fldrs{i};
        fprintf('\t%16s -> %s\n', f, settings.(f));
    end
end

function s = formatToggle(name, opt)
    if islogical(opt)
        if opt
            a = 'on';
        else
            a = 'off';
        end
        b = 'toggle';
        s = sprintf('%3s [<a href="matlab:mrstSettings(''set'',''%s'',''%s'');mrstSettings(''list'')">toggle</a>]', a, name, b);
    else
        error('Only logicals supported')
    end
end

function s = formatOption(opt)
    if islogical(opt)
        if opt
            s = 'on ';
        else
            s = 'off';
        end
    else
        error('Only logicals supported')
    end
end

% Storage / retrieval functions
function settings = loadSettings()
    pth = getSettingsPath();
    present = exist(pth, 'file');
    if present
        settings = load(pth);
    else
        settings = [];
    end
end

function settings = hardReset()
    settings = getDefaultMRSTSettings(true);
    storeSettings(settings);
end

function storeSettings(settings)
    pth = getSettingsPath();
    save(pth, '-struct', 'settings');
end

function pth = getSettingsPath()
    if exist('OCTAVE_VERSION', 'builtin') > 0
        % Use different set of settings for Octave in case someone is using
        % the same folder with both Matlab and Octave.
        fn = 'settings_octave.mat';
    else
        fn = 'settings.mat';
    end
    pth = fullfile(ROOTDIR(), fn);
end

function checkSetting(SETTINGS, setting)
    if ~isfield(SETTINGS, setting)
        error('I do not know of the option ''%s''. Note that options are case sensitive', setting);
    end
end

% Setting up defaults 
function opts = getDefaultMRSTSettings(setDefaults)
    % This should run once!
    if nargin == 0
        setDefaults = true;
    end
    out = default_output_dir();
    data = default_data_dir();
    ensure_directory_exists(out)
    ensure_directory_exists(data)

    opts = struct('outputDirectory', out, ...
                  'dataDirectory',   data);
    
    useMEX = checkCompiler(setDefaults);
    allowDL = checkDownload(setDefaults);
    useOpenMP = ~strcmpi(mrstPlatform('os'), 'macos');

    opts = setDefaultSetting(opts, 'useMEX', 'Use MEX-acceleration', ...
        ['Parts of MRST can be accelerated with compiled extensions (MEX-acceleration).', ...
        ' This requires a C/C++ compiler to be set up with Matlab (see ''mex -setup'' for more details). ', ...
        'Enabling this option will let MRST attempt to compile MEX files when needed. Building of such files can still be triggered manually.'], ...
        useMEX);
    opts = setDefaultSetting(opts, 'useOMP', 'Use OpenMP when using MEX-acceleration', ...
        ['Many MEX extensions include OpenMP support for parallelization. Enable this if you have OpenMP installed.'], ...
        useOpenMP);
    opts = setDefaultSetting(opts, 'promptMEX', 'Prompt before building MEX-extensions', ...
        'MRST normally asks before attempting to build MEX extensions. If disabled, building happens automatically without user input.', false);
    opts = setDefaultSetting(opts, 'allowDL', 'Download files from the internet', ...
        ['MRST may download data sets and other utilities from the internet.', ...
        ' Disabling this means that no files will be downloaded by default. You can still manually activate downloads through optional inputs.'], allowDL);
    opts = setDefaultSetting(opts, 'promptDL', 'Prompt before downloading files', ...
        ['MRST will ask for permission before downloading files.', ...
        ' You can disable these prompts (recommended if you are on an unmetered high-speed connection)'], true);
end

function opt = getDefaultSetting(label, description, default)
    opt = struct('description', description, ...   % Long description
                 'label',       label, ...         % Short label
                 'value',       default, ...       % Current value
                 'default',     default, ...       % Default value (on current system)
                 'defaulted',   true);             % Has not been explicitly set yet
end

function opts = setDefaultSetting(opts, name, varargin)
    opts.(name) = getDefaultSetting(varargin{:});
end

function allowDL = checkDownload(setDefaults)
    testurl = 'https://www.sintef.no/contentassets/efe341cda156406aa06bc8e6a149aa26/DL_TEST.zip';
    if setDefaults
        out = default_output_dir();
        try
            unzip(testurl, out);
            allowDL = exist(fullfile(out, 'DL_TEST.txt'), 'file') > 0;
        catch
            allowDL = false;
        end
    else
        allowDL = true;
    end
end

function useMEX = checkCompiler(setDefaults)
    useMEX = false;
    if setDefaults
        try
            compilers = mex.getCompilerConfigurations();
            for i = 1:numel(compilers)
                l = compilers(i).Language;
                useMEX = useMEX || strcmp(l, 'C') || strcmp(l, 'C++');
                if useMEX
                    % We found a C/C++ compiler. Abort.
                    break;
                end
            end
        catch
            % useMEX should be false
        end
    end
end

function ddir = default_data_dir()
    ddir = fullfile(ROOTDIR, 'examples', 'data');
end

function odir = default_output_dir()
    odir = fullfile(ROOTDIR, 'output');
end

function ensure_directory_exists(ddir)
    if ~isdir(ddir)                                                        %#ok
        [ok, msg, id] = mkdir(ddir);

        if ~ok
            error(id, 'Failed to create directory ''%s'': %s', ...
                  ddir, msg);
        end
    end
end
