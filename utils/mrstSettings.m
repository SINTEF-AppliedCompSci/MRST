function varargout = mrstSettings(verb, varargin)
    persistent SETTINGS
    isDesktop = usejava('desktop'); % Check for GUI etc
    need_save = false;
    if isempty(SETTINGS)
        SETTINGS = loadSettings();
        if isempty(SETTINGS)
            % Didn't manage to load them - we are restarting
            SETTINGS = firstTimeSetup(SETTINGS, isDesktop);
            need_save = true;
        end
    end
    if nargin == 0
        verb = 'list';
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
            if ischar(new) && islogical(SETTINGS.(sarg).value)
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
            SETTINGS.(sarg).value = new;
            SETTINGS.(sarg).defaulted = false;
            need_save = true;
        case 'get'
            if isempty(sarg)
                varargout{1} = SETTINGS;
            else
                checkSetting(SETTINGS, sarg);
                varargout{1} = SETTINGS.(sarg);
            end
        case 'list'
            % Default
            listSettings(SETTINGS, isDesktop)
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
        doWizard = prompt('MRST settings', ...
            ['MRST has several advanced configuration settings. Would you like to', ...
            ' set these up now? Otherwise, these settings will use reasonable defaults', ...
            ' and can be configured by calling ''mrstSettings'' later.'], isDesktop);
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
    fprintf('      Setting | Current value  | Description %s | Default\n', repmat(' ', n - 12, 1));
    fprintf('%s\n', repmat('-', n + 45, 1));
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
        fprintf(' %12s |  %12s  | %s%s | %5s \n', name, current, l, repmat(' ', n - numel(l), 1), formatOption(d));
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
    pth = fullfile(ROOTDIR(), 'settings.mat');
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
    opts = struct('outputDirectory', mrstOutputDirectory(), ...
                  'dataDirectory',   mrstDataDirectory());
              
    useMEX = checkCompiler(setDefaults);
    allowDL = checkDownload(setDefaults);

    opts = setDefaultSetting(opts, 'useMEX', 'Use MEX-acceleration', ...
        ['Parts of MRST can be accelerated with compiled extensions (MEX-acceleration).', ...
        ' This requires a C/C++ compiler to be set up with Matlab (see ''mex -setup'' for more details). ', ...
        'Enabling this option will let MRST attempt to compile MEX files when needed. Building of such files can still be triggered manually.'], ...
        useMEX);
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
        out = mrstOutputDirectory();
        unzip(testurl, out);
        allowDL = exist(fullfile(out, 'DL_TEST.txt'), 'file') > 0;
    else
        allowDL = true;
    end
end

function useMEX = checkCompiler(setDefaults)
    useMEX = false;
    if setDefaults
        compilers = mex.getCompilerConfigurations();
        for i = 1:numel(compilers)
            l = compilers(i).Language;
            useMEX = useMEX || strcmp(l, 'C') || strcmp(l, 'C++');
            if useMEX
                % We found a C/C++ compiler. Abort.
                break;
            end
        end
    end
end
