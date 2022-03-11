function out = mrstPlatform(arg)
%Retrieve Platform Dependent Settings for Current Session
%
% SYNOPSIS:
%   out = mrstPlatform
%   out = mrstPlatform(setting)
%
% PARAMETERS:
%   setting - Name of a specific setting.  If omitted, the return value
%             is a structure of all platform dependent settings.  Character
%             vector or, if available, a scalar string.
%
%             Note: Setting names may be shortened to a unique prefix and
%             are matched case insensitively for typing convenience, but
%             scripts/functions should use the full, unambigiuous setting
%             names listed below.
%
% RETURNS:
%   out - Either all features as a struct, or a specific setting, dependent
%         if a input argument is given. The possible features are:
%               'platform'  - Name of the platform (MATLAB or Octave)
%               'version'   - Version string from ver() output.
%               'major'     - Numerical value of major version
%               'minor'     - Numerical value of minor version
%               'octave'    - Boolean indicating if platform is Octave
%               'matlab'    - Boolean indicating if platform is MATLAB
%               'os'        - 'linux/unix', 'windows' or 'macos'
%               'desktop'   - If MATLAB/Octave is running in Desktop mode
%               'gui',      - If GUIs are supported (= QT backend for Oct)
%               'jvm',      - Availability of Java
%               'richtext'  - Command line supports rich text.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    persistent DATA

    if isempty(DATA)
        DATA = run_configuration_setup();
    end

    if nargin == 0
        out = DATA;
    else
        out = get_specific_setting(DATA, arg);
    end
end

%--------------------------------------------------------------------------

function config = run_configuration_setup()
   if exist('OCTAVE_VERSION', 'builtin') > 0
      config = octave_setup();
   else
      config = matlab_setup();
   end
end

%--------------------------------------------------------------------------

function setting = get_specific_setting(config, setting)
   if isfield(config, setting)
      setting = config.(setting);
   else
      setting = get_specific_setting_from_prefix(config, setting);
   end
end

%--------------------------------------------------------------------------

function setting = get_specific_setting_from_prefix(config, setting)
   settings = fieldnames(config);
   ix = find(strncmpi(settings, setting, numel(char(setting))));

   if numel(ix) == 1
      setting = config.(settings{ix});

   elseif isempty(ix)
      error('Unsupported setting name ''%s''.  Must be one of\n%s', ...
            setting, stringify_settings_list(settings));
   else
      error(['Ambiguous setting name ''%s''.  Could be one of\n%s', ...
             '\n\nPlease disambiguate'], setting, ...
             stringify_settings_list(settings(ix)));
   end
end

%--------------------------------------------------------------------------

function config = octave_setup()
    v = ver();
    isMain = arrayfun(@(x) strcmpi(x.Name, 'octave'), v);

    % Mostly where GUIs are supported
    hasGUI = strcmpi(graphics_toolkit(), 'qt');

    config = get_config(v(isMain), ...
                        'gui', hasGUI, ...
                        'desktop', isguirunning(), ...
                        'richtext', false);
end

%--------------------------------------------------------------------------

function config = matlab_setup()
    v = ver();
    isMain = arrayfun(@(x) strcmpi(x.Name, 'matlab'), v);
    isDesktop = usejava('desktop');

    config = get_config(v(isMain), ...
                        'desktop', isDesktop, ...
                        'richtext', isDesktop, ...
                        'gui', isDesktop);
end

%--------------------------------------------------------------------------

function s = get_config(ver_struct, varargin)
    name = ver_struct(1).Name;
    ver_char = ver_struct(1).Version;

    v = parse_version_numbers(ver_char);
    s = struct('platform',      name, ...
               'version',       ver_char, ...
               'major',         v(1), ...
               'minor',         v(2), ...
               'octave',        strcmpi(name, 'octave'), ...
               'matlab',        strcmpi(name, 'matlab'), ...
               'os',            infer_operating_system(), ...
               'desktop',       false, ...
               'gui',           false, ...
               'jvm',           usejava('jvm'), ...
               'richtext',      false);

   s = merge_options(s, varargin{:});
end

%--------------------------------------------------------------------------

function s = stringify_settings_list(settings)
   settings = sort(settings);
   s = sprintf('\n  * ''%s''', settings{:});
end

%--------------------------------------------------------------------------

function os = infer_operating_system()
    if ispc()
        os = 'windows';
    elseif ismac()
        os = 'macos';
    else
        % Must be Linux or Unix since we already
        % excluded mac os.
        os = 'linux/unix';
    end
end

%--------------------------------------------------------------------------

function v = parse_version_numbers(ver_char)
    v = sscanf(ver_char, '%d.%d.%d');
    if numel(v) < 3, v(3) = 0; end
end
