function out = mrstPlatform(arg)
% Get information about platform-dependent settings for current session
%
% SYNOPSIS:
%   out = mrstPlatform();
%   isOctave = mrstPlatform('octave');
%   isDesktop = mrstPlatform('desktop');
%
% REQUIRED PARAMETERS:
%   arg - May be omitted to get full platform struct. Otherwise, the name
%         of a specific setting to query. 
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

    persistent DATA
    if isempty(DATA)
        % Need setup
        if exist('OCTAVE_VERSION', 'builtin') > 0
            out = octave_setup();
        else
            out = matlab_setup();
        end
        DATA = out;
    end
    if nargin == 0
        out = DATA;
    else
        out = DATA.(arg);
    end
end

function s = octave_setup()
    v = ver();
    isMain = arrayfun(@(x) strcmpi(x.Name, 'octave'), v);
    hasGUI = strcmpi(graphics_toolkit(), 'qt'); % Mostly where GUIs are supported
    s = get_config(v(isMain), ...
                    'gui', hasGUI, ...
                    'desktop', isguirunning(), ...
                    'richtext', false);
end

function s = matlab_setup()
    v = ver();
    isMain = arrayfun(@(x) strcmpi(x.Name, 'matlab'), v);
    isDesktop = usejava('desktop');
    s = get_config(v(isMain), ...
                   'desktop', isDesktop, ...
                   'richtext', isDesktop, ...
                   'gui', isDesktop);
end

function s = get_config(ver_struct, varargin)
    if ispc()
        os = 'windows';
    elseif ismac()
        os = 'macos';
    else
        % Must be Linux or Unix since we already
        % excluded mac os.
        os = 'linux/unix';
    end
    ver_struct = ver_struct(1); % Just in case.
    name = ver_struct.Name;
    ver_char = ver_struct.Version;
    if exist('strsplit', 'file')
        tmp = strsplit(ver_char, '.');
        v1 = tmp{1};
        v2 = tmp{2};
    else
        dots = find(ver_char == '.');
        if numel(dots) > 1
            ver_char = ver_char(1:dots(2)-1);
        end
        v1 = ver_char(1:(dots(1)-1));
        v2 = ver_char(dots(1)+1:end);
    end
    major = str2double(v1);
    minor = str2double(v2);

    s = struct('platform',      name, ...
               'version',       ver_char, ...
               'major',         major, ...
               'minor',         minor, ...
               'octave',        strcmpi(name, 'octave'), ...
               'matlab',        strcmpi(name, 'matlab'), ...
               'os',            os, ...
               'desktop',       false, ...
               'gui',           false, ...
               'jvm',           usejava('jvm'), ...
               'richtext',      false);
   s = merge_options(s, varargin{:});
end
