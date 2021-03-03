function [path_string, present] = getCanonicalPath(path_string)
%Get absolute (canonical) path of file in MATLAB and Octave
%
% SYNOPSIS:
%   canonical_name = getCanonicalPath('/some/path')
%
% REQUIRED PARAMETERS:
%   path_string - String of potentially non-existent path that may contain
%                 relative indicators in format..
%
% RETURNS:
%   path_string - Best attempt at getting the canonical path.
%
%   present     - Indicator if file is actually present.

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
    persistent isOctave isfolder_safe
    if isempty(isOctave)
        isOctave = mrstPlatform('octave');
    end
    if isempty(isfolder_safe)
        isfolder_safe = exist('isfolder', 'builtin') || exist('isfolder', 'file');
    end
    % The FULLFILE(dname, '.') construct is to transparently handle presence
    % or absence of a terminating FILESEP on 'dname'.
    %
    pth = fileparts(fullfile(path_string, '.'));
    if isfolder_safe
       present = isfolder(pth); %#ok
    else
       present = isdir(pth); %#ok
    end

    if isOctave
        % OCTAVE
        % Use builtin routines
        if present
            path_string = canonicalize_file_name(pth);
        else
            path_string = make_absolute_filename(pth);
        end
    else
        % MATLAB
        if present
            % Use the side effect that function WHAT, in the case of an input
            % argument that names an existing directory, also resolves symbolic
            % links and special directories DOT and DOTDOT (i.e., '.' and '..').
            w     = what(pth);
            path_string = w.path;
        else
            try
                % Input directory (dname) does not (yet?) exist.  Try Java's
                % getCanonicalPath() to resolve any intermediate pathname
                % components.
                path_string = java.io.File(pth).getCanonicalPath();
            catch
                % Java not available (-nojvm?).  Punt back to caller.
                path_string = pth;
            end
        end
    end
end
