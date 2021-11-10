function [path_string, present] = getCanonicalPath(path_string)
%Get absolute (canonical) path of file or directory in MATLAB and Octave
%
% SYNOPSIS:
%   path_string            = getCanonicalPath(path_string)
%   [path_string, present] = getCanonicalPath(path_string)
%
% REQUIRED PARAMETERS:
%   path_string - String of potentially non-existent path that may contain
%                 relative indicators in format.
%
% RETURNS:
%   path_string - Best attempt at getting the canonical path.
%
%   present     - Whether or not file/directory is actually present on the
%                 local computer system.

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
        isfolder_safe = exist('isfolder', 'builtin') || ...
                        exist('isfolder', 'file');
    end

    [is_directory, canonicalise_path] = ...
       get_path_operations(isOctave, isfolder_safe);

    % FULLFILE(path_string, '.') transparently handles presence or absence
    % of terminating FILESEP on 'path_string'.
    %
    pth     = fileparts(fullfile(path_string, '.'));
    present = is_directory(pth);
    
    path_string = canonicalise_path(pth, present);
end

%--------------------------------------------------------------------------

function [is_directory, canonicalise] = ...
      get_path_operations(isOctave, isfolder_safe)

   is_directory = @isfolder;
   if ~isfolder_safe
      is_directory = @isdir;
   end

   canonicalise = @canonicalise_path_matlab;
   if isOctave
      canonicalise = @canonicalise_path_octave;
   end
end

%--------------------------------------------------------------------------

function path_string = canonicalise_path_matlab(pth, present)
   if present
      % Use the side effect that function WHAT, in the case of an input
      % argument that names an existing directory, also resolves symbolic
      % links and special directories DOT and DOTDOT (i.e., '.' and '..').
      w           = what(pth);
      path_string = w.path;
   else
      path_string = canonicalise_path_matlab_jvm(pth);
   end
end

%--------------------------------------------------------------------------

function path_string = canonicalise_path_octave(pth, present)
   canonicalise = @canonicalize_file_name;
   if ~present
      canonicalise = @make_absolute_filename;
   end

   path_string = canonicalise(pth);
end

%--------------------------------------------------------------------------

function path_string = canonicalise_path_matlab_jvm(pth)
   try
      % Input directory (pth) does not (yet?) exist.  Try Java's
      % getCanonicalPath() to resolve any intermediate pathname components.
      path_string = java.io.File(pth).getCanonicalPath();
   catch
      % Java not available (-nojvm?).  Punt back to caller.
      path_string = pth;
   end
end
