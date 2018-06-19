 function varargout = mrstOutputDirectory(varargin)
%Set or retrieve the current canonical data directory for MRST
%
% SYNOPSIS:
%   % Query directory
%   mrstOutputDirectory();
%   outdir = mrstOutputDirectory();
%
%   % Set directory
%   mrstOutputDirectory('/some/path');
%
%
% OPTIONAL PARAMETERS:
%
%   outdir  - If provided, the current output directory for the MRST session
%             will be set to this directory. If the directory itself does
%             not exist, a warning will be given and the directory will not
%             be changed.
%
% RETURNS:
%   outdir  - Current output dir. If not requested, this will instead be
%             printed to the command window.
%
% NOTE:
%   If you do not wish to use the default MRST output directory, consider
%   placing a call to mrstOutputDirectory in your startup_user.m file
%   located under the mrst root directory (see ROOTDIR())
%
% SEE ALSO:
%   `mrstModule`, `mrstOutputDirectory`, `mrstPath`

%{
Copyright 2009-2018 SINTEF Digital, Applied Mathematics.

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

    persistent OUTPUTDIR

    if isempty(OUTPUTDIR)
        OUTPUTDIR = initialize_output_directory(varargin{:});

        mlock
    end

    if nargin == 0 && nargout == 0
        fprintf('The current MRST output directory is set to:\n\t%s\n', ...
                 OUTPUTDIR);
    elseif nargin > 0
        newDir = varargin{1};
        assert(ischar(newDir), 'Output directory must be a string');

        if isdir(newDir)
            munlock
            OUTPUTDIR = newDir;
            mlock
        else
            warning(['Supplied directory ''', newDir, ''' is not a directory.', ...
                     ' Output directory has not been changed']);
        end
    end

    if nargout > 0
        varargout{1} = OUTPUTDIR;
    end
end

%--------------------------------------------------------------------------

function odir = initialize_output_directory(varargin)
    if (nargin > 0) && ischar(varargin{1})
        % Caller requested specific output directory.  Honour request.
        odir = varargin{1};
    else
        % No specific output directory requested.  Use default.
        odir = default_output_dir();
    end

    ensure_directory_exists(odir);
end

%--------------------------------------------------------------------------

function ensure_directory_exists(odir)
    if ~isdir(odir)
        [ok, msg, id] = mkdir(odir);

        if ~ok
            error(id, 'Failed to create Output Directory ''%s'': %s', ...
                  odir, msg);
        end
    end
end

%--------------------------------------------------------------------------

function odir = default_output_dir()
    odir = fullfile(ROOTDIR, 'output');
end
