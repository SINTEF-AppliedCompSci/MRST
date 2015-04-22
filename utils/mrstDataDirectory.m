function varargout = mrstDataDirectory(varargin)
%Set or retrieve the current canonical data directory for MRST
%
% SYNOPSIS:
%   % Query directory
%   mrstDataDirectory();
%   datadir = mrstDataDirectory();
%
%   % Set directory
%   mrstDataDirectory('/some/path');
%
%
% OPTIONAL PARAMETERS:
%
%   datadir - If provided, the current data directory for the MRST session
%             will be set to this directory. If the directory itself does
%             not exist, a warning will be given and the directory will not
%             be changed.
%
% RETURNS:
%   datadir - Current datadir. If not requested, this will instead be
%             printed to the command window.
%
% NOTE:
%   If you do not wish to use the default MRST data directory, consider
%   placing a call to mrstDataDirectory in your startup_user.m file located
%   under the mrst root directory (see ROOTDIR())
%
% SEE ALSO:
%   mrstModule, mrstDatasetGUI, mrstPath

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
    persistent DATADIR
    if isempty(DATADIR)
        % Default datadir
        DATADIR = fullfile(ROOTDIR, 'examples', 'data');
        if exist(DATADIR, 'dir') ~= 7
            mkdir(DATADIR)
        end
    end
    
    if nargin == 0 && nargout == 0
        fprintf('The current MRST data directory is set to:\n\t%s\n', ...
                 DATADIR);
    else
        newDir = varargin{1};
        assert(ischar(newDir), 'Data directory must be a string');
        if isdir(newDir)
            DATADIR = newDir;
        else
            warning(['Supplied directory ''', newDir, ''' is not a directory.', ...
                     ' Data directory has not been changed']);
        end
    end
    
    if nargout > 0
        varargout{1} = DATADIR;
    end
end
