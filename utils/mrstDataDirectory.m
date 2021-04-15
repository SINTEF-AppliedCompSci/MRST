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
%   `mrstModule`, `mrstDatasetGUI`, `mrstPath`, `getDatasetPath`

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

    if nargin == 0 && nargout == 0
        fprintf('The current MRST output directory is set to:\n\t%s\n', ...
                 mrstSettings('get', 'dataDirectory'));
    elseif nargin == 1
        mrstSettings('set', 'dataDirectory', varargin{1})
    elseif nargin > 1
        error('Function must be called with one or zero inputs.');
    end

    if nargout > 0
        varargout{1} = mrstSettings('get', 'dataDirectory');
    end
end
