function dir = opm_dir(varargin)

% Retrieve full path of opm-core root directory.
%
% SYNOPSIS:
%   dir = opm_dir()
%   dir = opm_dir(user)
%
% PARAMETERS:
%   user - if supplied, used to choose a directory.
%
% RETURNS:
%   root - Full path to opm-core root directory.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

switch nargin
    case 0
        user = 'default';
    case 1
        user = varargin{1};
    otherwise
        error('Too many arguments');
end

switch user
    case 'default'
        dir = 'if_it_is_possible_to_choose_a_useful_default_please_do';
    case 'atgeirr'
        dir = '/Users/atgeirr/opm/opm-core-debug';
    case 'hnil'
        dir = '/home/hnil/heim/SVN/OPM_related_hg/builds/db/opm-core/';
    case 'bska'
        dir = 'unknown_opm_dir_for_bska';
    otherwise
        error(['Excplicit user ', user, ' requested, but user is unknown'])
end
