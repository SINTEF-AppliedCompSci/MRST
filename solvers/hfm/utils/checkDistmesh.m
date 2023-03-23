function varargout = checkDistmesh()
% Checks if distmesh_2d file is installed in the MATLAB path
%
% SYNOPSIS:
%   hasLSI = checkLineSegmentIntersect();
%   checkLineSegmentIntersect();
%
% DESCRIPTION:
%   This function checks if the file 'distmesh_2d.m'
%   (http://people.sc.fsu.edu/~jburkardt/m_src/distmesh/distmesh.html) is
%   installed in the current MATLAB path and advises the user to install it
%   if it is not found.
%
%   If called with a return argument, it will produce a boolean indicating
%   the status of the distmesh_2d.m file. The application checking
%   for distmesh_2d.m can then fall back to other implementations
%   of graph algorithms as they please. If called without any return
%   arguments, it will throw an error if distmesh_2d.m was not
%   found, along with a link and install instructions for the user to
%   rectify the missing installation.
%
%
% RETURNS:
%   foundLSI (OPTIONAL) - If distmesh_2d.m was found. Note that
%                         the nature of the utility changes if called with
%                         return argument.
% 

%{
Copyright 2009-2015: SINTEF ICT, Applied Mathematics.

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

    foundLSI = exist('distmesh_2d.m', 'file') > 0 ;

    if nargout == 0 && ~foundLSI;
        % If called without return argument, act as a block and throw error
        error('mrst:Missingdistmesh_2d', ...
            ['You do not have the distmesh_2d.m file installed!\n\n'...
            'Please download it from <a href="http://people.sc.fsu.edu/~jburkardt/m_src/distmesh/distmesh.html">Distmesh 2D</a>'...
            ' and install it in your MATLAB path.\n'...
            'If you already have it installed somewhere but not on your path, use\n\n'...
            '\tmrstModule(''add'', <path-to-distmesh_2d.m>)' ...
            '\n\nto make it available for MRST.'])
    end
    varargout{1} = foundLSI;
end
