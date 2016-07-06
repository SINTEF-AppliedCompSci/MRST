function varargout = checkBGL()
%Checks if matlab boost graph library is installed.
%
% SYNOPSIS:
%   hasBGL = checkBGL();
%   checkBGL();
%
% DESCRIPTION:
%   This function checks if MatlabBGL
%   (http://www.mathworks.com/matlabcentral/fileexchange/10922) is
%   installed and advises the user to install it if it is not found.
%
%   If called with a return argument, it will produce a boolean indicating
%   the status of the BGL library. The application checking for BGL can
%   then fall back to other implementations of graph algorithms as they
%   please. If called without any return arguments, it will throw an error
%   if BGL was not found, along with a link and install instructions for
%   the user to rectify the missing installation.
%
%
% RETURNS:
%   foundBGL (OPTIONAL) - If BGL was found. Note that the nature of the
%                         utility changes if called with return argument.
% 

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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

    foundBGL = exist('matlab_bgl', 'dir') > 0 ||...
               exist('matlab-bgl', 'dir') > 0;

    if nargout == 0 && ~foundBGL;
        % If called without return argument, act as a block and throw error
        error('mrst:MissingBGL', ['You do not have the MATLAB BGL library installed!\n\n'...
                              'Please download it from <a href="http://www.mathworks.com/matlabcentral/fileexchange/10922">File Exchange</a>'...
                              ' and install it in your MATLAB path.\n'...
                              'If you already have it installed somewhere not on your path, use\n\n'...
                              '\tmrstModule(''add'', <path-to-matlab_bgl>)' ...
                              '\n\nto make it available for MRST.'])
    end
    varargout{1} = foundBGL;
end
