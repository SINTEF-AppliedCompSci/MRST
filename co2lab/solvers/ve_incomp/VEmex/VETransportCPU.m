function varargout = VETransportCPU(varargin)  
% build/link to mex file

d = fileparts(mfilename('fullpath'));
mexfile = ['VETransportCPU.', mexext];
mexfullfile = fullfile(d, 'build', 'mex', mexfile);
if exist(mexfullfile, 'file')
   tofile  = fullfile(d, mexfile);
   system(['ln -s ', mexfullfile, ' ', tofile]);
   [varargout{1:nargout}] = VETransportCPU(varargin{:});
else
   error(['You need to build VEmex manually. \n',...
         'For instructions see:' , fullfile(d, 'README.')]);
end

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