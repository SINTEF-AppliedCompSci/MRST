function varargout = interp1q_mex(varargin)
%Undocumented Utility Function

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

% Get path of this file
p = fileparts(mfilename('fullpath'));

% Build compile string
mexcmd = sprintf('mex -O CFLAGS="\\$CFLAGS -std=c99" %s %s', ...
    sprintf('-outdir "%s"', p), ...
    sprintf('"%s"', fullfile(p,'interp1q_mex.c')) );

% Run compile string
eval(mexcmd);

% Call MEX'ed edition.
[varargout{1:nargout}] = interp1q_mex(varargin{:});

end
