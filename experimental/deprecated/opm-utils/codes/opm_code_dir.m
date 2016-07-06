function root = opm_code_dir
%Retrieve full path of Toolbox installation directory.
%
% SYNOPSIS:
%   root = ROOTDIR
%
% PARAMETERS:
%   none.
%
% RETURNS:
%   root - Full path to MRST installation directory.

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

% $Date: 2010-03-16 11:01:57 +0100 (Tue, 16 Mar 2010) $
% $Revision: 4336 $

nm = mfilename('fullpath');
%nm /home/hnil/heim/SVN/OPM_related_hg
ix = strfind(nm, filesep);
if ~isempty(ix),
   root = nm(1 : ix(end)-1);
else
   root = nm;
end
root=fullfile(root,'codes');
root=fullfile('/home/hnil/heim/SVN/OPM_related_hg');

