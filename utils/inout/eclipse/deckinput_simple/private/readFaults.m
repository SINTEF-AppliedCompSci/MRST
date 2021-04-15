function grdecl = readFaults(fid,grdecl)
% Read FAULTS keyword
%
%  SYNOPSIS
%
%   grdecl = readFaults(fid,grdecl)
%
%  PARAMETERS:
%   grdecl - a grdecl structure
%
%  RETURN:
%   grdecl - a grdecl structure

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

kw='FAULTS';
if(~isfield(grdecl,kw))
    grdecl.(kw)=struct;
end
lin = fgetl(fid);
count = 1;
while isempty(regexp(lin, '^/', 'once'))
    % Skip blank lines and comments.
   if(~(isempty(lin) || all(isspace(lin)) || ~isempty(regexp(lin, '^--', 'match'))))
    split = regexp(lin, '(\w\.*\-*\**)+', 'match');
    name=char(split(1));
    if(~length(split)==8)
        error(['Wrong line in readOpertor ', kw ,' for ',name])
    end
    dir=char(split(8));
    cells=zeros(1,6);
    for i=2:7
       cells(i-1)=str2num(char(split(i)));
    end
    if(~isfield(grdecl.(kw),name))
        grdecl.(kw).(name)=struct('cells',[],'dir',[]);
    end
    grdecl.(kw).(name).dir=[grdecl.(kw).(name).dir;dir];
    grdecl.(kw).(name).cells=[grdecl.(kw).(name).cells;cells];
   end
   lin = fgetl(fid);
   count =count +1;
end
