function grdecl = readMultfelt(fid,grdecl)
% Read MULTFELT keyword
%
%  SYNOPSIS
%
%   grdecl = readMultFelt(fid,grdecl)
%
%  PARAMETERS:
%   grdecl - a grdecl structure which have to have defined the
%            faults which multfelt have values for
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
    if(~length(split)==2)
        error(['Wrong line in readMULTFELT ', kw ,' for ',name])
    end
    multflt =str2num(split{2});
    if(~isfield(grdecl.(kw),name))
        error('Fault multiplier for not named fault')
    else
        grdecl.(kw).(name).multflt=multflt;
    end
   end
   lin = fgetl(fid);
   count =count +1;
end
