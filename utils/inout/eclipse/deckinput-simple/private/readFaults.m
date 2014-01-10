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
%
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
