function grdecl = readOperator(fid,grdecl,kw)
% Read MULTFELT keyword
%
%  SYNOPSIS
%
%   grdecl = readFaults(fid,grdecl)
%
%  PARAMETERS:
%   grdecl - a grdecl structure which have to have defined the
%            faults which multfelt have values for
%
%  RETURN:
%   grdecl - a grdecl structure
%
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
    value=str2num(char(split(2)));
    region=zeros(1,6);
    for i=3:8
       region(i-2)=str2num(char(split(i)));
    end
    if(~isfield(grdecl.(kw),name))
        grdecl.(kw).(name)=struct('value',[],'region',[]);
    end
    grdecl.(kw).(name).value=[grdecl.(kw).(name).value;value];
    grdecl.(kw).(name).region=[grdecl.(kw).(name).region;region];
   end
   lin = fgetl(fid);
   count =count +1;
end
