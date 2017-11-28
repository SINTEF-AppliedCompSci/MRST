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
