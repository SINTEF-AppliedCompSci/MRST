function veroot = VEROOTDIR
%Retrieve full path of VE-root directory.
%
% SYNOPSIS:
%   veroot = VEROOTDIR
%
% PARAMETERS:
%   none.
%
% RETURNS:
%   veroot - Full path to VE module installation directory.

%{
#COPYRIGHT#
%}

nm = mfilename('fullpath');
ix = strfind(nm, filesep);
if ~isempty(ix),
   veroot = nm(1 : ix(end-2));
else
   veroot = nm;
end
