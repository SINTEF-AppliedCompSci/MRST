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
#COPYRIGHT#
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

