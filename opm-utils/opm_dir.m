function dir = opm_dir(varargin)

% Retrieve full path of opm-core root directory.
%
% SYNOPSIS:
%   dir = opm_dir()
%   dir = opm_dir(user)
%
% PARAMETERS:
%   user - if supplied, used to choose a directory.
%
% RETURNS:
%   root - Full path to opm-core root directory.

%{
#COPYRIGHT#
%}

switch nargin
    case 0
        user = 'default';
    case 1
        user = varargin{1};
    otherwise
        error('Too many arguments');
end

switch user
    case 'default'
        dir = 'if_it_is_possible_to_choose_a_useful_default_please_do';
    case 'atgeirr'
        dir = '/Users/atgeirr/opm/opm-core-debug';
    case 'hnil'
        dir = '/home/hnil/heim/SVN/OPM_related_hg/builds/db/opm-core/';
    case 'bska'
        dir = 'unknown_opm_dir_for_bska';
    otherwise
        error(['Excplicit user ', user, ' requested, but user is unknown'])
end
