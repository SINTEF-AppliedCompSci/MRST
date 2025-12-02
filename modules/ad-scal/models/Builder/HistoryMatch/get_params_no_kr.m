function params_no_kr = get_params_no_kr(model)  
%
% DESCRIPTION: get the number of the parameters for the relative
%              permeability saturation function, used during the history 
%              matching
%
% SYNOPSIS:
%   params_no_kr = get_params_no_kr(model)
%
% PARAMETERS:
%   - model: struct on which the history matching process is being run
%
% RETURNS:
%   params_no_kr - number of the parameters for the relative permeability
%   saturation function
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
if model.history_match.kr.status 
    if strcmp(model.history_match.kr.type,'MODIFIED-COREY')
        params_no_kr = nargin('modified_corey_kr');
    elseif strcmp(model.history_match.kr.type,'BROOKS-COREY')
        params_no_kr = nargin('brooks_corey_kr');
    elseif strcmp(model.history_match.kr.type,'MODIFIED-COREY-MASALMEH')
        params_no_kr = nargin('modified_corey_masalmeh_kr');
    elseif strcmp(model.history_match.kr.type,'LET')
        params_no_kr = nargin('LET_kr');
    elseif strcmp(model.history_match.kr.type,'BURDINE')
        params_no_kr = nargin('burdine_kr') - 2;
    elseif strcmp(model.history_match.kr.type,'VAN-GENUCHTEN-BURDINE')
        params_no_kr = nargin('van_genuchten_burdine_kr') - 2;
    else
        params_no_kr = [];
    end
else
    params_no_kr = [];
end