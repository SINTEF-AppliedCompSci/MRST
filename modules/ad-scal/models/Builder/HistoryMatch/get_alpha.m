function alpha = get_alpha(x, model)
%
% DESCRIPTION: builds the parameter input for the history matching of a
%              heterogeneous case, for the alpha parameter - this parameter
%              is between 0 an 1:
%              1 calculating the capillary scales all by pressure 
%              0 calculating the capillary scales all by saturation
%
% SYNOPSIS:
%   alpha = get_alpha(x, model)
%
% PARAMETERS:
%   - model: struct used for the history matching 
%   - x: the parameters used for the history match iterations
%
% RETURNS:
%   alpha: number between 0 and 1
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
x = x(:);
row_names = model.history_match.RowNames;
rock = model.experiment.rock;

% handle alpha
alpha_mask = strcmpi(row_names, 'alpha');
if any(alpha_mask)
    alpha = x(alpha_mask);
else
    if isfield(rock, 'alpha')
        alpha = rock.alpha.value;
    else
        alpha = 1;
    end
end
