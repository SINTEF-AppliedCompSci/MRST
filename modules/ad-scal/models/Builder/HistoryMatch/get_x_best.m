function [x_best, fval_sum_best] = get_x_best(model)
%
% DESCRIPTION: chooses the best combination with the lowest error sum from
%              the output of multi objective optimization
%
% SYNOPSIS:
%   [x_best, fval_sum_best] = get_x_best(model)
%
% PARAMETERS:
%   model - struct which includes the result of the multi objective history
%   matching in the history_match field
%
% RETURNS:
%   x_best - best iteration parameters
%   fval_sum_best - minimum error
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
x_list = model.history_match.x_list;
fval_sum = sum(model.history_match.fval_list,2);
sorted_fval_sum = sort(fval_sum);
fval_sum_best = fval_sum(fval_sum==sorted_fval_sum(1),:);
x_best = x_list(fval_sum==sorted_fval_sum(1),:);