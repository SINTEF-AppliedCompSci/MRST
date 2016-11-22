function tf = datasetHasCustomDownloadFcn(info)
%Predicate for whether or not a dataset provides a custom download function
%
% SYNOPSIS:
%   tf = datasetHasCustomDownloadFcn(info)
%
% PARAMETERS:
%   info - A dataset information structure as defined by function
%          datasetInfoStruct.
%
% RETURNS:
%   tf - Whether or not the dataset identified by 'info' does provide a
%        custom download function.

%{
#COPYRIGHT#
%}

   tf = (numel(info.downloadFcn) == 1) && ...
         isa(info.downloadFcn{1}, 'function_handle') && ...
        (nargin(info.downloadFcn{1}) == 0);
end
