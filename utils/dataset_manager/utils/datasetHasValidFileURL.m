function tf = datasetHasValidFileURL(info)
%Predicate for whether or not a dataset provides a valid file URL
%
% SYNOPSIS:
%   tf = datasetHasValidFileURL(info)
%
% PARAMETERS:
%   info - A dataset information structure as defined by function
%          datasetInfoStruct.
%
% RETURNS:
%   tf - Whether or not the dataset identified by 'info' does provide a
%        valid file URL from which to download the dataset.

%{
#COPYRIGHT#
%}

   tf = (~ isempty(info.fileurl)) && ischar(info.fileurl);
end
