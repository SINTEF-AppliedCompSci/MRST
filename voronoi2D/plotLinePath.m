function [fig] = plotLinePath(lines, varargin)
% plotLinePath plots the lines in a cell array
%
% SYNOPSIS:
%       plotGrid(lines)
%       plotGrid(lines, 'pn1', pv1, ...)
%
% PARAMETERS:
%   lines     - cell array of lines. Each element in the cell array is a 
%               n X 2 array of coordinates describing a line
%   'pn'/pv   - List of preperty names/property values. OPTIONAL.
%               This list will be passed directly on to function PLOT.
%
% Returns
%   fig       - Function handle to axis.
%
% EXAMPLE:
%   l = {[0,0;0.5,0;1,1],[0,1;1,0]};
%   plotLinePath(l)
% SEE ALSO
%   plotGrid

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
hold on
for i = 1:numel(lines)
  fig = plot(lines{i}(:, 1), lines{i}(:, 2) ,varargin{:});
end

end