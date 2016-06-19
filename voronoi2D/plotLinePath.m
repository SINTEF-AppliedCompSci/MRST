function [fig] = plotLinePath(lines, varargin)
% Plot the line paths stored in a cell

% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.

hold on
for i = 1:numel(lines)
  fig = plot(lines{i}(:, 1), lines{i}(:, 2) ,varargin{:});
end

end