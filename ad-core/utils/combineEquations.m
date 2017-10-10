function h = combineEquations(varargin)
% Combine equations. For doubles, this is equivialent to a vertical
% concatenation. Please note that the full implementation used is found as
% a part of the ADI class. This is the fallback for doubles.
    h = vertcat(varargin{:});
end
