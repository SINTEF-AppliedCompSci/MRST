function opt = merge_options_relaxed(opt, varargin)
% A less general version of merge_options focused on specific choices:
%   - Arguments must match the names of fields exactly
%   - No type checking
%   - Errors for unsupported fields
%
% INTENTIONALLY UNDERDOCUMENTED, SUBJECT TO CHANGE.
    if nargin == 2 && iscell(varargin{1})
        % Can pass in varargin cell array directly
        varargin = varargin{1};
        n = numel(varargin);
    else
        % varargin is a normal list
        n = nargin - 1;
    end
    assert(mod(n, 2) == 0);
    for i = 1:2:(n-1)
        f = varargin{i};
        v = varargin{i+1};
        assert(isfield(opt, f));
        opt.(f) = v;
    end
end
