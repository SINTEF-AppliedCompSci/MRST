function [Kup, report] = upAbsPerm(block, varargin)
opt = struct(...
    'method',     'pressure', ...
    'dims',       1:3, ...
    'psolver',    'tpfa', ...
    'dp',         1*barsa ...
    );
opt = merge_options(opt, varargin{:});

if strcmpi(opt.method, 'pressure')
    f = @() upAbsPermPres(block, 'dims', opt.dims, ...
        'psolver', opt.psolver, 'dp', opt.dp);
else
    f = @() upAbsPermAvg(block, 'dims', opt.dims, 'method', opt.method);
end

if nargout > 1
    [Kup, report] = f();
else
    Kup = f();
    report = [];
end

end