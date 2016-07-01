function [updata, report] = upAbsPerm(block, updata, varargin)
opt = struct(...
    'method',     'pressure', ...
    'dims',       1:3, ...
    'psolver',    'tpfa', ...
    'dp',         1*barsa ...
    );
opt = merge_options(opt, varargin{:});

if nargin==1
    updata = [];
end

if strcmpi(opt.method, 'pressure')
    f = @() upAbsPermPres(block, updata, 'dims', opt.dims, ...
        'psolver', opt.psolver, 'dp', opt.dp);
else
    f = @() upAbsPermAvg(block, updata, 'dims', opt.dims, ...
        'method', opt.method);
end

if nargout > 1
    [updata, report] = f();
else
    updata = f();
    report = [];
end

end