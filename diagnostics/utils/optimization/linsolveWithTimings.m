function out = linsolveWithTimings(A, x, linsolve)
persistent tocs
if isempty(tocs)
    tocs = [];
end
if nargin == 0
    out = tocs;
    tocs = [];
else
    if nargin < 3
        linsolve = @mldivide;
    end
    tic;
    out = linsolve(A,x);
    t   = toc;
    tocs = [tocs, t];
end
end

