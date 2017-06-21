function fn = getUniformInterpRegLinear(xs, fs, reg)
% Intermediate interpolator. Support for multiple regions, with the caveaet
% that all functions must be given on the same set of points.
    ns = numel(xs);
    slope = bsxfun(@rdivide, diff(fs, 1, 2), diff(xs, 1, 2));
    slope = slope(:, [1:(ns-1), ns-1]);
    
    % Dynamic part
    fs_unpack = reshape(fs', 1, [])';
    slope_unpack = reshape(slope', 1, [])';
    fn = @(x) interpReg(x, xs, fs_unpack, slope_unpack, reg);
end

function f = interpReg(x, xs, fs, slope, reg)
    xv = double(x);
    act = bsxfun(@gt, xv, [-inf, xs(1:end-1)]) & bsxfun(@le, xv, xs);
    [ii, jj] = find(act);
    jj = jj - 1;

    nx = numel(xv);
    ns = numel(xs);
    nreg = size(fs, 1)/ns;
    % Unpack regions, since the grid is assumed uniform
    pick = sparse((1:nx)' , jj + (reg-1)*ns, 1, nx, ns*nreg);
    f0 = pick*fs;
    dfdx = pick*slope;

    f =  f0 + dfdx.*(x - xs(jj)');
end