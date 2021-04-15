function [f, dfdx, dfdy] = interpPVT(T, x, v, flag, method, useMex)
% Interpolate function on the form f(x, v, flag)
% where flag indicates if the saturated curve is to be used for each entry
    % Interpolate PVT-type curves
    if nargin < 6
        useMex = false;
        if nargin < 5
            method = 'box';
        end
    end
    xAD = isa(x, 'ADI');
    vAD = isa(v, 'ADI');
    isAD = xAD || vAD;

    compDer = nargout > 1 || isAD;
    if numelValue(flag) == 1
        flag = repmat(flag, numelValue(x), 1);
    end
    [f, dfdx, dfdy] = interpolate(T, value(x), value(v), flag, method, compDer, useMex);
    % Deal with AD status
    if isAD
        if xAD
            % x is AD, v may be AD
            fa = x;
            fa.val = f;
            if vAD
                fa.jac = fa.timesJac(dfdx, dfdy, v.jac, x.jac); % note order of input
            else
                fa.jac = fa.lMultDiag(dfdx, x.jac);
            end
        else
            % v is AD, x is not
            fa = v;
            fa.val = f;
            fa.jac = fa.lMultDiag(dfdy, v.jac);
        end
        % We have made a AD copy with the correct value. Return that
        % instead.
        f = fa;
    end
end

function [yi, dyidxi, dyidvi] = interpolate(T, xi, vi, sat, method, compDer, useMex)
    if isa(T, 'function_handle')
        if compDer
            [yi, dyidxi, dyidvi] = T(xi, vi, sat);
        else
            yi = T(xi, vi, sat);
            dyidxi = [];
            dyidvi = [];
        end
    else
        yi = zeros(size(xi));
        dyidxi = zeros(size(xi));
        dyidvi = zeros(size(xi));
        if isempty(xi)
            [yi, dyidxi, dyidvi] = deal([]);
            return;
        else
            % Shouldn't extrapolate downwards on saturated curve - assume
            % that this is somehow undersaturated.
            sat(value(xi) < T.sat.x(1)) = false;
            usat = ~sat;
            [yi(sat), dyidxi(sat)] = interp1_internal(T.sat.x, T.sat.F, xi(sat), compDer, useMex);
            if any(usat)
                if compDer
                    [yi(usat), dyidxi(usat), dyidvi(usat)] = interp2DPVT_local(T, xi(usat), vi(usat), method, useMex);
                else
                    yi(usat) = interp2DPVT_local(T, xi(usat), vi(usat), method, useMex);
                end
            end
        end
    end
end

function varargout = interp2DPVT_local(T, xi, vi, method, useMex)
    varargout = cell(nargout, 1);
    switch method
        case 'box'
            [varargout{:}] = interp2DPVT_box(T, xi, vi, useMex);
        case 'parallel'
            [varargout{:}] = interp2DPVT_parallel(T, xi, vi, useMex);
        case 'linshift'
            [varargout{:}] = interp2DPVT_linshift(T, xi, vi, useMex);
        otherwise
            error('Invalid interpolation strategy');
    end
end

function [yi, dyidxi, dyidvi] = interp2DPVT_box(T, p, r, useMex)
    % Interpolate in a "box" - old default MRST interpolation
    compDer = (nargout>1);
    yil = zeros(size(p));
    yir = zeros(size(p));
    dyidxil = zeros(size(p));
    dyidxir = zeros(size(p));

    v = T.key; 
    [bin, w, dwdvi] = getBins(r, v, compDer);

    for tn = 1:numel(v)
        ixl = (bin==(tn-1));
        ixr = (bin==tn);    
        if ~(any(ixl) || any(ixr))
            continue
        end
        tab = T.expanded{tn};
        P = tab.x;
        F = tab.F;

        [yil(ixl), dyidxil(ixl)] = interp1_internal(P, F, p(ixl), compDer, useMex);
        [yir(ixr), dyidxir(ixr)] = interp1_internal(P, F, p(ixr), compDer, useMex);
    end

    yi = yil.*w + yir.*(1-w);
    if compDer
        dyidxi = dyidxil.*w + dyidxir.*(1-w);
        dyidvi = (yil-yir).*dwdvi;
    end
end

function [f, dfdp, dfdr] = interp2DPVT_parallel(T, p, r, useMex)
    % Interpolate in parallelogram. We first interpolate along the
    % saturated curve to estimate the bubble point pressure. Then, we
    % interpolate the pair of closest undersaturated 1D tables by the
    % relative distance to the bubble point.
    compDer = (nargout>1);
    f_l = zeros(size(p));
    f_r = zeros(size(p));
    if compDer
        dfdDp_l = zeros(size(p));
        dfdDp_r = zeros(size(p));
    end
    v = T.key;
    [bin, w, dwdr] = getBins(r, v, compDer);
    [ps, dpsdr] = interp1_internal(T.key, T.sat.x, r, compDer, useMex);
    Dp = p - ps;
    for tn = 1:numel(v)
        % Loop over each curve and compute the interpolated value and slope
        % along all points that have that line as either the left or the
        % right boundary.
        ixl = (bin==(tn-1));
        ixr = (bin==tn);
        active = ixl | ixr;
        if ~any(active)
            continue
        end
        left = bin(active) == tn-1;
        % Get the n-th undersaturated table, subtract the saturated
        % pressure and interpolate
        tab = T.expanded{tn};
        DP = tab.x;
        DP = DP - DP(1);
        F = tab.F;

        Dpa = Dp(active);
        [fi, df] = interp1_internal(DP, F, Dpa, compDer, useMex);
        f_l(ixl) = fi(left);
        f_r(ixr) = fi(~left);
        if compDer
            dfdDp_l(ixl) = df(left);
            dfdDp_r(ixr) = df(~left);
        end
    end

    f = f_l.*w + f_r.*(1-w);
    if compDer
        % Finalize the derivatives
        dfdp = dfdDp_l.*w + dfdDp_r.*(1-w);
        dfdr = (f_l-f_r).*dwdr;
        % Account for shifted coordinates for p table (chain rule)
        dfdr = dfdr - dfdp.*dpsdr;
    end
end


function [yi, dyidxi, dyidvi] = interp2DPVT_linshift(T, xi, vi, useMex)
    % Interpolate in box, then linear shift
    compDer = (nargout>1);

    yil = zeros(size(xi));
    yir = zeros(size(xi));
    dyidxil = zeros(size(xi));
    dyidxir = zeros(size(xi));
    v = T.key;
    [bin, w, dwdvi] = getBins(vi, v, compDer);
    x = T.sat.x;
    dx = (x(bin+1)-x(bin));
    a = (x(bin+1)-x(bin))./(v(bin+1)-v(bin));


    for tn = 1:numel(v)
        tab = T.expanded{tn};

        ixl = (bin==(tn-1));
        ixr = (bin==tn);
        % shift values such that intp is linear along (x,v)
        xils = xi(ixl) + (1-w(ixl)).*dx(ixl);
        xirs = xi(ixr) - w(ixr).*dx(ixr);
        
        P = tab.x;
        F = tab.F;
        [yil(ixl), dyidxil(ixl)] = interp1_internal(P, F, xils, compDer, useMex);
        [yir(ixr), dyidxir(ixr)] = interp1_internal(P, F, xirs, compDer, useMex);
    end
    yi = yil.*w + yir.*(1-w);
    if compDer
        dyidxi = dyidxil.*w + dyidxir.*(1-w);
        dyidvi = (yil-yir).*dwdvi - dyidxi.*a;
    end
end

function [f, df] = interp1_internal(X, F, x, compDer, useMex)
    if useMex
        if compDer
            [f, df] = interpTableMEX(X, F, x);
        else
            f = interpTableMEX(X, F, x);
            df = zeros(numel(x), 1);
        end
    else
        f  = fastInterpTable(X, F, x);
        if compDer
            df = dinterpq1(X, F, x);
        else
            df = zeros(numel(x), 1);
        end
    end
end

function [bin, w, dwdr] = getBins(r, v, compDer)
    lims = v; lims(1)=-inf; lims(end)=inf;
    [bin,bin] = histc(r, lims);                                     %#ok
    % Distance between saturated points
    w = (r - v(bin))./(v(bin+1) - v(bin));
    if compDer
        dwdr = 1./(v(bin+1)-v(bin));
    else
        dwdr = [];
    end
end
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

