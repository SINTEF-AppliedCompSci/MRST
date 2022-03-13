function [s, f, fv, fg] = buckleyLeverettProfile(varargin)
    % Get analytical Buckely-Leverett solution for a simple Riemann problem

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('vT'  , 1  , ...
                 'g'   , 0  , ...
                 'nW'  , 2  , ...
                 'nN'  , 2  , ...
                 'muW' , 1  , ...
                 'muN' , 1  , ...
                 'rhoW', 0.5, ...
                 'rhoN', 1  , ...
                 'sL'  , 1  , ...
                 'sR'  , 0  );
             
    opt = merge_options(opt, varargin{:});
    
    vT   = opt.vT;
    g    = opt.g;
    nW   = opt.nW;
    nN   = opt.nN;
    muW  = opt.muW;
    muN  = opt.muN;
    rhoW = opt.rhoW;
    rhoN = opt.rhoN;
    sL   = opt.sL;
    sR   = opt.sR;
    
    a  = @(s) s.^nW/muW;
    da = @(s) nW.*s.^(nW-1)/muW;
    b  = @(s) (1-s).^nN/muN;
    db = @(s) -nN.*(1-s).^(nN-1)/muN;
    
    fv = @(s) a(s)./(a(s) + b(s)).*vT;
    fg = @(s) -fv(s).*b(s).*(rhoW - rhoN)*g;
    f  = @(s) fv(s) + fg(s);
    
    dfv = @(s) (da(s).*(a(s) + b(s)) - a(s).*(da(s) + db(s)))./(a(s) + b(s)).^2.*vT;
    dfg = @(s) dfv(s).*fg(s) + fv(s).*nN.*(1-s).^(nN-1)/muN.*(rhoW - rhoN)*g;
    df  = @(s) dfv(s) + dfg(s);
    
    res = @(s) (f(s)./s - df(s)).*(s > 0) + 1*(s <= 0);

    sMin = 0.005;
    s0 = sMin;
    while s0 == sMin || s0 == 1-sMin
        sMin = sMin*2;
        s = linspace(sMin,1-sMin,10000);
        r = res(s);
        [~, ix] = min(abs(r));
        s0 = s(ix);
    end
    
    s = linspace(0,1,1000)';
    
    slope = (f(s0) - f(0))/s0;
    
    v = df(s).*(s > s0) + (linspace(0,1,numel(s))' + max(df(s))).*(s <= s0);
    dfinv = @(x,t) interpTable(v, s, x./t);
    
    if all(diff(v) == 0)
        dfinv = @(x,t) (sL + sR)/2 + x;
    end
    
    s = @(x,t) (x./t <= df(sL)).*sL + ...
               dfinv(x,t).*(x./t > df(sL) & x./t < slope) + ...
               (x./t >= slope)*sR;

end
