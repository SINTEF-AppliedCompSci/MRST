function [x, f] = regularizeSaturationFunction(x, f, varargin)
%Undocumented Utility Function

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

    opt = struct('extrapolateEnd', true, ...
                 'monotoneMethod', 'remove', ...
                 'maxValue',    1, ...
                 'monotone',    true);
    opt = merge_options(opt, varargin{:});
    act = x >= 0 & isfinite(x) & isfinite(f);
    
    x = x(act);
    f = min(f(act), opt.maxValue);
    [x, ix] = unique(x);
    f = f(ix);

    if numel(x) <= 1
        return
    end
    
    if opt.monotone && any(diff(f) < 0)
        switch lower(opt.monotoneMethod)
            case 'remove'
                [f, x] = fixByRemoval(f, x);
            case 'integral'
                [f, x] = fixByIntegral(f, x);
            case 'slope'
                [f, x] = fixBySlope(f, x);
            case 'upper'
                [f, x] = fixByUpper(f, x);
            case 'lower'
                [f, x] = fixByLower(f, x);
            otherwise
                error('Unknown')
        end
    end

    if numel(x) <= 1
        return
    end
    
    xmin = min(x);
    xmax = max(x);
    
    if xmin > 0 
        slope = (f(2) - f(1))./(x(2) - x(1));
        x_added = x(1) - f(1)./slope;
        
        if x_added <= 0 || ~isfinite(x_added)
            x = [0; x];
            f = [0; f];
        else
            x = [0; x_added; x];
            f = [0; 0; f];
        end
    end
    
    if xmax < opt.maxValue
        if opt.extrapolateEnd
            slope = (f(end) - f(end-1))./(x(end) - x(end-1));
            x_added = x(end) + (1 - f(end))./slope;
            if x_added >= 1 || ~isfinite(x_added)
                
                f_added = f(end) + (1 - x(end))*slope;
                x_added = 1;
            else
                f_added = [1; 1];
                x_added = [x_added; 1];
            end
        else
            x_added = 1;
            f_added = f(end);
        end
        x = [x; x_added];
        f = [f; f_added];
    end
    [x, ix] = unique(x);
    f = f(ix);
end

function [f, x] = fixByUpper(f, x)
    for i = 2:numel(f)
        f(i) = max(f(i-1), f(i));
    end
end
function [f, x] = fixByLower(f, x)
    for i = (numel(f)-1):-1:1
        f(i) = min(f(i+1), f(i));
    end
end
function [f, x] = fixByIntegral(f, x)
    h = fixByUpper(f, x);
    % Get lower bound
    g = fixByUpper(f, x);
    % Set value as sum of upper and lower bound, with weights
    % chosen to ensure the same integral
    F = quadrature(f, x);
    G = quadrature(g, x);
    H = quadrature(h, x);

    a = (F-G)/(H-G);
    if ~isfinite(a)
        a = 1;
    end
    assert(a <= 1 && a >= 0);
    f = a.*h + (1-a).*g;
    assert(all(diff(f)>=0))
end

function [f, x] = fixBySlope(f, x)
    h = fixByUpper(f, x);
    % Get lower bound
    g = fixByLower(f, x);
    f = f(g == h);
    x = x(g == h);
end


function [f, x] = fixByRemoval(f, x)
    ok = [true; diff(f) >= 0];
    while any(~ok)
        f = f(ok);
        x = x(ok);
        ok = [true; diff(f) >= 0];
    end
end

function F = quadrature(f, x)
    f = f(1:end-1) + diff(f);
    F = sum(f.*diff(x));
end
