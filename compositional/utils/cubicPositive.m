function r = cubicPositive(a, b, c)
% Straightforward implementation of a cubic root solver for vectorized
% problems. See https://en.wikipedia.org/wiki/Cubic_function
    if nargin == 1
        % We recieved a matrix. Unpack into vectors.
        c = a(:, 3);
        b = a(:, 2);
        a = a(:, 1);
    end
    Q = (a.^2 - 3*b)./9;
    R = (2*a.^3 - 9.*a.*b + 27.*c)./54;
    
    M = R.^2 - Q.^3;
    r = nan(numel(c), 3);

    neg = M < 0;
    pos = ~neg;
    
    if any(neg)
        r(neg, :) = threeRoots(a(neg), R(neg), Q(neg));
    end
    if any(pos)
        r(pos, 1) = singleRoot(a(pos), R(pos), Q(pos), M(pos));
    end
end

function r = threeRoots(a, R, Q)
    theta = acos(R./sqrt(Q.^3));
    r1 = -(2.*sqrt(Q).*cos(theta/3)) - a./3;
    r2 = -(2.*sqrt(Q).*cos((theta + 2*pi)/3)) - a./3;
    r3 = -(2.*sqrt(Q).*cos((theta - 2*pi)/3)) - a./3;
    r = [r1, r2, r3];
end

function r = singleRoot(a, R, Q, M)
    S = -sign(R).*(abs(R) + sqrt(M)).^(1/3);
    T = Q./S;
    T(~isfinite(T)) = 0;
    r = S + T - a./3;
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
