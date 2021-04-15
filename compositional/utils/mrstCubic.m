function r = mrstCubic(a, b, c, d)
% Straightforward implementation of a cubic root solver for vectorized
% problems. See https://en.wikipedia.org/wiki/Cubic_function
%
% Vectorized version of the Matlab builtin `roots` for cubic and lower
% order equations.
    if nargin == 1
        % We recieved a matrix. Unpack into vectors.
        b = a(:, 2);
        c = a(:, 3);
        d = a(:, 4);
        a = a(:, 1);
    end

    d0 = b.^2 - 3*a.*c;
    d1 = 2*b.^3 - 9*a.*b.*c + 27*d.*a.^2;
    
    C = ((0.5) *( d1 + sqrt(d1.^2 - 4*d0.^3))).^(1/3);
    
    % Corner case: Handle division by zero
    bad = C == 0;
    if any(bad)
        C(bad) = ((0.5) *( d1(bad) - sqrt(d1(bad).^2 - 4*d0(bad).^3))).^(1/3);
    end
    % Handle d == 0 and d1 == 0
    C1 = C;
    C1(d == 0 & d1 == 0) = 1;
    
    u1 = 1;
    u2 = (-1 + 1i*sqrt(3))/2;
    u3 = (-1 - 1i*sqrt(3))/2;
    
    % Roots parametrized by three inputs
    v = @(u) -(1./(3*a)).*(b + u.*C + d0./(u.*C1));
    r = [v(u1), v(u2), v(u3)];
    % Filter out any imaginary noise (assuming that the numbers are nicely
    % scaled)
    noise = abs(imag(r)) < 1e-12;
    r(noise) = real(r(noise));
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
