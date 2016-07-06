%% Two-point inner product
% Illustrate the choice of inner product for a simple rectangular cell
% [-1,1]x[-1,1]
K   = eye(2);
C   = [-1 0; 1 0; 0 -1; 0 1]; 
N   = 2*C; 
vol = 4;

% Two-point transmissibility matrix
T = diag(diag(N*K*C')./sum(C.*C,2))

% Inverse mimetic inner product: diagonal tensor
for i=1:2
   if i==1,
      K = eye(2);
   else
      K = [1 .5; .5 1];
   end
   Q = orth(C);
   P = eye(size(C,1)) - Q*Q';
   W = (N * K * N') ./ vol                                     %#ok<*NOPTS>
   St = diag(diag(W));
   P*St*P
   T = W + 2*P*St*P
end

%% Raviart-Thomas inner product
C   = .5*[-1 0; 1 0; 0 -1; 0 1];
N   = 2*C; 
vol = 1;
K   = eye(2)

%
Q  = orth(N);
P  = eye(size(C,1)) - Q*Q';
Sm = diag(1 ./ diag(N*K*N'))*vol;
M1 = C*(K\C')./vol
M2 = P*Sm*P

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
