function C=GenC(Grid)
%Undocumented function

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

Nx=Grid.Nx; Ny=Grid.Ny; Nz=Grid.Nz;
C=sparse(0,0);                                             % Empty sparse matrix
Nxy=Nx*Ny; N=Nxy*Nz;                                       % Number of grid-points
vx=ones(Nx,1); vy=ones(Nxy,1); vz=ones(N,1);               % Diagonals

for i=1:Ny*Nz                                              % vx-block of C
  Cx=spdiags([vx,-vx],[-1,0]-(i-1)*Nx,N,Nx-1);             % create bidiagonal block
  C=[C,Cx];                                                % append to C
end

for i=1:Nz                                                 % vy-block of C
  Cy=spdiags([vy,-vy],[-Nx,0]-(i-1)*Nxy,N,Nxy-Nx);         % create bidiagonal block
  C=[C,Cy];                                                % append to C
end

C = [C, spdiags([vz,-vz],[-Nxy,0],N,N-Nxy)];               % vz-block of C
