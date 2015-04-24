[Nx,Ny,Nz] = deal(4,3,3);
Grid.Nx=Nx; Grid.hx=1/(Grid.Nx+1);       % Dimension in x-direction
Grid.Ny=Ny; Grid.hy=1/(Grid.Ny+1);       % Dimension in x-direction
Grid.Nz=Nz; Grid.hz=1/(Grid.Nz+1);       % Dimension in x-direction

K=ones(3,Nx,Ny,Nz);                      % Make unit permeability tensor

B=GenB(Grid,K);                          % Compute B-block of matrix
C=GenC(Grid);                            % Compute C-blocks of matrix
N = Nx*Ny*Nz;
A=[B,C';-C,sparse(N,N)];                % Assemble matrix

%%
spy(A,7);
hold on
m = size(A,1); n=size(B,1);
plot([0 n+.5; m n+.5], [n+.5 0; n+.5 m],'k-');
hold off
set(gca,'XTick',[],'YTick',[]); xlabel([]);

%%
spy(B,7);
hold on
mx = (Nx-1)*Ny*Nz+.5;
plot([.5 mx mx .5 .5],[.5 .5 mx mx .5],'-k');
my = Nx*(Ny-1)*Nz + mx;
plot([mx my my mx mx],[mx mx my my mx],'-k');
mz = Nx*Ny*(Nz-1) + my;
plot([my mz mz my my],[my my mz mz my],'-k');
hold off
set(gca,'XTick',[],'YTick',[]); xlabel([]);

%%
spy(C,7);
hold on
mx = (Nx-1)*Ny*Nz+.5;
my = Nx*(Ny-1)*Nz + mx;
plot([mx my; mx my], [0 0; N+1 N+1],'-k');
hold off
set(gca,'XTick',[],'YTick',[]); xlabel([]);

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
