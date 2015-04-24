function B=GenB(Grid,K)
%Undocumented function

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

Nx=Grid.Nx; Ny=Grid.Ny; Nz=Grid.Nz; N=Nx*Ny*Nz;
hx=Grid.hx; hy=Grid.hy; hz=Grid.hz;
L = K.^(-1);
ex=N-Ny*Nz; ey=N-Nx*Nz; ez=N-Nx*Ny;                      % Number of edges
tx=hx/(6*hy*hz); ty=hy/(6*hx*hz); tz=hz/(6*hx*hy);       % Transmissibilities

X1=zeros(Nx-1,Ny,Nz); X2=zeros(Nx-1,Ny,Nz);              % Preallocate memory
X0=L(1,1:Nx-1,:,:)+L(1,2:Nx,:,:);   x0=2*tx*X0(:);       % Main diagonal
X1(2:Nx-1,:,:)=L(1,2:Nx-1,:,:);     x1=tx*X1(:);         % Upper diagonal
X2(1:Nx-2,:,:)=L(1,2:Nx-1,:,:);     x2=tx*X2(:);         % Lower diagonal

Y1=zeros(Nx,Ny-1,Nz); Y2=zeros(Nx,Ny-1,Nz);              % Preallocate memory
Y0=L(2,:,1:Ny-1,:)+L(2,:,2:Ny,:); y0=2*ty*Y0(:);         % Main diagonal
Y1(:,2:Ny-1,:)=L(2,:,2:Ny-1,:);   y1=ty*Y1(:);           % Upper diagonal
Y2(:,1:Ny-2,:)=L(2,:,2:Ny-1,:);   y2=ty*Y2(:);           % Lower diagonal

Lz1=L(3,:,:,1:Nz-1); z1=tz*Lz1(:);                       % Upper diagonal
Lz2=L(3,:,:,2:Nz);   z2=tz*Lz2(:);                       % Lower diagonal
z0=2*(z1+z2);                                            % Main diagonal

B=[spdiags([x2,x0,x1],[-1,0,1],ex,ex),sparse(ex,ey+ez);...
   sparse(ey,ex),spdiags([y2,y0,y1],[-Nx,0,Nx],ey,ey),sparse(ey,ez);...
   sparse(ez,ex+ey),spdiags([z2,z0,z1],[-Nx*Ny,0,Nx*Ny],ez,ez)];
