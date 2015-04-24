function C=GenC(Grid)

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
