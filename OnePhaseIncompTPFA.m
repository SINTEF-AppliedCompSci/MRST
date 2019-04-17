function [state,matrixInfo] = OnePhaseIncompTPFA(G,T,src)
%Single phase incompressible flow using TPFA
%   Applicable to general 2D and 3D grids
%   Grid - Grid structure of MRST
%   T - half-face transmissibility matrix including boundary conditions
%   src - Source structure of MRST

nc=G.cells.num;
nf=G.faces.num;
N=G.faces.neighbors;
TT=zeros(nf,1);
ind=all(N~=0,2);
TT(ind)=T(ind,1).*T(ind,2)./(T(ind,1)+T(ind,2));
ncf=max(diff(G.cells.facePos));
% Assemble Ax=b------------------------------------------------------------
b=zeros(nc,1);
[I,J,V]=deal(zeros(ncf*nc,1));k=1;
for i_face=1:nf
    c1=N(i_face,1);c2=N(i_face,2);
    if(all([c1 c2]~=0))
        I(k)=c1;J(k)=c1;V(k)=TT(i_face);k=k+1;
        I(k)=c1;J(k)=c2;V(k)=-TT(i_face);k=k+1;
        I(k)=c2;J(k)=c2;V(k)=TT(i_face);k=k+1;
        I(k)=c2;J(k)=c1;V(k)=-TT(i_face);k=k+1;
    else
        c=max(c1,c2);
        I(k)=c;J(k)=c;V(k)=T(i_face,1);k=k+1;
        b(c)=b(c)-T(i_face,2);
    end
end
I(k:end)=[];J(k:end)=[];V(k:end)=[];
A=sparse(I,J,V,nc,nc);
b(src.cell(:))=b(src.cell(:))+src.rate(:);
if(condest(A)>1e10),A(1)=2*A(1);end
u=A\b;
if(nargout>1)
    matrixInfo.A=A;
    matrixInfo.b=b;
end
% Flux computation---------------------------------------------------------
flux=zeros(nf,1);
c1=N(ind,1);c2=N(ind,2);
flux(ind)=TT(ind).*(u(c1)-u(c2));
c=max(N(~ind,:),[],2);
flux(~ind)=T(~ind,1).*u(c)+T(~ind,2);
ind=N(:,1)==0;
flux(ind)=-flux(ind);

state.pressure=u;
state.flux=flux;
end

