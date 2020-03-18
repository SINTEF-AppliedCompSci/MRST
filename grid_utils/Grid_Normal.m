function [Grid]=Grid_Normal(Lx,Ly,Nx,Ny,perturb,isMove)
%Generate normal quadrilateral grid grid
% isMove: whether to perturb mid-lines
%      'moveYes': perturb mid-lines
%      'moveNo': do not perturb mid-lines

Grid=cartGrid([Nx Ny],[Lx Ly]);
c=Grid.nodes.coords;
if(strcmpi(isMove,'moveYes'))
    ind=any(c==0,2)|any(c(:,1)==Lx,2)|any(c(:,2)==Ly,2);
else
    ind=any(c==0,2)|any(c(:,1)==Lx,2)|any(c(:,2)==Ly,2)|...
        any(c(:,1)==0.5*Lx,2)|any(c(:,2)==0.5*Ly,2);
end
rn=-0.5+rand(sum(~ind),2);
c(~ind,1)=c(~ind,1)+perturb*rn(:,1)*Lx/Nx;
c(~ind,2)=c(~ind,2)+perturb*rn(:,2)*Ly/Ny;
Grid.nodes.coords=c;
Grid=computeGeometry(Grid);
end