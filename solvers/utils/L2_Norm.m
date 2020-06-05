function [ep,ef] = L2_Norm(G,ua,fa,uh,fh)
%Compute L2 norm of pressure errors and flux errors. used for convergence
%study
%   G - Grid structure of MRST
%   ua - analytical pressure solution at cell centroids
%   fa - analytical flux solution at face centroids
%   uh - numerical pressure solution at cell centroids
%   fh - numerical flux solution at face centroids
%   ep - L2 norm of pressure error
%   ef - L2 nomr of flux error

A=G.cells.volumes;
ep=sqrt(sum(A.*(ua-uh).^2)/sum(A));

fa=abs(fa./G.faces.areas);
fh=abs(fh./G.faces.areas);
A=zeros(G.faces.num,1);
c1=G.faces.neighbors(:,1);
c2=G.faces.neighbors(:,2);
ind=all(G.faces.neighbors~=0,2);
A(ind)=0.5*(G.cells.volumes(c1(ind))+G.cells.volumes(c2(ind)));
A(~ind)=G.cells.volumes(max(G.faces.neighbors(~ind,:),[],2));
ef=sqrt(sum(A.*(fa-fh).^2)/sum(A));

% A=G.cells.volumes;
% ep=sqrt(sum(A.*(ua-uh).^2)/sum(A.*ua.^2));
% fa=abs(fa./G.faces.areas);
% fh=abs(fh./G.faces.areas);
% A=zeros(G.faces.num,1);
% c1=G.faces.neighbors(:,1);
% c2=G.faces.neighbors(:,2);
% ind=all(G.faces.neighbors~=0,2);
% A(ind)=0.5*(G.cells.volumes(c1(ind))+G.cells.volumes(c2(ind)));
% A(~ind)=G.cells.volumes(max(G.faces.neighbors(~ind,:),[],2));
% ef=sqrt(sum(A.*(fa-fh).^2)/sum(A.*fa.^2));
end
