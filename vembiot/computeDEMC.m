function [CC,op] = computeDEMC(G,C,varargin)
%[CC,op] = computeDEMC(G,C,varargin)
% compute mdem CC intended to use to modify C to make DEM type hooks
% parameters

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

opt=struct('cells',1:G.cells.num,'use_dem',false);
opt=merge_options(opt,varargin{:});

rcell=true(G.cells.num,1);
rcell(opt.cells)=false;
if(any(~rcell))
    [G,cellmap,facemap,nodemap]=removeCells(G,find(rcell));
    
    G=computeGeometry(G);
    G=mrstGridWithFullMappings(G);
else
   nodemap=[1:G.nodes.num]';%#ok
end
[S,extra] =VEM_mrst_vec(G,C(opt.cells,:));%#ok
% mdem



%% 2D
node_num=G.nodes.num;
lindim=3;
u={};
u{1}=[G.nodes.coords(:,1),zeros(node_num,1)];
u{2}=[zeros(node_num,1),G.nodes.coords(:,2)]; 
u{3}=[G.nodes.coords(:,2),G.nodes.coords(:,1)];
% rot
u{4}=[ones(node_num,1),zeros(node_num,1)];
u{5}=[zeros(G.nodes.num,1),ones(node_num,1)]; 
u{6}=[G.nodes.coords(:,2),-G.nodes.coords(:,1)];
%%
%uu=u{4};
nfaces=G.cells.faces(:,1);
ifnodes=mcolon(G.faces.nodePos(nfaces),G.faces.nodePos(nfaces+1)-1);
ee=reshape(G.faces.nodes(ifnodes),2,[])';
barvec=G.nodes.coords(ee(:,2),:)-G.nodes.coords(ee(:,1),:);
el=sqrt(sum(barvec.^2,2));
barvecn=bsxfun(@rdivide,barvec,el);
%%
%
u2uu =@(u) reshape(u',[],1);
node2dofs =@(nodes) [mcolon(G.griddim*(nodes-1)+1,G.griddim*(nodes-1)+G.griddim)]';%#ok
reldis=sparse(node2dofs([1:size(ee,1)]),node2dofs(ee(:,2)),1,lindim*G.griddim*G.cells.num,G.nodes.num*G.griddim)...
    +sparse(node2dofs([1:size(ee,1)]),node2dofs(ee(:,1)),-1,lindim*G.griddim*G.cells.num,G.nodes.num*G.griddim);%#ok
[i,j]=blockDiagIndex(repmat(1,size(ee,1),1),repmat(G.griddim,size(ee,1),1));%#ok
rdis2rlen=sparse(i,j,reshape(barvecn',[],1));
%%
M=nan(lindim*G.cells.num,lindim);
u2rlen=rdis2rlen*reldis;
for i=1:3
    M(:,i)=u2rlen*u2uu(u{i});
end
%%
[i,j]=blockDiagIndex(repmat(lindim,G.cells.num,1),repmat(lindim,G.cells.num,1));
Mbl=sparse(i,j,reshape(M',[],lindim*lindim))';
%%
invM=invv(reshape(M',[],1),repmat(lindim,G.cells.num,1));
invMbl=sparse(i,j,reshape(invM',[],lindim*lindim))';
%%
k_mdem=invMbl'*extra.D*invMbl;
if(opt.use_dem)
   k_mdem=diag(diag(k_mdem)); 
end    
%{
if(break_edge & false)
   assert(all(barvecn(1,:)==[1 0 ]));
   kd =diag(k_mdem);
   kd(1)=0;
   k_dem=diag(kd);
else
    k_dem=diag(diag(k_mdem));
end
%k_dem=k_mdem;
D_dem=Mbl'*k_dem*Mbl
%}
%%
D_dem=Mbl'*k_mdem*Mbl;
[i,j]=blockDiagIndex(repmat(lindim,G.cells.num,1),repmat(lindim,G.cells.num,1));
ind=sub2ind(size(D_dem),i,j);
DD=reshape(D_dem(ind),lindim^2,[])';
%full(DD)
%%
%Dmat(k,:) = reshape(D_dem,1,[]);
CC=C2D(full(DD),G,'inv',true);
op = struct('u2rlen',u2rlen,'k_mdem',k_mdem,'ind',ind,'cells',opt.cells,'Mbl',Mbl,'invMbl',invMbl,'lindim',lindim,'C',CC,'nodemap',nodemap,'cellmap',cellmap,'facemap',facemap);
