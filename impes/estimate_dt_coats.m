function dt=estimate_dt_coats(G,trans,pv,mob,dmob,p,rho,pc,dpc,clf_fac)
% Implemantation of IMPES stability step from SPE 69225 written by Coats.
% For three phases we have assumed %aqua   = 1;liquid = 2; vapor  = 3;
% For two phases we have assumed %liquid = 1; vapor  = 2;
% Implemention could been made simpler with out assuming any of the
% assumption of black oil model i.e. water do not mix with gas.

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


  if(size(mob,2)==2)
     pc=zeros(numel(p),2);
     dpc=pc;
     phi=repmat(p,1,2)-pc;%+bsxfun(@times,rho,sum(bsxfun(@times,G.cells.centroids,gravity(),2)));
     internal=sum(G.faces.neighbors>0,2)==2;
     cell1=G.faces.neighbors(internal,1);cell2=G.faces.neighbors(internal,2);
     dphi=phi(cell1,:)-phi(cell2,:);
     ff_int=  trans(internal).*...
        ((mob(internal,1).*dmob(internal,4).*abs(dphi(:,2))...
        -mob(internal,2).*(-dmob(internal,1)).*abs(dphi(:,1))...%forteikn?
        -mob(internal,2).*mob(internal,1).*(dpc(cell1,2)+dpc(cell2,2))...
        )./sum(mob(internal,:),2));
        %-mob(internal,2).*mob(internal,1).*(dpc(2,cell1)+dpc(2,cell2))...%derivitives of cappilary terms eq 6
     nf     = diff(G.cells.facePos);
     cellno = rldecode(1 : G.cells.num, nf, 2) .';
     cellfaces=G.cells.faces(:,1);
     ff=zeros(G.faces.num,1);
     ff(internal)=ff_int;
     F=accumarray(double(cellno),ff(cellfaces),[double(G.cells.num),1]);
     assert(all(F>=0));
     dt=clf_fac.*pv./F;
     dt=min(dt);
  elseif(size(mob,2)==3)
     %error('three phases not implemented jet')
     pc=zeros(numel(p),3);
     dpc=pc;
     phi=repmat(p,1,3)-pc;%+bsxfun(@times,rho,sum(bsxfun(@times,G.cells.centroids,gravity(),2)));
     internal=sum(G.faces.neighbors>0,2)==2;
     cell1=G.faces.neighbors(internal,1);cell2=G.faces.neighbors(internal,2);
     dphi=phi(cell1,:)-phi(cell2,:);
     nf     = diff(G.cells.facePos);
     cellno = rldecode(1 : G.cells.num, nf, 2) .';
     cellfaces=G.cells.faces(:,1);
     %% f22
     ff22_int=  trans(internal).*...
        (((mob(internal,1)+mob(internal,2)).*dmob(internal,9).*abs(dphi(:,3))...
        -mob(internal,3).*(-dmob(internal,8)).*abs(dphi(:,2))...%forteikn?
        -mob(internal,3).*dmob(internal,7).*abs(dphi(:,1))...
        +mob(internal,3).*(mob(internal,2)+mob(internal,1)).*(dpc(cell1,3)+dpc(cell2,3))...
        ));%./sum(mob(internal,:),2));
        %-mob(internal,2).*mob(internal,1).*(dpc(2,cell1)+dpc(2,cell2)-dpc(
        %1,cell1)+dpc(1,cell2))...%derivitives of cappilary terms eq 6
     no_zero=~(ff22_int==0);
     tmob_int=sum(mob(internal,:),2);
     ff22_int(no_zero)=ff22_int(no_zero)./tmob_int(no_zero);
     ff22=zeros(G.faces.num,1);
     ff22(internal)=ff22_int;
     F22=accumarray(double(cellno),ff22(cellfaces),[double(G.cells.num),1]);
     %% f11
     ff11_int=  trans(internal).*...
        (((mob(internal,2)+mob(internal,3)).*dmob(internal,1).*abs(dphi(:,1))...
        -mob(internal,1).*(-dmob(internal,2)).*abs(dphi(:,2))...%forteikn?
        +mob(internal,1).*(mob(internal,2)+mob(internal,3)).*(dpc(cell1,1)+dpc(cell2,1))...
        ));%./sum(mob(internal,:),2));
     no_zero=~(ff11_int==0);
     ff11_int(no_zero)=ff11_int(no_zero)./tmob_int(no_zero);
        %-mob(internal,2).*mob(internal,1).*(dpc(2,cell1)+dpc(2,cell2)-dpc(
        %1,cell1)+dpc(1,cell2))...%derivitives of cappilary terms eq 6
     ff11=zeros(G.faces.num,1);
     ff11(internal)=ff11_int;
     F11=accumarray(double(cellno),ff11(cellfaces),[double(G.cells.num),1]);
     %% f12
     ff12_int=  -trans(internal).*...
        ((mob(internal,1).*dmob(internal,8).*abs(dphi(:,2))...
        +mob(internal,1).*(dmob(internal,9)).*abs(dphi(:,3))...%forteikn?
        -(mob(internal,2)+mob(internal,3)).*dmob(internal,7).*abs(dphi(:,1))...
        +mob(internal,1).*mob(internal,3).*(dpc(cell1,3)+dpc(cell2,3))...
        ));%./sum(mob(internal,:),2));
        %-mob(internal,2).*mob(internal,1).*(dpc(2,cell1)+dpc(2,cell2)-dpc(
        %1,cell1)+dpc(1,cell2))...%derivitives of cappilary terms eq 6
     no_zero=~(ff12_int==0);
     ff12_int(no_zero)=ff12_int(no_zero)./tmob_int(no_zero);
     ff12=zeros(G.faces.num,1);
     ff12(internal)=ff12_int;
     F12=accumarray(double(cellno),ff12(cellfaces),[double(G.cells.num),1]);
     %% f21
     ff21_int=  -trans(internal).*...
        ((mob(internal,3).*dmob(internal,1).*abs(dphi(:,1))...
        +mob(internal,3).*(-dmob(internal,2)).*abs(dphi(:,2))...%forteikn?
        -mob(internal,3).*mob(internal,1).*(dpc(cell1,3)+dpc(cell2,3))));
     no_zero=~(ff21_int==0);
     ff21_int(no_zero)=ff21_int(no_zero)./tmob_int(no_zero);
        %-mob(internal,2).*mob(internal,1).*(dpc(2,cell1)+dpc(2,cell2)-dpc(
        %1,cell1)+dpc(1,cell2))...%derivitives of cappilary terms eq 6
     ff21=zeros(G.faces.num,1);
     ff21(internal)=ff21_int;
     F21=accumarray(double(cellno),ff21(cellfaces),[double(G.cells.num),1]);
     assert(all(F11>=0));
     assert(all(F22>=0));
     fdet=(F11+F22).^2-4*(F11.*F22-F12.*F21);
     assert(all(fdet>-eps));
     fdet=min(0,fdet);
     F=0.5*abs(F11+F22+sqrt(fdet));
     dt=clf_fac.*pv./F;
     dt=min(dt);
  elseif(size(mob,2)==1)
     dt=inf;
  else
     error(['Fluid dimension is ',num2str(size(mob,2))])
  end
end



