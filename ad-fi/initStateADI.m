function state = initStateADI(G, f, deck,varargin)
%Undocumented Utility Function

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

opt = struct('dz', 50);
opt = merge_options(opt, varargin{:});
eql      = deck.SOLUTION.EQUIL;
max_z=max(G.cells.centroids(:,3));
min_z=min(G.cells.centroids(:,3));
state.s=nan(G.cells.num,3);
state.pressure=nan(G.cells.num,1);
op=nan(G.cells.num,1);
wp=nan(G.cells.num,1);
gp=nan(G.cells.num,1);
rs=nan(G.cells.num,1);
% set steping of generated tables
dz=opt.dz;
for nr=1:size(deck.SOLUTION.EQUIL,1)
    datDepth = eql(nr,1);
    datPress  = eql(nr,2);
    woc      = eql(nr,3);
    goc      = eql(nr,5);
    % fid pvt cell
    if(isfield(deck.REGIONS,'PVTNUM'))
        pvtnrs = deck.REGIONS.PVTNUM(deck.REGIONS.EQLNUM==nr);
        pvtnr = pvtnrs(1);
        cells=find(deck.REGIONS.PVTNUM(G.cells.indexMap)==pvtnr);
        cell=cells(1);
    else
        cell=1;
    end
    rsZ =@(z) interpTable(deck.SOLUTION.RSVD{nr}(:,1),deck.SOLUTION.RSVD{nr}(:,2), z);
    if(datDepth > goc   && datDepth <=woc)
        tabO=equilibriate_regionO(datDepth,datPress,f,rsZ,max_z,min_z, dz,cell);
        datDepthW= woc;
        datPressW= interpTable(tabO(:,1), tabO(:,2), woc)-eql(nr,4);
        tabW=equilibriate_regionW(datDepthW,datPressW,f,max_z,min_z, dz, cell);
        datDepthG= goc;
        datPressG= interpTable(tabO(:,1), tabO(:,2), goc)+eql(nr,7);
        tabG=equilibriate_regionG(datDepthG,datPressG,f,max_z,min_z, dz, cell);
    elseif(datDepth <= goc)
        tabG=equilibriate_regionG(datDepth,datPress,f,max_z,min_z, dz, cell);
        datDepthO=goc;
        datPressO= interpTable(tabG(:,1), tabG(:,2), goc)-eql(nr,7);
        tabO=equilibriate_regionO(datDepthO,datPressO,f,rsZ,max_z,min_z, dz, cell);

        datDepthW= woc;
        datPressW= interpTable(tabO(:,1), tabO(:,2), woc)-eql(nr,4);
        tabW=equilibriate_regionW(datDepthW,datPressW,f,max_z,min_z, dz, cell);
    elseif(datDepth>= woc)
        tabW=equilibriate_regionW(datDepth,datPress,f,max_z,min_z, dz, cell);
        datDepthO=woc;
        datPressO= interpTable(tabW(:,1), tabW(:,2), woc)+eql(nr,4);
        tabO=equilibriate_regionO(datDepthO,datPressO,f,rs,max_z,min_z, dz, cell);
        datDepthG= goc;
        datPressG= interpTable(tabO(:,1), tabO(:,2), wog)+eql(nr,7);
        tabG=equilibriate_regionG(datDepthG,datPressG,f,max_z,min_z, dz, cell);
    else
        error('No proper value for datDepth')
    end

    if(isfield(deck.REGIONS,'EQLNUM'))
        eqnum_cells=nr==deck.REGIONS.EQLNUM(G.cells.indexMap);
    else
        eqnum_cells=true(G.cells.num,1);
    end

    indo=eqnum_cells & G.cells.centroids(:,3)<woc & G.cells.centroids(:,3)>goc;
    indw=eqnum_cells & G.cells.centroids(:,3)>woc & G.cells.centroids(:,3)>goc;
    indg=eqnum_cells & G.cells.centroids(:,3)<woc & G.cells.centroids(:,3)<goc;
    op(eqnum_cells)=interpTable(tabO(:,1), tabO(:,2), G.cells.centroids(eqnum_cells,3));
    wp(eqnum_cells)=interpTable(tabW(:,1), tabW(:,2), G.cells.centroids(eqnum_cells,3));
    gp(eqnum_cells)=interpTable(tabG(:,1), tabG(:,2), G.cells.centroids(eqnum_cells,3));
    rs(eqnum_cells)= rsZ(G.cells.centroids(eqnum_cells,3));



    state.snc(indw,:) = repmat([1 0 0],sum(indw),1);
    state.snc(indo,:) = repmat([0 1 0],sum(indo),1);
    state.snc(indg,:) = repmat([0 0 1],sum(indg),1);



end
state.flux = zeros(G.faces.num, 1);
state.ph_press=[wp,op,gp];
% need sto limit if not gass or oil??

for i=1:max(deck.REGIONS.SATNUM)
    ind=i==deck.REGIONS.SATNUM(G.cells.indexMap);
    ind_cells=find(ind);
    if(sum(ind)>0)
        state.pressure(ind)=op(ind);

        tabW=deck.PROPS.SWOF{i};
        tabG=deck.PROPS.SGOF{i};
        dp=op(ind)-wp(ind);

        %sw=interpTable(tabW(end:-1:1,end),tabW(end:-1:1,1),dp);
        sw=interp1q(tabW(end:-1:1,end),tabW(end:-1:1,1),dp);

        sw(dp>tabW(1,end))=tabW(1,1);
        sw(dp<tabW(end,end))=tabW(end,1);
        % there is not flowing oil or gas
        state.pressure(ind_cells(dp<tabW(end,end)))=wp(ind_cells(dp<tabW(end,end)))+tabW(end,end);

        dp =gp(ind)-op(ind);
        %sg=interpTable(tabG(:,end),tabG(:,1),op(ind)-gp(ind));
        sg=interp1q(tabG(:,end),tabG(:,1),op(ind)-gp(ind));
        sg(dp<tabG(1,end))=tabG(1,1);
        sg(dp>tabG(end,end))=tabG(end,1);

        % there is no flowing oil or water
        state.pressure(ind_cells(dp>tabG(end,end)))=gp(ind_cells(dp>tabG(end,end)))+tabG(end,end);

        % correct ?
        sw=max(min(tabW(:,1)),sw);
        sw=min(max(tabW(:,1)),sw);

        sg=max(min(tabG(:,1)),sg);
        sg=min(max(tabG(:,1)),sg);
        swg=sw+sg;
        swg(swg>1)=1;
        sg(swg==1)=1-sw(swg==1);

        %sw(swg==1)=1-sg(swg==1);
        state.s(ind,:)=[sw, 1-(sw+sg), sg];

    end
end
rs_max=f.rsSat(state.pressure);
rs=min(rs,rs_max);

rs(state.s(:,2)==0)=rs_max(state.s(:,2)==0);
state.rs = rs;
end
function tab=equilibriate_regionO(datDepth,datPress,f, rsZ, max_z, min_z, dz, cell)

    g = norm(gravity);
    isSat=1;
    rsSat =@(z,p) max(0,min(rsZ(z),f.rsSat(p,'cellInx',1)));
    grho =@(z,p) g*f.rhoOS./f.BO(p,rsSat(z,p),isSat,'cellInx',cell);

    tab = makePressureTab(dz,max_z,min_z,datDepth,datPress, grho);


end
function tab=equilibriate_regionW(datDepth,datPress,f, max_z, min_z, dz, cell)

    g = norm(gravity);
    grho =@(z, p) g*f.rhoWS./f.BW(p,'cellInx', cell);

    tab = makePressureTab(dz,max_z,min_z,datDepth,datPress, grho);

end
function tab=equilibriate_regionG(datDepth,datPress,f, max_z, min_z, dz, cell)
    g = norm(gravity);

    grho =@(z,p) g*f.rhoGS./f.BG(p,'cellInx', cell);

    tab = makePressureTab(dz,max_z,min_z,datDepth,datPress, grho);


end
function tab = makePressureTab(dz,max_z,min_z,datDepth,datPress,grho)
    tab=[];
     options = odeset('OutputFcn',@odeplot)
    %options=odeset();
    if(datDepth<max_z)

        [zd,pd]=ode45(@(z,p) grho(z,p),[0:dz:(max_z+dz-datDepth)], datPress,options);%#ok
        tab=[datDepth+zd,pd];
    end
    if(datDepth>min_z)
        %dz=min((datDepth-min_z),dz);
        [zup,pup]=ode45(@(z,p) -grho(z,p),[0:dz:(datDepth+dz-min_z)], datPress,options);%#ok
        if(~isempty(tab))
            tab=[datDepth-zup(end:-1:2),pup(end:-1:2);tab];
        else
            tab=[datDepth-zup(end:-1:1),pup(end:-1:1)];
        end
    end
end
