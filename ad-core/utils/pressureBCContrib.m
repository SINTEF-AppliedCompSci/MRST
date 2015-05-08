function  [bXqXbc,bc_cell] = pressureBCContrib(G,s,pX,rhoX,mobX,bX,bc)        
        assert(all(strcmp(bc.type,'pressure')),'only pressure bc allowed');
        Tbc=s.T_all(bc.face);
        assert(all(sum(G.faces.neighbors(bc.face,:)>0,2)==1),'bc on internal boundary');
        bc_cell=sum(G.faces.neighbors(bc.face,:),2);        
        dzbc=(G.cells.centroids(bc_cell,end)-G.faces.centroids(bc.face,end));
        g=norm(gravity);
        bXqXbc=cell(numel(pX),1);
        for i=1:numel(pX)
            pXbc=pX{i}(bc_cell);
            rhoXbc=rhoX{i}(bc_cell);
            dpbc=bc.value-pXbc+g*(rhoXbc.*dzbc);
            bXmobXbc = bX{i}(bc_cell).*mobX{i}(bc_cell);
            bXmobXbc(dpbc>0)=0;
            %NB some definition of in flow has to be done
            %if(any(dpbc_w>0))
            %    bWmobWbc(dpbc_w>0)=bW(bc_cell(dpbc_w>0)).*(mobW(bc_cell(dpbc_w>0))+mobO(bc_cell(dpbc_w>0)));     
            %end
            bXqXbc{i}  = -bXmobXbc.*Tbc.*(dpbc);
        end      
end
