function [hm,hf,hw]=plotFracSystem(G,fracplanes,wells,varargin)
% Plots fracplanes with matrix grid as backdrop

    opt=struct('fracturelist',1:length(fracplanes),'label',true,'fracfacealpha',0.5,'fracfacecolor','y');
    opt=merge_options(opt,varargin{:});


    hf={};
    for i = opt.fracturelist
        X=fracplanes(i).points(:,1);
        Y=fracplanes(i).points(:,2);
        Z=fracplanes(i).points(:,3);
        hf{end+1}=fill3(X,Y,Z,opt.fracfacecolor);
        hold on;
        set(hf{end},'facealpha',opt.fracfacealpha);

        if opt.label
            [~,index]=max(X);
            text(X(index),Y(index),Z(index),['\leftarrow Fracture ',num2str(i)],'Color','r');
        end
    end

    hm=plotGrid(G,'facealpha',0.0,'edgealpha',0.5);
    
    hw={};
    for i = 1:numel(wells)
        if(isempty(fieldnames(wells))) %Empty wells
            break
        end
        X=wells(i).points(:,1);
        Y=wells(i).points(:,2);
        Z=wells(i).points(:,3);
        hw{end+1}=plot3(X,Y,Z,'b-','linewidth',2);
        hold on;
        %Well Frac Intersection
        if(isfield(wells(i),'WellXFracPts') && numel(wells(i).WellXFracPts)>0)
            X=wells(i).WellXFracPts(:,1);
            Y=wells(i).WellXFracPts(:,2);
            Z=wells(i).WellXFracPts(:,3);
            scatter3(X,Y,Z,'MarkerEdgeColor','k',...
                'MarkerFaceColor',[0 .75 .75]);
            hold on;
        end
        
        %Intersected Frac Cell Centroid
        if(isfield(wells(i),'XFracCellIDs') && numel(wells(i).XFracCellIDs)>0)
            Centers=[];
            XFracCellIDs_global=wells(i).XFracCellIDs(:,1);
            XFracIDs=wells(i).XFracCellIDs(:,2);
            for i = 1:numel(XFracIDs)
                Frac=G.FracGrid.(['Frac' num2str(XFracIDs(i))]);
                CellID_local=XFracCellIDs_global(i)-Frac.cells.start+1;
                Centers=[Centers; Frac.cells.centroids(CellID_local,:)];
            end
            Centers;
            X=Centers(:,1);
            Y=Centers(:,2);
            Z=Centers(:,3);
            scatter3(X,Y,Z,100,'Marker','s',...
                'MarkerEdgeColor',[0 0.4470 0.7410],...
                'MarkerFaceColor',[0 0.4470 0.7410]);
            hold on;
        end
        
    end




    view(15,20);
    axis equal tight;
    hold off;



end