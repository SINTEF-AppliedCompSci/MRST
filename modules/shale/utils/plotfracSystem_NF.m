function [hm,hf,hw]=plotfracSystem_NF(G,fracplanes,numHFplanes,wells,varargin)
% Plots fracplanes with matrix grid as backdrop
% This will create Hydraulic fracture,natural fractures and wells in
% different colors
% Author: Harun Rashid, Date modified: 01/03/2020
opt=struct('fracturelist',1:length(fracplanes),'label',false,'fracfacealpha',0.5,'fracfacecolor1','y','fracfacecolor2','b');
% opt=merge_options(opt,varargin{:});


hf={};
for i = opt.fracturelist(1:numHFplanes)
    X=fracplanes(i).points(:,1);
    Y=fracplanes(i).points(:,2);
    Z=fracplanes(i).points(:,3);
    hf{end+1}=fill3(X,Y,Z,opt.fracfacecolor1);
    hold on;
    set(hf{end},'facealpha',opt.fracfacealpha);
    
    if opt.label
        [~,index]=max(X);
        text(X(index),Y(index),Z(index),['\leftarrow Fracture ',num2str(i)],'Color','r');
    end
end
for i = opt.fracturelist(numHFplanes+1:end)
    X=fracplanes(i).points(:,1);
    Y=fracplanes(i).points(:,2);
    Z=fracplanes(i).points(:,3);
    hf{end+1}=fill3(X,Y,Z,opt.fracfacecolor2);
    hold on;
    set(hf{end},'facealpha',opt.fracfacealpha);
    
    if opt.label
        [~,index]=max(X);
        text(X(index),Y(index),Z(index),['\leftarrow Fracture ',num2str(i)],'Color','r');
    end
end
hm=plotGrid(G,'facealpha',0);

hw={};
for i = 1:numel(wells)
    if(isempty(fieldnames(wells))) %Empty wells
        break
    end
    X=wells(i).points(:,1);
    Y=wells(i).points(:,2);
    Z=wells(i).points(:,3);
    hw{end+1}=plot3(X,Y,Z,'r-','linewidth',2);
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
        X=Centers(:,1);
        Y=Centers(:,2);
        Z=Centers(:,3);
        scatter3(X,Y,Z,200,'Marker','s',...
            'MarkerEdgeColor',[0 0.4470 0.7410],...
            'MarkerFaceColor',[0 0.4470 0.7410]);
        hold on;
    end

end

view(15,20);
axis equal tight;
hold off;

end
