function plotfracongrid(G_matrix,fracplanes,varargin)
% Plots fracplanes with matrix grid as backdrop

opt=struct('fracturelist',1:length(fracplanes),'label',true,'fracfacealpha',0.5);
opt=merge_options(opt,varargin{:});



for i = opt.fracturelist
    X=fracplanes(i).points(:,1);
    Y=fracplanes(i).points(:,2);
    Z=fracplanes(i).points(:,3);
    h=fill3(X,Y,Z,'y');
    hold on;
    set(h,'facealpha',opt.fracfacealpha);
    
    if opt.label
        [~,index]=max(X);
        text(X(index),Y(index),Z(index),['\leftarrow Fracture ',num2str(i)],'Color','r');
    end
end

plotGrid(G_matrix,'facealpha',0);

view(15,20);
axis equal tight;
hold off;



end