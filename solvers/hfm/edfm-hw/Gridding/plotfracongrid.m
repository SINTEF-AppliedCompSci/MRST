function [hm,hf]=plotfracongrid(G_matrix,fracplanes,varargin)
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

hm=plotGrid(G_matrix,'facealpha',0);

view(15,20);
axis equal tight;
hold off;



end