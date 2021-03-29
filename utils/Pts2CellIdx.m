function [Ids]=Pts2CellIdx(Pts,G)
    %Convert Points coord into cartesian gridblock global index
    Pts=Pts*1.0000001; %small pertubation to avoid Pts just on the face
    
    NX=G.cartDims(1);
    NY=G.cartDims(2);
    if(numel(G.cartDims)==3), NZ=G.cartDims(3);end
    %Row order in tensor grid
    G_coordX=G.nodes.coords(1:NX+1,1); 
    G_coordY=G.nodes.coords(1:NX+1:(NX+1)*(NY+1),2);
    if(numel(G.cartDims)==3)
        G_coordZ=G.nodes.coords(1:(NX+1)*(NY+1)+1:(NX+1)*(NY+1)*(NZ+1),3);
    end
    
    %Find the index
    Ids=zeros(size(Pts,1),numel(G.cartDims));
    for i=1:size(Pts,1)
        I=find(G_coordX> Pts(i,1)); Ids(i,1)=min(I)-1;
        J=find(G_coordY> Pts(i,2)); Ids(i,2)=min(J)-1;
        if(numel(G.cartDims)==3)
          K=find(G_coordZ> Pts(i,3)); Ids(i,3)=min(K)-1;
        end
    end
    
    %Convert IJK into global 1D index
    if(numel(G.cartDims)==3)
        Ids = sub2ind(G.cartDims, Ids(:,1), Ids(:,2),Ids(:,3));
    else
        Ids = sub2ind(G.cartDims, Ids(:,1), Ids(:,2));
    end
end