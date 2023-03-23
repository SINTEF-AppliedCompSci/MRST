function edges=findAABBedges(AABBnodes,tol)
% FIND AABB EDGES. Takes in the nodal positions of an AABB (each row
% containing xyz coordinates). The output is a 12x2 matrix with each row
% containing the row indices of vertices that form an edge.

% METHODOLOGY: We make use of the fact that an AABB is axis-aligned. So for
% every vertex pair that forms an edge, only one of the coordinates can be
% different. We start with vertex1, compare with vertices 2 through 8. Then
% move on to vertex2, compare with vertices 3 through 8. This continues
% until vertex7, which is compared with vertex8.

edges=zeros(12,2); % initialize edges
count=1;
for i=1:7
    xi=AABBnodes(i,1); yi=AABBnodes(i,2); zi=AABBnodes(i,3);
    for j=(i+1):8
        xj=AABBnodes(j,1); yj=AABBnodes(j,2); zj=AABBnodes(j,3);
        xsame=(abs(xi-xj)<tol);
        ysame=(abs(yi-yj)<tol);
        zsame=(abs(zi-zj)<tol);
        if (~xsame && ysame && zsame) || (xsame && ~ysame && zsame) || (xsame && ysame && ~zsame)
            edges(count,:)=[i,j];   
            count=count+1;
        end
    end
end


end

