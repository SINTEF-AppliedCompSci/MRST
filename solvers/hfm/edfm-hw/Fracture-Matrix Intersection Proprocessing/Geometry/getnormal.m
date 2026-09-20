function normal=getnormal(nodes,tol)
% get unit normal of a polygon from its first three vertices.
%
% REVISED: the original had its collinear-vertex guard commented out, so if
% nodes 1-3 were collinear the cross product was zero and the normalised
% result was NaN, which then propagated silently. The guard is restored: it
% advances to further vertices until a non-degenerate normal is found.

dir1=nodes(2,:)-nodes(1,:);
dir2=nodes(3,:)-nodes(1,:);
count=3;
normal=cross(dir1,dir2);
while abs(norm(normal))<tol   % nodes 1..count collinear, try the next vertex
    count=count+1;
    if count>size(nodes,1)
        disp(nodes);
        error('getnormal: unable to find a normal to the polygon. Check nodes.')
    end
    dir2=nodes(count,:)-nodes(1,:);
    normal=cross(dir1,dir2);
end
normal=normal/norm(normal); % make unit

end
