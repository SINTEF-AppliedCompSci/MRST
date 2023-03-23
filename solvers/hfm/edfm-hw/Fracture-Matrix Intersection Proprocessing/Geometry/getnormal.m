function normal=getnormal(nodes,tol)
% get normal
dir1=nodes(2,:)-nodes(1,:);
dir2=nodes(3,:)-nodes(1,:);
% count=3;
normal=cross(dir1,dir2);
% while abs(norm(normal))<tol % parallel directions, use another vertex
%     count=count+1;
%     if count>size(nodes,1)
%         disp(nodes);
%         error('convexpolygonarea unable to find a normal to the polygon. Check nodes.')
%     end
%     dir2=nodes(count,:)-nodes(1,:);
%     normal=cross(dir1,dir2);
% end
normal=normal/norm(normal); % make unit

end