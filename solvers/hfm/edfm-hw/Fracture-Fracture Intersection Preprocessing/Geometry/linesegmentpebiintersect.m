function truth=linesegmentpebiintersect(line,pebi,tol)
% INCOMPLETE. DON'T USE.
% Function checks if a line segment [x1 y1 z1; x2 y2 z2] intersects with a
% pebi grid. pebi is an nx3 array with nodes in each row.
disp('Using linesegmentpebiintersect. Code is incomplete. Careful!');
numnodes=size(pebi,1);

normal=getnormal(pebi,tol);

positivedir=cross(normal,line(2,:)-line(1,:));

relpos=pebi-repmat(line(1,:),numnodes,1);

loc=dot(relpos,repmat(positivedir,numnodes,1),2);

loc=loc./abs(loc);

check=abs(sum(loc));

if abs(check-numnodes)<tol
    truth=false;
else
    truth=true;
end

end