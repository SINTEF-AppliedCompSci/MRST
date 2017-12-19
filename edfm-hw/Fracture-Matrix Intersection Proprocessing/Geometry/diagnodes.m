function [nodepairs,direction] = diagnodes(nodes,tol)
% DIAGONAL NODES - Takes in nodes of a cuboid (each row containing x,y,z
% coordinates) and generates diagonal pairings in nodepairs (contains the
% row indices in nodes of diagonal pairs). Also generates direction which
% contains unit vectors in each row corresponding to nodepairs. The unit
% vector points from the first to second element in nodepairs.
%
% Diagonal pairings generated in this function always pass through the
% center of the cube.

numnodes=size(nodes,1);
nodepairs=zeros(4,2);
direction=zeros(4,3);
count=0;
for i=1:(numnodes-1)
    for j=(i+1):numnodes
        % if none of the x,y,z coordinates are the same, we have found our
        % diagonal node. Save the node pairs, calculate and save the unit
        % direction vector. Then break out of the current inner loop and 
        % proceed to the next node with the outer loop.
        xdifferent=abs(nodes(i,1)-nodes(j,1))>tol;
        ydifferent=abs(nodes(i,2)-nodes(j,2))>tol;
        zdifferent=abs(nodes(i,3)-nodes(j,3))>tol;
        if xdifferent && ydifferent && zdifferent
            count=count+1;
            nodepairs(count,:)=[i,j];
            dir=(nodes(j,:)-nodes(i,:))/norm(nodes(j,:)-nodes(i,:));
            direction(count,:)=dir;
            break;
        end
    end
end

end

