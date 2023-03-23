function [maxdiag,mindiag,maxdiaglength,mindiaglength]=findquaddiag3D(quadpoints)
% finds the row indices in quadpoints that form the max and min diagonals.
% Also calculates the max and min diagonal lengths.

possiblediags=[1,2;
               1,3;
               1,4;
               2,3;
               2,4;
               3,4];
possiblediaglength=zeros(6,1);

for i=1:6
    point1=quadpoints(possiblediags(i,1),:);
    point2=quadpoints(possiblediags(i,2),:);
    possiblediaglength(i)=norm(point1-point2);
end

[maxdiaglength,maxind]=max(possiblediaglength);

maxdiag=possiblediags(maxind,:);

temp1=and((possiblediags(:,1)~=maxdiag(1)),(possiblediags(:,1)~=maxdiag(2)));
temp2=and((possiblediags(:,2)~=maxdiag(1)),(possiblediags(:,2)~=maxdiag(2)));
temp=(temp1 & temp2);
minind=find(temp);
mindiag=possiblediags(minind,:);
mindiaglength=possiblediaglength(minind,:);

end