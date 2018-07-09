function [truth,line]=truncatecolinearlines(line1,line2,tol)
% This function takes colinear lines line1 and line2 and determines the
% intersection line between the two. The input lines are 2x3 matrices with
% each row representing an endpoint's x,y,z coordinates.
%
% If lines do not overlap, truth==false and line=[].
% If lines overlap, truth==true and line is a 2x3 matrix containing
% endpoint coordinates.
%
% Note: error returned if lines are not co-linear, or if one of the lines
% has zero length (meaning there is either no intersection or a point
% intersection, both of which have no contribution to flow)



refpoint=line1(1,:);
refdir=(line1(2,:)-line1(1,:))/norm(line1(2,:)-line1(1,:)); % reference direction vector

% Line 1 endpoints
P1_1=line1(1,:);
P2_1=line1(2,:);
length1=norm(P2_1-P1_1);

% Line 2 endpoints
P1_2=line2(1,:);
P2_2=line2(2,:);
length2=norm(P2_2-P1_2);

% pre-process points
r1_1=(P1_1-refpoint)/refdir;
r2_1=(P2_1-refpoint)/refdir;
r1_2=(P1_2-refpoint)/refdir;
r2_2=(P2_2-refpoint)/refdir;

rline1=[r1_1;r2_1];
rline2=[r1_2;r2_2];

%% Check if lines have zero length

if length1<tol || length2<tol
    error('One of the lines has zero length!');
end

%% Check if co-linear by using area calculation
check1=heronsformula(P1_1,P2_1,P1_2,tol);
check2=heronsformula(P1_1,P2_1,P2_2,tol);

if ~(check1<tol && check2<tol)
    disp(['Line 1 is [',num2str(P1_1),']',' [',num2str(P1_2),'].']);
    disp(['Line 2 is [',num2str(P2_1),']',' [',num2str(P2_2),'].']);
    error('Lines are not co-linear!');
end

%% Locate intersection line endpoints along reference direction vector

% Check if no overlap
if max(rline1)<min(rline2) || min(rline1)>max(rline2)
    truth=false;
    line=[];
    return;
end

rcombi=[rline1;rline2];
lcombi=[line1;line2];

[~,ind]=sort(rcombi);
lcombi=lcombi(ind,:);

truth=true;
line=lcombi(2:3,:);
line=uniquetol(line,tol,'ByRows',true,'DataScale',1); % if point overlap, then line will just have one output

end

