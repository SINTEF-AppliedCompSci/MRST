function [X,Y, segInt] = lineLineInt(L1, L2)
% Calculates the intersections of line segments L1 with line segments L2
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
TOL = 1000 * (max(L1(:)) - min(L1(:))) * eps;
numL2 = size(L2,1);
numL1 = size(L1,1);
L1    = repmat(L1,  1 ,numL2);
L2    = repmat(L2', numL1, 1); 

L1x1  = L1(:,1:4:end);
L1x2  = L1(:,3:4:end);
L1y1  = L1(:,2:4:end);
L1y2  = L1(:,4:4:end);

L2x1  = L2(1:4:end,:);
L2x2  = L2(3:4:end,:);
L2y1  = L2(2:4:end,:);
L2y2  = L2(4:4:end,:);

A1    = L1y2 - L1y1;
B1    = L1x1 - L1x2;
C1    = A1.*L1x1 + B1.*L1y1;
A2    = L2y2 - L2y1;
B2    = L2x1 - L2x2;
C2    = A2.*L2x1 + B2.*L2y1;

det   = A1.*B2 - A2.*B1;

X     = (B2.*C1 - B1.*C2)./det;
Y     = (A1.*C2 - A2.*C1)./det;

segInt = min(L1x1, L1x2) - X<=TOL & X - max(L1x1, L1x2)<=TOL ...
       & min(L2x1, L2x2) - X<=TOL & X - max(L2x1, L2x2)<=TOL ...
       & min(L1y1, L1y2) - Y<=TOL & Y - max(L1y1, L1y2)<=TOL ...
       & min(L2y1, L2y2) - Y<=TOL & Y - max(L2y1, L2y2)<=TOL;
end