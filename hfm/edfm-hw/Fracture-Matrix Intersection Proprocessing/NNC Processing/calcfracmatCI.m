function [CI,davg]=calcfracmatCI(nodes,area,planenormal,planepoint,davg,tol)
% CI CALCULATOR - Calculates CI between fracture and matrix. Inputs are
% matrix cell nodes (8x3 matrix), the area of intersection, fracture plane
% normal, fracture plane point, existing davg calculation and a numerical
% tolerance.
% If davg is negative, then the function reads this as davg has not been
% calculated before. However, if davg is positive, the function does not
% calculate davg again.

if davg<(-tol) % this is to help save computational time since davg calculation is quite cumbersome
    davg=calcdavg(nodes,planenormal,planepoint,tol);
end

CI=area/davg;
end