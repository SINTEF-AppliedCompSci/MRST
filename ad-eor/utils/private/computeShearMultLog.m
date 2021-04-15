function zSh = computeShearMultLog(fluid, vW, muWMultf)
    % The current version handles all the values one by one
    % It needs to be improved for better performance.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    refConcentration = fluid.plyshlog.refcondition(1);
    refViscMult = fluid.muWMult(refConcentration);

    plyshlogTable = fluid.plyshlog.data{1};

    % the water velocity/shear rate in the PLYSHLOG table
    waterVel = plyshlogTable(:,1);
    % the viscosity factors in the PLYSHLOG table
    refM = plyshlogTable(:,2);

    % converting the table using the reference condition
    refM = (refViscMult * refM -1)/(refViscMult-1);

    % the minimum velocity/shear rate value specified in the table
    % it should be the first entry value
    minWaterVel = min(waterVel);

    % only calcuate shear factors when polymer exists and velocity/shear rate
    % is big enough
    iShear = find((vW > minWaterVel) & (muWMultf > 1.));

    P = muWMultf(iShear);

    % convert the velocity/shear rate to their logarithms
    V = log(waterVel);
    vW0 = log(vW(iShear));

    % function to decide the relative location of a point to a line
    f = @(x,y, x0) x+y -x0;

    % the shear factors remain to be calculated
    z = ones(numel(iShear),1);

    for i = 1:numel(iShear)

       % for each value from P, we need to generate a table to calculate
       % the shear factor
       Z = (1 + (P(i) - 1) * refM)/ P(i);
       % covert the shear factors to their logarithms
       Z = log(Z);

       % to flag the relative location of the line and line segments
       % to help to determine the intersection point faster
       sign = f(V, Z, vW0(i));

       % to check if there is sign change
       temp = sign(1:end-1).*sign(2:end);

       % to find the index that sign changed
       j = find(temp <= 0);

       % make sure one and only one intersection point is found.
       assert(numel(j) <= 1);

       if (numel(j) == 1)
           [~,z(i)] = findIntersection([V(j), Z(j); V(j+1), Z(j+1)], [0, vW0(i); vW0(i), 0]);
       end

       if (numel(j) == 0)
           % out of the table range
           % since we handled the small value already, this must be a big
           % value.
           assert(vW0(i) - Z(end) > V(end));
           z(i) = Z(end);
       end

    end

    z = exp(z);

    % the shear factors
    zSh = ones(numel(vW), 1);
    zSh(iShear) = z;

end

%--------------------------------------------------------------------------

% finding the intersection of one line segment l1 and one straight line l2
% l1 is defined with the beginning and ending points.
% l2 is defined with two points along the line.
function [x, y] = findIntersection(l1, l2)

    assert(all(size(l1) == [2 2]));
    assert(all(size(l2) == [2 2]));

    x1 = l1(1,1);
    y1 = l1(1,2);
    x2 = l1(2,1);
    y2 = l1(2,2);

    x3 = l2(1,1);
    y3 = l2(1,2);
    x4 = l2(2,1);
    y4 = l2(2,2);

    d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);

    assert(d ~= 0);

    x = ((x3-x4)*(x1*y2-y1*x2)-(x1-x2)*(x3*y4-y3*x4))/d;
    y = ((y3-y4)*(x1*y2-y1*x2)-(y1-y2)*(x3*y4-y3*x4))/d;

    assert(x >= min(x1,x2) && x <= max(x1,x2));
end
