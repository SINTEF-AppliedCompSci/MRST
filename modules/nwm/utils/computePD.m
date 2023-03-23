function f = computePD(x, y, a, b, xw, yw)
% Compute dimensionless pressure for the flow problem that a well is 
% arbitrarily located inside a rectangular box with width a and height b. 
% The distance from well to right boundary is xw and to the lower boundary 
% is yw
%       ------------------------------------------
%      |                               xw         |
%      |                 (Well) o ............... |
%      |                        .                 |
%      |                    yw  .                 |
%      |                        .                 |
%      |                        .                 |
%       -------------------------------------------

    % A constant for computing the infinite series
    N = 10; 
    nn = (-N : N)';
    p1 = pi/b * (x - 2*nn*a);
    p2 = pi/b * (x - 2*nn*a - 2*xw);

    q1 = pi/b *  y;
    q2 = pi/b * (y + 2*yw);
    s1 = cosh(p1) - cos(q2);
    s2 = cosh(p2) - cos(q1);
    s3 = cosh(p1) - cos(q1);
    s4 = cosh(p2) - cos(q2);

    f = log( (s1.*s2)./(s3.*s4) );
    f = sum(f);
end