function [Xq, w, V, vol] = triangleQuadRule(k)
%   Returns quadrature rule for the reference triangle with vertices (0,0),
%   (1,0) and (0,1).
%
%   SYNOPSIS:
%       [Xq, w, V, vol] = triangleQuadRule(k)
%
%   DESCRIPTION:
%       Returns quadrature rule of precission k. k = 1 returns the centroid
%       rule, while k >= 2 returns rules named STRANG k, as described in
%       [1]. Usage of the rule is as follows:
%
%           \int_T f \dx = vol\sum_{i = 1}^n w_i*f(Xq_i),
%
%       where T is the reference triangle with vertices (0,0), (1,0) and
%       (0,1), f is the funtion to be integrated, vol is the area of T, and
%       w_i and Xq_i is the ith wheight and quadrature point, respectively.
%
%   REQUIRED PARAMETERS:
%       k       - Quadrature rule precision. supported precisions
%                 are 1,2,3 and 7.
%
%   RETURNS:
%       Xq      - nq x 2 matrix of quadrature points.
%       w       - nq x 1 vector of quadrature wheights.
%       V       -  3 x 2 matrix of reference tirangle vertices.
%       vol     - Area of reference triangel.
% 
%   REFERENCES:
%       [1]     - http://people.sc.fsu.edu/~jburkardt/datasets/...
%                        quadrature_rules_tri/quadrature_rules_tri.html
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

assert(k >= 1 & (k <= 3 | k == 7), ...
      'Only supported quadrature rules are of precision 1, 2, 3 and 7');

V = [ 0.0, 0.0;
      1.0, 0.0;
      0.0, 1.0 ];

vol = 1/2;

if k == 1
    
    Xq = [  0.33333333333333333333  0.33333333333333333333 ];
    
    w = 1;

elseif k == 2

    Xq = [ ...
           0.50000000000000000000  0.00000000000000000000;
           0.50000000000000000000  0.50000000000000000000;
           0.00000000000000000000  0.50000000000000000000 ];

     w = [ ...
           0.33333333333333333333, ...
           0.33333333333333333333, ...
           0.33333333333333333333     ];

elseif k == 3
    
    Xq = [ ...
           0.33333333333333333333  0.33333333333333333333;
           0.60000000000000000000  0.20000000000000000000;
           0.20000000000000000000  0.60000000000000000000;
           0.20000000000000000000  0.20000000000000000000 ];
   
    w =  [ ...
          -0.56250000000000000000, ...
           0.52083333333333333333, ...
           0.52083333333333333333, ...
           0.52083333333333333333     ];
       
elseif k == 7
    
    Xq = [ ...
           0.333333333333333  0.333333333333333;
           0.479308067841923  0.260345966079038;
           0.260345966079038  0.479308067841923;
           0.260345966079038  0.260345966079038;
           0.869739794195568  0.065130102902216;
           0.065130102902216  0.869739794195568;
           0.065130102902216  0.065130102902216;
           0.638444188569809  0.312865496004875;
           0.638444188569809  0.048690315425316;
           0.312865496004875  0.638444188569809;
           0.312865496004875  0.048690315425316;
           0.048690315425316  0.638444188569809;
           0.048690315425316  0.312865496004875];
       
     w = [ ...
          -0.149570044467670, ...
           0.175615257433204, ...
           0.175615257433204, ...
           0.175615257433204, ...
           0.053347235608839, ...
           0.053347235608839, ...
           0.053347235608839, ...
           0.077113760890257, ...
           0.077113760890257, ...
           0.077113760890257, ...
           0.077113760890257, ...
           0.077113760890257, ...
           0.077113760890257];
       
end