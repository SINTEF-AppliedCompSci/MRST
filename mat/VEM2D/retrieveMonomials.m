function [m, grad_m, int_m] = retrieve2DMonomials(k)
%--------------------------------------------------------------------------
%   Returns all 2D monomials of degree <= k, along with their gradients and
%   their anti-derivatives with respect to baricentric coordinates.
%
%   SYNOPSIS:
%       [m, grad_m, int_m] = retrieveMonomials(k)
%
%   DESCRIPTION:
%       Returns all 2D monomials of degree <= k, along with their
%       gradients, expressed in baricentric coordinates, and their
%       anti-derivatives with respect to the first baricentric coordinate,
%       also expressed in baricentric coordinates. Let xK and hK be the
%       centroid and diameter of cell K. Denote the monomials expressed in
%       cartesian coordinates by M. Then,
%
%           M(x)            = m((x-xK)/hK),
%           \nabla M(x)     = grad_m((x-xK)/hK)/hK,
%           \int M(x) \dx_1 = int_m((x-xK)/hK)*hK.
%
%       See [1] for detail.
%
%       The function handels are ordered in the following way:
%
%           m      = [m^(0,0)     , m^(1,0), m^(0,1), ...]
%           grad_m = [grad_m^(1,0); grad_m^(0,1)    ; ...]
%           int_m  = [int_m^(0,0) , int_m^(1,0)     , ...]   
%
%   REQUIRED PARAMETERS:
%       k       - Order of monomial space. Supported orders are 1 and 2.
%
%   RETURNS:
%       m       - Funciton handle of monomials.
%       grad_m  - Function handle of monomial gradients.
%       int_m   - Function handle of anti-derivative wrt first baricentric
%                 coordinate of monomials.
%
%   REFERENCES:
%       [1]     - Thesis title.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

assert(k == 1 | k == 2);

if k == 1  
    m      = @(X)                ...
             [ones(size(X,1),1), ...                        % (0,0)
              X(:,1)           , ...                        % (1,0)
              X(:,2)           ];                           % (0,1)

    grad_m = @(X) ...
             [ones(size(X,1),1) , zeros(size(X,1),1); ...   % (1,0)
              zeros(size(X,1),1), ones(size(X,1),1) ];      % (0,1)

    int_m  = @(X)              ...
             [X(:,1)         , ...                          % (0,0)
              X(:,1).^2/2    , ...                          % (1,0)
              X(:,1).*X(:,2)];                              % (0,1)
end

if k == 2 
    m      = @(X)                ...
             [ones(size(X,1),1), ...                        % (0,0)
              X(:,1)           , ...                        % (1,0)
              X(:,2)           , ...                        % (0,1)
              X(:,1).^2        , ...                        % (2,0)
              X(:,1).*X(:,2)   , ...                        % (1,1)
              X(:,2).^2        ];                           % (0,2)

    grad_m = @(X) ...
             [ones(size(X,1),1) , zeros(size(X,1),1); ...   % (1,0)
              zeros(size(X,1),1), ones(size(X,1),1) ; ...   % (0,1)
              2*X(:,1)          , zeros(size(X,1),1); ...   % (2,0)
              X(:,2)            , X(:,1)            ; ...   % (1,1)
              zeros(size(X,1),1), 2*X(:,2)          ];      % (0,2)

    int_m  = @(X)                  ...
             [X(:,1)             , ...                      % (0,0)
              X(:,1).^2/2        , ...                      % (1,0)
              X(:,1).*X(:,2)     , ...                      % (0,1)
              X(:,1).^3/3        , ...                      % (2,0)
              X(:,1).^2.*X(:,2)/2, ...                      % (1,1)
              X(:,1).*X(:,2).^2  ];                         % (0,2)
end