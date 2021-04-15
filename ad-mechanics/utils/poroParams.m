function [res, rel_hist, residuals, eqs] = poroParams(phi, uniform, varargin)
%
% From a set of initially known poroelastic parameters, compute as many as
% possible of the rest. 
%
% The full list of supported parameters are:
% 
% 'K'         - drained bulk modulus
% 'H'         - inverse of poroelastic expansion coefficient
% 'R'         - inverse of S_gamma
% 'M'         - inverse of S_epsilon
% 'K_s'       - unjacketed bulk modulus
% 'K_p'       - inverse of drained pore compressibility
% 'K_f'       - inverse of fluid compressibility
% 'K_phi'     - inverse of unjacketed pore compressibility
% 'K_v'       - uniaxial drained bulk modulus
% 'S'         - uniaxial specific storage coefficient
% 'S_sigma'   - unconstrained specific storage coefficient
% 'S_epsilon' - constrained specific storage coefficient
% 'S_gamma'   - unjacketed specific storage
% 'alpha'     - Biot-Willis coef. (eff. stress coef. for bulk volume)
% 'beta'      - eff. stress coef. for pore volume
% 'gamma'     - loading efficency
% 'eta'       - poroelastic stress coefficient
% 'c_m'       - Geertsma's parameter
% 'B'         - Skempton's coefficient
% 'K_u'       - undrained bulk modulus
% 'K_vu'      - uniaxial undrained bulk modulus
% 'E'         - Young's modulus (drained)
% 'E_u'       - Young's modulus (undrained)
% 'lambda'    - Lame's parameter (drained)
% 'lambda_u'  - Lame's parameter (undrained)
% 'nu'        - Poisson's ratio (drained)
% 'nu_u'      - Poisson's ratio (undrained)
% 'G'         - Shear modulus (drained and undrained)
%   
% If the equation solver was not able to satisfy the full set of poroelastic
% relations, a warning is shown.  Usually, this is because the user has
% over-specified the problem with incompatible parameter values.
% 
% SYNOPSIS:
%   function [res, rel_hist, residuals, eqs] = poroParams(phi, uniform, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   phi      - porosity of the medium (must be provided)
%   uniform - 'true' if the porous medium is assumed to have uniform grain
%              compressibility (in which case the unjacketed bulk modulus equals
%              the solid-grain modulus)
%   varargin - list on the form: 'param', value
%              where 'param' is the name of one of the poroelastic parameters
%              listed above, and 'value' is its respective value               
%
% RETURNS:
%   res       - structure containing all the supported poroelastic
%               parameters.  Those that could not be computed from the input
%               provided are left as 'NaN's.
%   rel_hist  - Lists the ordered sequence of relations and poroelastic
%               parameters computed during the execution of the algorithm
%   residuals - the residuals associated with each of the poroelastic
%               relationships available to the algorithm to compute the full
%               set of parameter values
%   eqs       - a cell array that lists (in order) all the equations that
%               were used during execution of the algorithm in order to
%               compute the poroelastic parameters
%
% EXAMPLE:
%   poroParams(0.25, false, 'K', 1e9, 'H', 1.1e9, 'R', 1.2e9, 'G', 1.1e9) 

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

   % Specifying names of poroelastic parameters and capturing those provided
   % by user
   pnames = parameter_names();
   nans = cell(size(pnames)); [nans{:}] = deal(nan);
   params = merge_options(cell2struct(nans, pnames, 2), varargin{:});
   
   params = rescale_params(params, 'reverse', false); % rescale for numerical robustness
   
   % Specifying relations
   relations = setup_relations(phi);
   
   % If uniform grains, we know that K_phi = K_s, so we add this relation
   if uniform
      relations = [relations;
                   relfun({'K_s', 'K_phi'},   @(K_s, K_phi) 1/K_s - 1/K_phi)];
   end
   
   % Compute all possible parameters based on data provided
   rel_hist = {};
   while true
   
      % Identify all relations that can be solved with what we have 
      [rel_ixs, found_pnames] = find_relations_with_single_unknown(relations, params);
         
      % Removing multiple relations for computing the same variable
      [~, ix] = unique(found_pnames);
      found_pnames = found_pnames(ix);
      rel_ixs = rel_ixs(ix);
      
      if isempty(rel_ixs)
         % no more relations with a single unknown found
         break;
      end
      
      % Solving the identified relations
      for i = 1:numel(rel_ixs)
         rel_ix = rel_ixs(i);
         pname = found_pnames{i};
         params.(pname) = solve_relation(relations(rel_ix), pname, params);

         % Keeping track of history
         rel_hist = [rel_hist, {rel_ix, pname}]; %#ok
      end
   end
   
   residuals = compute_residuals(relations, params);
   if any(abs(residuals) > 1e-3) % @@ This should really be scaled with the
                                 % maginitude of the involved parameters
      warning(['Some residuals did not vanish.  Input parameters might overspecify ' ...
               'system.']); 
   end
   
   params = rescale_params(params, 'reverse', true); % rescale for numerical robustness
   
   eqs       = relations;
   res       = params;
   
end

% ----------------------------------------------------------------------------
function params = rescale_params(params, varargin)
   opt.reverse = false;
   opt = merge_options(opt, varargin{:});
   
   modfac = 1/1e7; % scaling factor for moduli
   comprfac = 1/modfac; % scaling factor for compressibilities

   if opt.reverse
      modfac = 1/modfac;
      comprfac = 1/comprfac;
   end
   
   modnames = {'K', 'H', 'R', 'M', 'K_s', 'K_p', 'K_f', 'K_phi', 'K_v', 'K_u', ...
               'K_vu', 'E', 'E_u', 'lambda', 'lambda_u', 'G'};
   comprnames = {'S', 'S_sigma', 'S_epsilon', 'S_gamma', 'c_m'};
   
   for m = modnames
      params.(m{:}) = params.(m{:}) * modfac;
   end

   for m = comprnames
      params.(m{:}) = params.(m{:}) * comprfac;
   end
   
end

% ----------------------------------------------------------------------------

function res = compute_residuals(rels, params)
   res = nan(numel(rels),1);
   
   for i = 1:numel(rels)
      res(i) = feval2(rels(i).f, extract_vals_from_struct(rels(i).args, params));
   end
end

% ----------------------------------------------------------------------------
function carray = change_entry(carray, ix, val)
   carray{ix} = val;
end

% ----------------------------------------------------------------------------

function val = feval2(f, args)
% wrapper to 'feval', allowing arguments to be passed as a single cell array
   val = feval(f, args{:});
end

% ----------------------------------------------------------------------------

function val = solve_relation(rel, pname, params)

   if mrstVerbose 
      fprintf('Now solving for: %s\n', pname);
   end
   
   vals = extract_vals_from_struct(rel.args, params);
   p_ix = find(isnan([vals{:}]));
   assert(numel(p_ix) == 1); % this function should only be called when exactly one
                             % parameter is missing
   
   g = @(x) feval2(rel.f, change_entry(vals, p_ix, x));
   
   % Search for zero.  Naive implementation as for now @@
   %init_guesses = 4.^(-1:20);
   init_guesses = 10.^(-1:14);
   init_guesses = [init_guesses;-init_guesses];
   init_guesses = init_guesses(:)';
   %init_guesses = [-1 0 init_guesses];
   
   for guess = init_guesses
      status = 0;
      try
         %[val, ~, status] = fzero(g, guess);
         [val, ~, status] = fzero(g, guess, struct('TolX', 1e-40));
      catch
         % do nothing, keep guessing
      end
      if status == 1
         % minimum found
         return;
      end
   end
   %Minimum not found.  Check if an infinite value gives a zero
   if g(Inf) < 1e-12
      val = Inf;
      return
   end
   error(['Could not identify minimum when computing ' pname '.']);
end

% ----------------------------------------------------------------------------

function [rel_ixs, pnames] = find_relations_with_single_unknown(rels, params)

   rel_ixs = [];
   pnames = {};
   
   % searching for a relation where all parameters except one are known
   for i = 1:numel(rels)
      vals = extract_vals_from_struct(rels(i).args, params);
      
      if sum(isnan([vals{:}])) == 1
         % We found a relation with exactly one unknown
         rel_ixs = [rel_ixs, i]; %#ok
         p_ix = find(isnan([vals{:}]));
         pname = rels(i).args(p_ix);
         pnames = {pnames{:}, pname{:}}; %#ok
         %return;
      end
   end
end

% ----------------------------------------------------------------------------

function res = extract_vals_from_struct(fields, str)

   res = cell(numel(fields),1);
   for i = 1:numel(res)
      res{i} = str.(fields{i});
   end
end

% ----------------------------------------------------------------------------
function pnames = parameter_names()
   
   % Fundamental poroelastic constant
   pnames = {'K'         , ...  % drained bulk modulus
             'H'         , ...  % inverse of poroelastic expansion coefficient
             'R'         , ...  % inverse of S_gamma
             'M'};              % inverse of S_epsilon
   
   % Other compressibilities
   pnames = [pnames      , ... 
            {'K_s'       , ... % unjacketed bulk modulus
             'K_p'       , ... % inverse of drained pore compressibility
             'K_f'       , ... % inverse of fluid compressibility
             'K_phi'     , ... % inverse of unjacketed pore compressibility
             'K_v'}];           % uniaxial drained bulk modulus
   
   % Storativities
   pnames = [pnames      , ...
            {'S'         , ... % uniaxial specific storage coefficient
             'S_sigma'   , ... % unconstrained specific storage coefficient
             'S_epsilon' , ... % constrained specific storage coefficient
             'S_gamma'}];       % unjacketed specific storage

   % Special parameters
   pnames = [pnames      , ...
            {'alpha'     , ... % Biot-Willis coef. (eff. stress coef. for bulk volume)
             'beta'      , ... % eff. stress coef. for pore volume
             'gamma'     , ... % loading efficency
             'eta'       , ... % poroelastic stress coefficient
             'c_m'       , ... % Geertsma's parameter
             'B'}];             % Skempton's coefficient
   
   % Other parameters from linear elasticity, drained and undrained             
   pnames = [pnames      , ...
            {'K_u'       , ... % undrained bulk modulus
             'K_vu'      , ... % uniaxial undrained bulk modulus
             'E'         , ... % Young's modulus (drained)
             'E_u'       , ... % Young's modulus (undrained)
             'lambda'    , ... % Lame's parameter (drained)
             'lambda_u'  , ... % Lame's parameter (undrained)
             'nu'        , ... % Poisson's ratio (drained)
             'nu_u'      , ... % Poisson's ratio (undrained)
             'G'}];             % Shear modulus (drained and undrained)
end

% ----------------------------------------------------------------------------

function f = relfun(carray, fun)
   f.args = carray;
   f.f = fun;
end


% ----------------------------------------------------------------------------

function rel = setup_relations(phi)

% The relations are all taken from:
% Wang, H.F., "Theory of linear poroelasticity", Princeton Series in
% Geophysics, Princeton University Press, NJ (2000).
   
   %relfun = @(carray, fun) struct('args', {carray}, 'f', fun);
   
   % Standard relations in linear elasticity (drained)
   rel = [relfun({'K',   'lambda', 'G'},  @(K,   lambda, G) 1/K - 1/(lambda + 2*G/3));                ...
          relfun({'E',   'lambda', 'G'},  @(E,   lambda, G) E - (G*(3*lambda + 2 * G)/(lambda + G))); ...
          relfun({'nu',  'lambda', 'G'},  @(nu,  lambda, G) nu - (lambda/(2*(lambda+G))));            ...
          relfun({'K_v', 'lambda', 'G'},  @(K_v, lambda, G) 1/K_v - 1/(lambda + 2*G));                ...
          relfun({'E', 'K', 'lambda'},    @(E, K, lambda) E - (9*K * (K-lambda)/(3*K-lambda)));       ...
          relfun({'nu', 'K', 'lambda'},   @(nu,K, lambda) nu - (lambda/(3*K-lambda)));                ...
          relfun({'K_v', 'K', 'lambda'},  @(K_v, K, lambda) 1/K_v - 1/(3*K-2*lambda));                ...
          relfun({'G', 'K', 'lambda'},    @(G, K, lambda) G - (3/2*(K-lambda)));                      ...
          relfun({'E', 'K', 'G'},         @(E, K, G) E - (9*K*G/(3*K+G)));                            ...
          relfun({'lambda', 'K', 'G'},    @(lambda, K, G) lambda - (K-2*G/3));                        ...
          relfun({'nu', 'K', 'G'},        @(nu, K, G) nu - ((3*K-2*G)/(2*(3*K+G))));                  ...
          relfun({'K_v', 'K', 'G'},       @(K_v, K, G) 1/K_v - 1/(K+4*G/3));                          ...
          relfun({'K', 'E', 'G'},         @(K, E, G) 1/K - 1/(E*G/(3*(3*G-E))));                      ...
          relfun({'lambda', 'E', 'G'},    @(lambda, E, G) lambda - (G*(E-2*G)/(3*G-E)));              ...
          relfun({'nu', 'E', 'G'},        @(nu, E, G) nu - (E/(2*G) - 1));                            ...
          relfun({'K_v', 'E', 'G'},       @(K_v, E, G) 1/K_v - 1/(G*(4*G-E)/(3*G-E)));                ...
          relfun({'lambda', 'K', 'E'},    @(lambda, K, E) lambda - (3*K*(3*K-E)/(9*K-E)));            ... 
          relfun({'nu', 'K', 'E'},        @(nu, K, E) nu - ((3*K-E)/(6*K)));                          ...
          relfun({'K_v', 'K', 'E'},       @(K_v, K, E) 1/K_v - 1/(3*K*(3*K+E)/(9*K-E)));              ...
          relfun({'G', 'K', 'E'},         @(G, K, E) G - (3*K*E/(9*K-E)));                            ...
          relfun({'K', 'lambda', 'nu'},   @(K, lambda, nu) 1/K - 1/(lambda * (1+nu)/(3*nu)));         ...
          relfun({'E', 'lambda', 'nu'},   @(E, lambda, nu) E - (lambda * ((1+nu)*(1-2*nu))/(nu)));    ...
          relfun({'K_v', 'lambda', 'nu'}, @(K_v, lambda, nu) 1/K_v - 1/(lambda * (1-nu)/nu));         ...
          relfun({'G', 'lambda', 'nu'},   @(G, lambda, nu) G - (lambda * (1-2*nu)/(2*nu)));           ...
          relfun({'K', 'G', 'nu'},        @(K, G, nu) 1/K - 1/(G * (2*(1+nu))/(3*(1-2*nu))));         ...
          relfun({'E', 'G', 'nu'},        @(E, G, nu) E - (2*G*(1+nu)));                              ...
          relfun({'lambda', 'G', 'nu'},   @(lambda, G, nu) lambda - (G * 2 * nu / (1-2*nu)));         ...
          relfun({'K_v', 'G', 'nu'},      @(K_v, G, nu) 1/K_v - 1/(G * (2-2*nu)/(1-2*nu)));           ...
          relfun({'E', 'K', 'nu'},        @(E, K, nu) E - (3*K*(1-2*nu)));                            ...
          relfun({'lambda', 'K', 'nu'},   @(lambda, K, nu) lambda - (3 * K * nu / (1+nu)));           ...
          relfun({'K_v', 'K', 'nu'},      @(K_v, K, nu) 1/K_v - 1/(3*K * (1-nu)/(1+nu)));             ...
          relfun({'G', 'K', 'nu'},        @(G, K, nu) G - (3*K*(1-2*nu)/(2+2*nu)));                   ...
          relfun({'K', 'E', 'nu'},        @(K, E, nu) 1/K - 1/(E/(3*(1-2*nu))));                      ...
          relfun({'lambda', 'E', 'nu'},   @(lambda, E, nu) lambda - ((E*nu)/((1+nu)*(1-2*nu))));      ...
          relfun({'K_v', 'E', 'nu'},      @(K_v, E, nu) 1/K_v - 1/((E*(1-nu))/((1+nu)*(1-2*nu))));    ...
          relfun({'G', 'E', 'nu'},        @(G, E, nu) G - (E/(2+2*nu)));                              ...
          ];

   % Standard relations in linear elasticity (undrained)
   rel = [rel; ...
          relfun({'K_u',   'lambda_u', 'G'},   @(K_u,   lambda_u, G) 1/K_u - 1/(lambda_u + 2*G/3));                     ...
          relfun({'E_u',   'lambda_u', 'G'},   @(E_u,   lambda_u, G) E_u - (G*(3*lambda_u + 2 * G)/(lambda_u + G)));    ...
          relfun({'nu_u',  'lambda_u', 'G'},   @(nu_u,  lambda_u, G) nu_u - (lambda_u/(2*(lambda_u+G))));               ...
          relfun({'K_vu', 'lambda_u', 'G'},    @(K_vu, lambda_u, G) 1/K_vu - 1/(lambda_u + 2*G));                       ...
          relfun({'E_u', 'K_u', 'lambda_u'},   @(E_u, K_u, lambda_u) E_u - (9*K_u * (K_u-lambda_u)/(3*K_u-lambda_u)));  ...
          relfun({'nu_u', 'K_u', 'lambda_u'},  @(nu_u,K_u, lambda_u) nu_u - (lambda_u/(3*K_u-lambda_u)));               ...
          relfun({'K_vu', 'K_u', 'lambda_u'},  @(K_vu, K_u, lambda_u) 1/K_vu - 1/(3*K_u-2*lambda_u));                   ...
          relfun({'G', 'K_u', 'lambda_u'},     @(G, K_u, lambda_u) G - (3/2*(K_u-lambda_u)));                           ...
          relfun({'E_u', 'K_u', 'G'},          @(E_u, K_u, G) E_u - (9*K_u*G/(3*K_u+G)));                               ...
          relfun({'lambda_u', 'K_u', 'G'},     @(lambda_u, K_u, G) lambda_u - (K_u-2*G/3));                             ...
          relfun({'nu_u', 'K_u', 'G'},         @(nu_u, K_u, G) nu_u - ((3*K_u-2*G)/(2*(3*K_u+G))));                     ...
          relfun({'K_vu', 'K_u', 'G'},         @(K_vu, K_u, G) 1/K_vu - 1/(K_u+4*G/3));                                 ...
          relfun({'K_u', 'E_u', 'G'},          @(K_u, E_u, G) 1/K_u - 1/(E_u*G/(3*(3*G-E_u))));                         ...
          relfun({'lambda_u', 'E_u', 'G'},     @(lambda_u, E_u, G) lambda_u - (G*(E_u-2*G)/(3*G-E_u)));                 ...
          relfun({'nu_u', 'E_u', 'G'},         @(nu_u, E_u, G) nu_u - (E_u/(2*G) - 1));                                 ...
          relfun({'K_vu', 'E_u', 'G'},         @(K_vu, E_u, G) 1/K_vu - 1/(G*(4*G-E_u)/(3*G-E_u)));                     ...
          relfun({'lambda_u', 'K_u', 'E_u'},   @(lambda_u, K_u, E_u) lambda_u - (3*K_u*(3*K_u-E_u)/(9*K_u-E_u)));       ... 
          relfun({'nu_u', 'K_u', 'E_u'},       @(nu_u, K_u, E_u) nu_u - ((3*K_u-E_u)/(6*K_u)));                         ...
          relfun({'K_vu', 'K_u', 'E_u'},       @(K_vu, K_u, E_u) 1/K_vu - 1/(3*K_u*(3*K_u+E_u)/(9*K_u-E_u)));           ...
          relfun({'G', 'K_u', 'E_u'},          @(G, K_u, E_u) G - (3*K_u*E_u/(9*K_u-E_u)));                             ...
          relfun({'K_u', 'lambda_u', 'nu_u'},  @(K_u, lambda_u, nu_u) 1/K_u - 1/(lambda_u * (1+nu_u)/(3*nu_u)));        ...
          relfun({'E_u', 'lambda_u', 'nu_u'},  @(E_u, lambda_u, nu_u) E_u - (lambda_u * ((1+nu_u)*(1-2*nu_u))/(nu_u))); ...
          relfun({'K_vu', 'lambda_u', 'nu_u'}, @(K_vu, lambda_u, nu_u) 1/K_vu - 1/(lambda_u * (1-nu_u)/nu_u));          ...
          relfun({'G', 'lambda_u', 'nu_u'},    @(G, lambda_u, nu_u) G - (lambda_u * (1-2*nu_u)/(2*nu_u)));              ...
          relfun({'K_u', 'G', 'nu_u'},         @(K_u, G, nu_u) 1/K_u - 1/(G * (2*(1+nu_u))/(3*(1-2*nu_u))));            ...
          relfun({'E_u', 'G', 'nu_u'},         @(E_u, G, nu_u) E_u - (2*G*(1+nu_u)));                                   ...
          relfun({'lambda_u', 'G', 'nu_u'},    @(lambda_u, G, nu_u) lambda_u - (G * 2 * nu_u / (1-2*nu_u)));            ...
          relfun({'K_vu', 'G', 'nu_u'},        @(K_vu, G, nu_u) 1/K_vu - 1/(G * (2-2*nu_u)/(1-2*nu_u)));                ...
          relfun({'E_u', 'K_u', 'nu_u'},       @(E_u, K_u, nu_u) E_u - (3*K_u*(1-2*nu_u)));                             ...
          relfun({'lambda_u', 'K_u', 'nu_u'},  @(lambda_u, K_u, nu_u) lambda_u - (3 * K_u * nu_u / (1+nu_u)));          ...
          relfun({'K_vu', 'K_u', 'nu_u'},      @(K_vu, K_u, nu_u) 1/K_vu - 1/(3*K_u * (1-nu_u)/(1+nu_u)));              ...
          relfun({'G', 'K_u', 'nu_u'},         @(G, K_u, nu_u) G - (3*K_u*(1-2*nu_u)/(2+2*nu_u)));                      ...
          relfun({'K_u', 'E_u', 'nu_u'},       @(K_u, E_u, nu_u) 1/K_u - 1/(E_u/(3*(1-2*nu_u))));                       ...
          relfun({'lambda_u', 'E_u', 'nu_u'},  @(lambda_u, E_u, nu_u) lambda_u - ((E_u*nu_u)/((1+nu_u)*(1-2*nu_u))));   ...
          relfun({'K_vu', 'E_u', 'nu_u'},      @(K_vu, E_u, nu_u) 1/K_vu - 1/((E_u*(1-nu_u))/((1+nu_u)*(1-2*nu_u))));   ...
          relfun({'G', 'E_u', 'nu_u'},         @(G, E_u, nu_u) G - (E_u/(2+2*nu_u)));                                   ...
         ];
   
   % other relations
   rel = [rel; ...
          relfun({'S_sigma', 'R'},                            @(S_sigma, R) S_sigma - (1/R));                                                              ...
          relfun({'B', 'R', 'H'},                             @(B, R, H) B - (R/H));                                                                       ...
          relfun({'S_epsilon', 'M'},                          @(S_epsilon, M) S_epsilon - (1/M));                                                          ...
          relfun({'S_epsilon', 'S_sigma', 'K', 'H'},          @(S_epsilon, S_sigma, K, H) S_epsilon - (S_sigma - K/(H^2)));                                ...
          relfun({'alpha', 'K', 'H'},                         @(alpha, K, H) alpha - K/H);                                                                 ...
          relfun({'S_sigma', 'alpha', 'K', 'B'},              @(S_sigma, alpha, K, B) S_sigma - (alpha/(K*B)));                                            ...
          relfun({'K_u', 'K', 'alpha', 'B'},                  @(K_u, K, alpha, B) 1/K_u - ((1-alpha*B)/K));                                                ...
          relfun({'S_epsilon', 'alpha', 'K_u', 'B'},          @(S_epsilon, alpha, K_u, B) S_epsilon - (alpha/(K_u*B)));                                    ...
          relfun({'K_s', 'K', 'H'},                           @(K_s, K, H) 1/K_s - (1/K - 1/H));                                                           ...
          relfun({'K_s', 'K', 'alpha'},                       @(K_s, K, alpha) 1/K_s - (1/K - alpha/K));                                                   ...
          relfun({'K_p', 'alpha', 'K'},                       @(K_p, alpha, K) 1/K_p - (alpha/(phi*K)));                                                   ...
          relfun({'K_phi', 'alpha', 'K', 'B', 'K_f'},         @(K_phi, alpha, K, B, K_f) 1/K_phi + 1/phi*(alpha/(K*B) - phi/K_f - alpha/K));               ...
          relfun({'K_phi', 'K_s', 'K', 'B', 'K_f'},           @(K_phi, K_s, K, B, K_f) 1/K_phi + 1/phi*((1/K - 1/K_s)*(1/B - 1) - phi/K_f));               ...
          relfun({'beta', 'K_p', 'K_phi'},                    @(beta, K_p, K_phi) beta - (1-K_p/K_phi));                                                   ...
          relfun({'nu_u', 'nu', 'alpha', 'B'},                @(nu_u, nu, alpha, B) nu_u - ((3*nu + alpha*B*(1-2*nu))/(3-alpha*B*(1-2*nu))));              ...
          relfun({'S_sigma', 'K', 'K_s', 'K_f', 'K_phi'},     @(S_sigma, K, K_s, K_f, K_phi) S_sigma-((1/K-1/K_s) + phi*(1/K_f - 1/K_phi)));               ...
          relfun({'S_epsilon', 'S_sigma', 'alpha', 'K'},      @(S_epsilon, S_sigma, alpha, K) S_epsilon - (S_sigma - alpha^2/K));                          ...
          relfun({'S_epsilon', 'S_sigma', 'K', 'K_u'},        @(S_epsilon, S_sigma, K, K_u) S_epsilon/S_sigma - K/K_u);                                    ...
          relfun({'S_epsilon', 'S_sigma', 'alpha', 'B'},      @(S_epsilon, S_sigma, alpha, B) S_epsilon - (S_sigma * (1-alpha*B)));                        ...
          relfun({'S_epsilon', 'alpha', 'K_u', 'K'},          @(S_epsilon, alpha, K_u, K) S_epsilon - alpha^2/(K_u-K));                                    ...
          relfun({'S_epsilon', 'K', 'K_s', 'K_f', 'K_phi'},   @(S_epsilon, K, K_s, K_f, K_phi) S_epsilon - (1/K_s * (1-K/K_s) + phi * (1/K_f - 1/K_phi))); ...
          relfun({'M', 'G', 'nu', 'nu_u', 'alpha'},           @(M, G, nu, nu_u, alpha) M - (2*G*(nu_u-nu)/(alpha^2*(1-2*nu_u)*(1-2*nu))));                 ...
          relfun({'eta', 'nu', 'alpha'},                      @(eta, nu, alpha) eta - ((1-2*nu)*alpha/(2*(1-nu))));                                        ...
          relfun({'S', 'S_sigma', 'eta', 'B'},                @(S, S_sigma, eta, B) S - (S_sigma * (1 - 4*eta*B/3)));                                      ...
   % relfun({'S', 'alpha', 'K_v', 'K_f', 'K_phi'},       @(S, alpha, K_v, K_f, K_phi) S - ((alpha^2)/K_v + phi*(1/K_f - 1/K_phi))); % @@ apparent error in formula
          relfun({'S', 'S_epsilon', 'nu', 'nu_u'},            @(S, S_epsilon, nu, nu_u) S - (S_epsilon * ((1-nu_u)/(1-nu)) * ((1-2*nu)/(1-2*nu_u))));      ...
          relfun({'S', 'S_epsilon', 'K_vu', 'K_v'},           @(S, S_epsilon, K_vu, K_v) S - S_epsilon*K_vu/K_v);                                          ...
          relfun({'K_vu', 'K_u', 'nu_u'},                     @(K_vu, K_u, nu_u) 1/K_vu - ((1/(3*K_u)) * ((1+nu_u)/(1-nu_u)) ));                           ...
          relfun({'S', 'S_sigma', 'nu', 'nu_u'},              @(S, S_sigma, nu, nu_u) S - (S_sigma * ((1-nu_u)/(1+nu_u)) * ((1+nu)/(1-nu))));              ...
          relfun({'S', 'S_sigma', 'K', 'K_u', 'K_v', 'K_vu'}, @(S, S_sigma, K, K_u, K_v, K_vu) S - (S_sigma * (K/K_u) * (K_vu/K_v)));                      ...
          relfun({'S', 'alpha', 'K_v', 'gamma'},              @(S, alpha, K_v, gamma) S - alpha/(K_v*gamma));                                              ...
          relfun({'gamma', 'B', 'nu_u'},                      @(gamma, B, nu_u) gamma - (B*(1+nu_u))/(3*(1-nu_u)));                                        ...
          relfun({'gamma', 'eta', 'G', 'S'},                  @(gamma, eta, G, S) gamma - eta/(G*S));                                                      ...
          relfun({'S_gamma', 'K_f', 'K_phi'},                 @(S_gamma, K_f, K_phi) S_gamma - phi*(1/K_f - 1/K_phi));                                     ...
          relfun({'S_gamma', 'S_epsilon', 'alpha', 'K'},      @(S_gamma, S_epsilon, alpha, K) S_gamma - (S_epsilon - alpha*(1-alpha)/K));                  ...
          relfun({'alpha', 'K', 'K_s'},                       @(alpha, K, K_s) alpha - (1-K/K_s));                                                         ...
          relfun({'H', 'K_p'},                                @(H, K_p) 1/H - phi/K_p);                                                                    ...
          relfun({'c_m', 'alpha', 'K_v'},                     @(c_m, alpha, K_v) c_m - alpha/K_v);                                                         ...
          relfun({'c_m', 'eta', 'G'},                         @(c_m, eta, G) c_m - eta/G);                                                                 ...
          relfun({'B', 'K', 'K_s', 'K_f', 'K_phi'},           @(B, K, K_s, K_f, K_phi) B - ((1/K - 1/K_s)/(1/K - 1/K_s + phi*(1/K_f - 1/K_phi))));         ...
          relfun({'B', 'K', 'K_u', 'K_s'},                    @(B, K, K_u, K_s) B - ((1-K/K_u)/(1-K/K_s)));                                                ...
          relfun({'B', 'nu', 'nu_u', 'alpha'},                @(B, nu, nu_u, alpha) B - ((3*(nu_u-nu))/(alpha*(1+nu_u)*(1-2*nu))));                        ...
          ];
   
end


