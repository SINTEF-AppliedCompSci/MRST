function [HT_cg, T_cg, cgwells, report] = ...
      upscaleTransNew(cg, T_fine, varargin)
% Calculate upscaled transmissibilities for a coarse model
%
% SYNOPSIS:
%   [HT_cg, T_cg, cgwells, upscaled] = upscaleTrans(cg, T_fine)
%   [HT_cg, T_cg, cgwells, upscaled] = ...
%       upscaleTrans(cg, T_fine, 'pn1', pv1, ...)
%
% PARAMETERS:
%   cg       - coarse grid using the new-coarsegrid module structure
%
%   T_fine   - transmissibility on fine grid
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%               - Verbose -- Whether or not to display informational
%                            messages throughout the computational process.
%                            Logical.  Default value: Verbose = false
%                            (don't display any informational messages).
%               - bc_method  choise of global boundary condition method.
%                            Valid choises is
%                            'bc_simple' (default),
%                            'bc' need  option 'bc'
%                            'wells_simple', need option 'wells' with one
%                            well configuration and give an upscaled well
%                            cgwells as output
%                            'wells' need 'wells as' cellarray and do not
%                            upscale wells
%             - match_method how to match to the set of global boundary
%                            conditions options are
%                            'max_flux'
%                             lsq_flux'
%             - fix_trans     'true/false' set negative and zero
%                             transmissibility to lowest positive found
%             - opt_trans_alg use optimization of trans base valide values
%                             'none','local','global'
%          - opt_trans_method method for optimizing trans 'linear_simple','
%                              'convecs_simple'
%
% RETURNS:
%   HT_cg    - One-sided, upscaled transmissibilities
%
%   T_cg     - Upscaled interface transmissibilities
%
%   cgwells  - Coarse grid well with upscaled well indices
%
%   upscaled - Detailed infomation from the upscaling procedure.
%              For testing purposes.
%
% SEE ALSO:
%   coarsegrid module, `grid_structure`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('verbose',              mrstVerbose,     ...
             'bc_method',            'bc_simple',     ...
             'bc',                   [],              ...
             'wells',                [],              ...
             'match_method',         'max_flux',      ...
             'fix_trans',            false,           ...
             'opt_trans_alg',        'none',          ...
             'opt_trans_method',     'linear_simple', ...
             'check_trans',          true,            ...
             'plot_opt_grid',        true,            ...
             'use_average_pressure', true,            ...
             'use_trans',            true,            ...
             'LinSolve',             @mldivide);

opt = merge_options(opt, varargin{:});

if ~isempty(opt.opt_trans_alg) && strcmp(opt.opt_trans_alg, 'global'),
   error('global optimization not implemented')
end

% Choose global boundary conditions to use
[bc, well_cases, cgwells_cases] = setupCases(cg, opt);

% Number of cases
nr_global = numel(bc);

% Use simple, single phase fluid
fluid = initSingleFluid('mu' , 1*centi*poise, ...
                        'rho', 0*kilogram/meter^3);

% Initialise structures for global solve
init_state = initState(cg.parent, [], 100*barsa);

% Loop over global boundary conditions and wells
psolve = @(bcase, wcase) ...
   incompTPFA(init_state, cg.parent, T_fine, fluid, 'bc', bcase, ...
              'wells', wcase, 'use_trans', opt.use_trans, ...
              'LinSolve', opt.LinSolve);
states = cellfun(psolve, bc, well_cases, 'UniformOutput', false);

% Calculate upscaled properties
upscale  = @(x, cgw) calculateUpscaledSolution(cg, x, cgw, opt);
upscaled = cellfun(upscale, states, cgwells_cases, 'UniformOutput', false);

% Calulate transmissibility in place in upscaled{i}.trans, upscaled{i}.WI
utrans = @(k) calculateUpscaledTransSimple(upscaled{k}, cg, ...
                                           cgwells_cases{k}, fluid);
if strcmp(opt.opt_trans_alg, 'local'),
   utrans = @(k) localLinOptTrans(utrans(k), states{k}, ...
                                  cgwells_cases{k}, fluid, opt);
end

ext = any(cg.faces.neighbors == 0, 2);

for i = 1 : nr_global,
   upscaled{i} = utrans(i);

   if isfield(upscaled{i}, 'facePressure'),
      % Use one-sided transmissibilities to derive boundary
      % transmissibility
      t = 1 ./ accumarray(cg.cells.faces(:,1), 1 ./ upscaled{i}.htrans);

      upscaled{i}.trans(ext) = t(ext);
   end

   % Fix transmissibilities that are not finite
   atrans    = abs(upscaled{i}.trans);
   min_trans = min(atrans(atrans > 0));

   upscaled{i}.trans(~ isfinite(upscaled{i}.trans))  = min_trans * 1e-6;
   upscaled{i}.trans(~ (abs(upscaled{i}.trans) > 0)) = min_trans * 1e-6;

   if opt.check_trans,
      [upscaled{i}, ok] = checkUpscaledTrans(upscaled{i}, cg, bc{i}, ...
                                             cgwells_cases{i}, fluid);

      dispif(~ok, ['Pressure solve %d/%d with upscaled ', ...
                   'transmissibility does not give comparable ', ...
                   'results\n'], i, nr_global);
   end

   assert (all(isfinite(upscaled{i}.trans)))
end

% Make a choice of which transmisiblities to use on each face. For now use
% the one generated with largest flux
cgwells = cgwells_cases{1};
[T_cg, HT_cg, cgwells] = selectTrans(cg, upscaled, cgwells, fluid, opt);

% Fix some transmissibilities
% Fixes upcaled values for where it is undefined or if numerical error can
% make nonvalid trans NB this need a refinement.
[T_cg, HT_cg, cgwells] = fixTrans(cg, T_cg, HT_cg, cgwells, opt);

assert (all(isfinite(T_cg)))

report = struct('cg_state', upscaled, 'states', states);
end

%==========================================================================
% below is functions for subtasks
%==========================================================================

%--------------------------------------------------------------------------

% make initial CG well structure
function cgwells = makeCGWells(cg, wells)
   cgwells = wells;
   for i = 1 : numel(wells),
      fcells  = wells(i).cells;
      cgcells = cg.partition(fcells);

      tab        = sortrows([cgcells, fcells, (1 : numel(fcells)) .']);
      [cells, n] = rlencode(tab(:,1));
      fcellspos  = cumsum([1 ; n]);

      if cg.griddim > 2,
         pno = rldecode(1 : numel(cells), n, 2) .';
         cc  = cg.parent.cells.centroids(tab(:,2), 3);
         cv  = cg.parent.cells.volumes  (tab(:,2));

         hpos = sparse(pno, 1 : numel(pno), cv) ...
                * [ cc, ones([numel(pno), 1]) ];

         hpos = hpos(:,1) ./ hpos(:,2);         clear pno cc cv
      else
         hpos = 0;
      end

      cgwells(i).cells     = cells;
      cgwells(i).WI        = nan([numel(cells), 1]);
      cgwells(i).dZ        = hpos - wells(i).refDepth;
      cgwells(i).fcellspos = fcellspos;
      cgwells(i).fcells    = tab(:,2);
      cgwells(i).fperf     = tab(:,3);
   end

   cgwells(i).parent = wells(i);
end

%--------------------------------------------------------------------------

% calculate upscaled solution by averages
function upscaled = calculateUpscaledSolution(cg, state, cgwells, opt)
    cfacesno = rldecode(1:cg.faces.num, diff(cg.faces.connPos), 2) .';

    facePressure = [];

    if opt.use_average_pressure,

       % normal average pressure
       pressure = sparse(cg.partition, 1 : cg.parent.cells.num, ...
                         cg.parent.cells.volumes) ...
                  * [ state.pressure, ones([cg.parent.cells.num, 1]) ];
       pressure = pressure(:,1) ./ pressure(:,end);

       % calculate total flux over coarse face
       cfsign = fineToCoarseSign(cg);
       flux   = accumarray(cfacesno, state.flux(cg.faces.fconn) .* cfsign);

       if isfield(state, 'facePressure'),
          fp = state.facePressure   (cg.faces.fconn);
          fa = cg.parent.faces.areas(cg.faces.fconn);

          facePressure = sparse(cfacesno, 1 : numel(fa), fa) ...
                         * [ fp, ones([numel(fp), 1]) ];
          facePressure = facePressure(:,1) ./ facePressure(:,end);

          clear fp fa
       end

    else
       cg = coarsenGeometry(cg);

       % Define coarse block pressure as the fine-scale cell pressure in
       % the cell closest to the centroid of each coarse block.
       dist = @(x, y) sum((x - y) .^ 2, 2);

       df = dist(cg.cells.centroids(cg.partition,:), ...
                 cg.parent.cells.centroids);

       t  = sortrows([cg.partition, df, (1 : numel(cg.partition)) .']);
       p  = cumsum([1 ; accumarray(t(:,1), 1)]);

       pressure = state.pressure(t(p(1:end-1), 3));  clear df t p

       % Accumulate total flux over coarse face
       cfsign = fineToCoarseSign(cg);
       flux   = accumarray(cfacesno, state.flux(cg.faces.fconn) .* cfsign);

       if isfield(state, 'facePressure'),
          % Define coarse interface pressure as the fine-scale interface
          % pressure of the interface closest to the centroid of each
          % coarse face.
          fpart = rldecode(1:cg.faces.num, diff(cg.faces.connPos), 2) .';

          df = dist(cg.faces.centroids(fpart,:), ...
                    cg.parent.faces.centroids(cg.faces.fconn,:));

          t  = sortrows([fpart, df, (1 : numel(fpart)) .']);
          p  = cumsum([1 ; accumarray(t(:,1), 1)]);

          facePressure = state.facePressure(t(p(1:end-1), 3));
          clear df t p
       end
    end

    % make upscaled well and well quantities
    if ~isempty(cgwells),
       wellSol = state.wellSol;
       for i = 1 : numel(cgwells),
          % need more refinement
          pno = rldecode(1 : numel(cgwells(i).cells), ...
                         diff(cgwells(i).fcellspos), 2) .';

          wellSol(i).flux = ...
             accumarray(pno, state.wellSol(i).flux(cgwells(i).fperf));
       end
    else
       wellSol = [];
    end

    upscaled = struct('pressure',     pressure,     ...
                      'flux',         flux,         ...
                      'facePressure', facePressure, ...
                      'wellSol',      wellSol);
end

%--------------------------------------------------------------------------

% set up all the global problems to be solved
function [bc, well_cases, cgwells_cases] = setupCases(cg, opt)
   switch opt.bc_method
      case 'bc_simple'
         [pmin, pmax] = deal(100*barsa, 200*barsa);

         % Simple pressure drop in 2 (or 3) directions
         bc    = cell([2 + (cg.griddim > 2), 1]);
         bc{1} = pside([]   , cg.parent, 'West' , pmin);
         bc{1} = pside(bc{1}, cg.parent, 'East' , pmax);
         bc{2} = pside([]   , cg.parent, 'North', pmin);
         bc{2} = pside(bc{2}, cg.parent, 'South', pmax);

         if cg.griddim > 2,
            bc{3} = pside([]   , cg.parent, 'Top'   , pmin);
            bc{3} = pside(bc{3}, cg.parent, 'Bottom', pmax);
         end

         well_cases    = cell([numel(bc), 1]);
         cgwells_cases = cell([numel(bc), 1]);

      case 'bc'
         % given bc conditions
         if isempty(opt.bc),
            error('bc_method==bc need non empty bc  option')
         end

         bc            = opt.bc;
         well_cases    = cell([numel(bc), 1]);
         cgwells_cases = cell([numel(bc), 1]);

      case 'wells_simple'
         % use all linearly independent well configurations for
         % incompressible flow
         if isempty(opt.wells),
            error('bc_method==wells need non empty wells option')
         end

         assert (iscell(opt.wells))
         assert (numel(opt.wells) == 1)

         nw         = numel(opt.wells{1});
         rates      = null(ones(nw)) * (10*meter^3/day);
         w_cases    = size(rates, 2);
         well_cases = repmat(opt.wells, [w_cases, 1]);

         for c = 1 : w_cases,
            for w = 1 : (nw - 1),
               well_cases{c}(w).type = 'rate';
               well_cases{c}(w).val  = rates(w, c);
            end

            % make solution unique
            well_cases{c}(nw).type = 'bhp';
            well_cases{c}(nw).val  = 300*barsa;
         end

         bc            = cell([numel(well_cases), 1]);
         cgwells_cases = cellfun(@(wcase) makeCGWells(cg, wcase), ...
                                 well_cases, 'UniformOutput', false);

      case 'wells'
         % given wells conditions, need to be valid for incompressible flow
         if isempty(opt.wells),
            error('bc_method==''wells'' needs non-empty ''wells'' option')
         end

         if ~iscell(opt.wells),
            error('wells option must be cell array for bc_method=wells')
         end

         bc            = cell([numel(opt.wells), 1]);
         well_cases    = opt.wells;
         cgwells_cases = cellfun(@(wcase) makeCGWells(cg, wcase), ...
                                 well_cases, 'UniformOutput', false);
      otherwise
         error('This bc_method is not implemented')
   end

   assert (numel(bc) == numel(well_cases)   , 'Internal error');
   assert (numel(bc) == numel(cgwells_cases), 'Internal error');

   bc            = reshape(bc           , [], 1);
   well_cases    = reshape(well_cases   , [], 1);
   cgwells_cases = reshape(cgwells_cases, [], 1);
end

%--------------------------------------------------------------------------

% select trans from all global solves
function [T_cg, HT_cg, cgwells] = ...
      selectTrans(cg, upscaled, cgwells, fluid, opt)

   internal  = ~ any(cg.faces.neighbors == 0, 2);
   nr_global = numel(upscaled);

   dpress = @(u) diff(u.pressure(cg.faces.neighbors(internal, :)), [], 2);
   alloc  = @(nrows) NaN([nrows, nr_global]);

   flux   = alloc(cg.faces.num);
   dpf    = alloc(cg.faces.num);
   trans  = alloc(cg.faces.num);
   htrans = alloc(numel(cg.cells.faces(:,1)));
   dphf   = alloc(numel(cg.cells.faces(:,1)));

   for i = 1 : nr_global,
      trans(:,i)      = upscaled{i}.trans;
      flux(:,i)       = upscaled{i}.flux;
      htrans(:,i)     = upscaled{i}.htrans;
      dpf(internal,i) = dpress(upscaled{i});

      for j = 1 : numel(cgwells),
         fluxwells{j}(:,i) = upscaled{i}.fluxwells{j}; %#ok
         dpwells{j}(:,i)   = upscaled{i}.dpwells{j};   %#ok
      end
   end

   switch opt.match_method,
      case 'max_flux'
         % use the transmisibility of the case with largest flux over the
         % face
         [mm, choice] = max(abs(flux), [], 2); %#ok
         ch_ind = sub2ind(size(htrans), ...
                          (1 : numel(cg.cells.faces(:,1))) .',...
                          choice(cg.cells.faces(:,1)));

         HT_cg  = htrans(ch_ind);
         ch_ind = sub2ind(size(trans), (1 : cg.faces.num) .', choice);
         T_cg   = trans(ch_ind);

         cfun = @(func) cellfun(func, upscaled, 'UniformOutput', false);

         if ~isempty(cgwells),
            for j = 1 : numel(cgwells),
               bb   = cfun(@(x) x.WI{j});
               wind = [ bb{:} ];

               bb    = cfun(@(x) x.wellSol(j).flux);
               wflux = [ bb{:} ];

               [mm, choice] = max(abs(wflux), [], 2); %#ok
               ch_ind       = sub2ind(size(wind), ...
                                      (1 : size(wind,1)) .', choice);
               cgwells(j).WI = wind(ch_ind);
            end
         end

      case 'lsq_flux'
         % minimize least square flux error over the set of cases
         mu    = fluid.properties();
         HT_cg = - mu * sum(flux(cg.cells.faces(:,1),:), 2) ./ sum(dphf, 2);
         T_cg  = - mu * sum(flux, 2) ./ sum(dpf, 2);

         if ~isempty(cgwells),
            for j = 1 : numel(cgwells),
               cgwells(j).WI = sum(fluxwells{j}, 2) ./ sum(dpwells{j}, 2);
               cgwells(j).WI = cgwells(j).WI * mu;
            end
         end

         T_tmp = 1 ./ accumarray(cg.cells.faces(:,1), ...
                                 1 ./ HT_cg, [cg.faces.num, 1]);
         T_cg(~internal) = T_tmp(~internal);

         if any(~isfinite(HT_cg)) || any(~isfinite(T_cg)),
            warning(msgid(''), ...
                    'Some transmisibilities undecided set to zero');
            HT_cg(~ isfinite(HT_cg)) = 0;
            T_cg (~ isfinite(T_cg))  = 0;
         end

      otherwise
         error('This match method is not implemented')
   end
end

%--------------------------------------------------------------------------

% check solution by solving the coarse system with upscaled trans
function [upscaled, ok] = ...
      checkUpscaledTrans(upscaled, cg, bc, well_cases, fluid)
   if isempty(bc),
      cgbc = [];
   else
      % Extract subfaces on coarse interfaces
      [nsub, sub]     = subFaces(cg.parent, cg);
      [sgn, coarse_f] = signOfFineFacesOnCoarseFaces(cg.parent, cg, ...
                                                     nsub, sub);

      % Coarse boundary conditions
      cgbc = convertBC2Coarse(bc, cg.parent, cg, nsub, sub, coarse_f, sgn);
   end

   cgwells_tmp = well_cases;
   if ~isempty(cgwells_tmp) && isfield(upscaled, 'WI') && ...
         iscell(upscaled.WI) && (numel(upscaled.WI) == numel(cgwells_tmp)),
      [ cgwells_tmp.WI ] = upscaled.WI{:};
   end

   cgtrans = upscaled.trans;

   % do not manage to to the check for all cases.

   internal = ~any(cg.faces.neighbors == 0, 2);

   cgstate = struct('flux',     upscaled.flux,     ...
                    'pressure', upscaled.pressure, ...
                    's',        zeros([cg.cells.num, 1]));

   cgstate = incompTPFA(cgstate, cg, cgtrans, fluid, 'bc', cgbc, ...
                        'wells', cgwells_tmp, 'use_trans', true);
%{
   if(max(cgstate.pressure)>min(cgstate.pressure))
      assert(all(abs((cgstate.pressure-upscaled.pressure)/(max(cgstate.pressure) - min(cgstate.pressure)))<1e-6))
   end
   if(any(upscaled.flux(internal)))
      assert(all(abs(cgstate.flux(internal)-upscaled.flux(internal))./max(abs(cgstate.flux))< 1e-6))
   end
%}
    ok = all(abs(cgstate.flux(internal) - upscaled.flux(internal)) ./ ...
             max(abs(cgstate.flux)) < 1e-6);

    if ~isempty(cgwells_tmp),
       ok = ok && ...
          all(abs(vertcat(upscaled.wellSol.flux) - ...
                  vertcat(cgstate.wellSol.flux)) / ...
              max(abs(vertcat(cgstate.wellSol.flux))) < 1e-6);

       for k = 2 : numel(cgwells_tmp),
          dc = cgstate.wellSol(k) .pressure - cgstate .wellSol(1).pressure;
          df = upscaled.wellSol(k).pressure - upscaled.wellSol(1).pressure;

          ok = ok && (abs(dc - df) / abs(df) < 1e-6);
       end
    end

    dispif(~ ok, 'something is not ok.\n')
end

%--------------------------------------------------------------------------

% Calulate effective transmissibility from cg fluxes and pressures
function upscaled = ...
      calculateUpscaledTransSimple(upscaled, cg, cgwells, fluid)

   mu = fluid.properties();  % Used throughout function.

   int = ~ any(cg.faces.neighbors == 0, 2);
   dp  = diff(upscaled.pressure(cg.faces.neighbors(int, :)), [], 2);

   upscaled.trans      = nan([cg.faces.num, 1]);
   upscaled.trans(int) = - mu .* (upscaled.flux(int) ./ dp);

   clear dp int

   if isfield(upscaled, 'facePressure'),
      % Calculate one-sided transmissibilities if possible
      cellno = rldecode(1 : cg.cells.num, diff(cg.cells.facePos), 2) .';
      cf     = cg.cells.faces(:,1);

      dp    = upscaled.facePressure(cf) - upscaled.pressure(cellno);
      sgn   = 2*(cg.faces.neighbors(cf, 1) == cellno) - 1;
      hflux = sgn .* upscaled.flux (cf);

      upscaled.htrans = - mu .* (hflux ./ dp);

      clear cellno cf dp sgn hflux
   end

   if ~isempty(cgwells),
      assert (numel(upscaled.wellSol) == numel(cgwells), ...
              'Internal error defning upscaled well quantities');

      % Calculate upscaled well indices
      nperf = cellfun('prodofsize', { cgwells.cells });
      pos   = cumsum([ 1 ; reshape(nperf, [], 1) ]);

      press = rldecode(vertcat(upscaled.wellSol.pressure), nperf);
      dp    = press - upscaled.pressure(vertcat(cgwells.cells));
      flux  = vertcat(upscaled.wellSol.flux);

      WI    = mu .* (flux ./ dp);

      for j = 1 : numel(cgwells),
         ix = pos(j) : pos(j + 1) - 1;

         upscaled.fluxwells{j} = flux(ix);
         upscaled.dpwells{j}   = dp  (ix);
         upscaled.WI{j}        = WI  (ix);
      end
   end
end

%--------------------------------------------------------------------------

% fix some tansmissibilities
function [T_cg, HT_cg, cgwells] = fixTrans(cg, T_cg, HT_cg, cgwells, opt)
   ext = any(cg.faces.neighbors == 0, 2);

   if any(strcmp(opt.bc_method, {'wells', 'wells_simple'})),
      warning('upscaleTransNew:ZeroBoundary', ...
          'Set boundary face trans to zero');
      hf_ext        = ext(cg.cells.faces(:,1));
      HT_cg(hf_ext) = 0;
      T_cg(ext)     = 0;
   end

   if opt.fix_trans,
      % Set negative and zero transmissibility to the minimum positive
      % (one-sided) transmissibility.
      HT_min = min(HT_cg(HT_cg > 0));
      HT_cg  = max(HT_cg, HT_min);
      T_cg   = max(T_cg , HT_min);

      if ~isempty(cgwells),
         WI = cellfun(@(x) max(x, 0), { cgwells.WI }, ...
                      'UniformOutput', false);
         [ cgwells.WI ] = WI{:};
      end
   end

   if ~ (isempty(opt.wells) || iscell(opt.wells)),
      assert (numel(opt.wells) == numel(cgwells));

      [ cgwells.type ] = opt.wells.type;
      [ cgwells.val  ] = opt.wells.val;
   end
end
