function W = processWells(G, rock, control, varargin)
%Construct MRST well structure from ECLIPSE input deck control.
%
% SYNOPSIS:
%   W = processWells(G, rock, control)
%   W = processWells(G, rock, control, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   rock    - Rock data structure.  Must contain valid field 'perm'.
%
%   control - Well control mode.  Assumed to be defined from function
%             'readSCHEDULE'.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%           - InnerProduct -- The inner product with which to define
%                             the mass matrix.
%                             String.  Default value = 'ip_tpf'.
%             Supported values are:
%               - 'ip_simple'
%               - 'ip_tpf'
%               - 'ip_quasitpf'
%               - 'ip_rt'
% 
%           - Verbose      -- Whether or not to emit informational messages
%                            if any requested well control mode is
%                            unsupported.
%                            Logical.  Default value = false.
%           - DepthReorder -- Boolean indicating if we should attempt to
%                             reorder the perforations by depth. Should
%                             only be used for vertical wells. Default off.
%           - strictParsing -- Throw errors if unsupported well controls
%                              are encountered. If disabled, wells with
%                              unsupported controls will recieve NaN values
%                              for control and phase composition.
%           - createDefaultWell -- If enabled, this function will not parse
%                              any wells, but simply return a default well
%                              structure suitable for further modification.
%           - cellDims -- optional input of accurate cellDims
%
% RETURNS:
%   W  - Updated (or freshly created) well structure, each element of which
%        has the following fields:
%        cells   -- Grid cells perforated by this well.
%        type    -- Well control type.
%        val     -- Target control value.
%        r       -- Well bore radius.
%        dir     -- Well direction.
%        WI      -- Well index.
%        dZ      -- Vertical displacement of each well perforation measured
%                   from 'highest' horizontal contact (i.e., the 'TOP'
%                   contact with the minimum 'Z' value counted amongst all
%                   cells perforated by this well).
%       name     -- Well name.
%       compi    -- Fluid composition--only used for injectors.
%       sign     -- injection (+) or production (-) flag.
%       status   -- Boolean indicating if the well is open or shut.
%       cstatus  -- One entry per cell, indicating if the completion is
%                   open.
%       lims     -- Limits for the well. Contains subfields for the types
%                   of limits applicable to the well (bhp, rate, orat, ...)
%                   Injectors generally have upper limits, while producers
%                   have lower limits.
%       c        -- Component concentrations. Unit convention may vary.
%
% SEE ALSO:
%   `readSCHEDULE`, `readWellKW`, `addWell`.

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

   opt = struct('InnerProduct',     'ip_tpf', ...
                'Verbose',           mrstVerbose(), ...
                'DepthReorder',      false, ...
                'strictParsing',     true, ...
                'createDefaultWell', false, ...
                'outputDefaulted',   false, ...
                'cellDims',         []);
   opt = merge_options(opt, varargin{:});

   if opt.createDefaultWell
      W = createDefaultWell(G, rock);
      return
   end

   well_id      = enumerateWells(control);
   control      = insertDefaults(control, well_id, G.cartDims(3));
   [control, p] = orderCompletions(control);

   inj  = struct('WCONINJ' , @process_wconinj , ...
                 'WCONINJE', @process_wconinje, ...
                 'WCONINJH', @process_wconinjh);

   prod = struct('WCONHIST', @process_wconhist, ...
                 'WCONPROD', @process_wconprod);

   post = struct('WPOLYMER', @process_wpolymer, ...
                 'WSURFACT', @process_wsurfact, ...
                 'WSOLVENT', @process_wsolvent, ...
                 'WTEMP'   , @process_wtemp);

   % ----------------------------------------------------------------------
   % Well processing stages
   %
   stages = { inj, prod, post };

   W = [];

   c2a = make_cart_to_active(G);

   for stage = reshape(stages, 1, [])
      fn = reshape(sort(fieldnames(stage{1})), 1, []);

      for kw = fn(isfield(control, fn))
         W = stage{1}.(kw{1}) (W, control, G, rock, c2a, well_id, p, opt);
      end
   end

   if isempty(W)
       if opt.createDefaultWell
           W = createDefaultWell(G, rock);
       else
           W = [];
       end
      return
   end

   if isempty(opt.DepthReorder) 
      active = arrayfun(@(x) strcmpi(x.dir(1), 'z'), W);
      W = reorderWellPerforationsByDepth(W, active);
   elseif opt.DepthReorder
      W = reorderWellPerforationsByDepth(W);
   end

   % Close all completions if well inactive (W.status == false).
   cstat = arrayfun(@(w) w.status & w.cstatus, ...
                    W, 'UniformOutput', false);
   [W.cstatus] = cstat{:};

   % Close all wells for which all completions are inactive.
   stat = arrayfun(@(w) w.status && any(w.cstatus), ...
                   W, 'UniformOutput', false);
   [W.status] = stat{:};
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function c2a = make_cart_to_active(G)
   c2a = zeros([prod(G.cartDims), 1]);
   c2a(G.cells.indexMap) = 1 : G.cells.num;
end

%--------------------------------------------------------------------------

function W = process_wconinj(W, control, G, rock, c2a, well_id, p, opt)
   for i = 1 : size(control.WCONINJ, 1)
      nm = control.WCONINJ{i,1};
      status = strcmp(control.WCONINJ{i,3}, 'OPEN');

      type = lower(control.WCONINJ{i,4});
      switch type
        case 'rate', val = control.WCONINJ{i, 5};
        case 'resv', val = control.WCONINJ{i, 6};
        case 'bhp' , val = control.WCONINJ{i, 9};
        case 'thp' , val = control.WCONINJ{i, 10};

        otherwise
          dispif(opt.Verbose && status, ...
                 ['Control mode ''%s'' unsupported for injector ', ...
                  '''%s''.  Ignored.\n'], upper(type), nm);
          continue
      end

      % Define injection fluid.  This is a hack.
      switch lower(control.WCONINJ{i,2}(1))
        case 'w', compi = [1, 0, 0];  % Water, 1st phase
        case 'o', compi = [0, 1, 0];  % Oil, 2nd phase
        case 'g', compi = [0, 0, 1];  % Gas, 3rd phase
        otherwise
          dispif(opt.Verbose, ...
                 ['Injection phase ''%s'' is unknown.  ', ...
                  'Well ignored.\n'], control.WCONINJ{i,2});
          continue
      end

      sizeW = numel(W);
      W = buildWell(W, G, rock, c2a, control, well_id(nm), ...
                    p, type, val, compi, opt.InnerProduct, 1, opt);

      if numel(W) > sizeW
         W(end).lims.rate = control.WCONINJ{i, 5};
         W(end).lims.bhp  = control.WCONINJ{i, 9};
         W(end).lims.thp  = control.WCONINJ{i, 10};
         W(end).status    = status;
      end
      W(end).vfp_index = control.WCONINJ{i, 11};
   end
end

%--------------------------------------------------------------------------

function W = process_wconinje(W, control, G, rock, c2a, well_id, p, opt)
   for i = 1 : size(control.WCONINJE, 1)
      nm = control.WCONINJE{i,1};
      status = strcmp(control.WCONINJE{i,3}, 'OPEN');

      type = lower(control.WCONINJE{i,4});
      switch type
        case 'rate'
          val = control.WCONINJE{i, 5};
        case 'resv', val = control.WCONINJE{i, 6};
        case 'bhp' , val = control.WCONINJE{i, 7};
        case 'thp' , val = control.WCONINJE{i, 8};
        otherwise
          dispif(opt.Verbose, ...
                 ['Control mode ''%s'' unsupported for injector ', ...
                  '''%s''.  Ignored.\n'], upper(type), nm);
          %continue
          val = 0;
      end


      % Define injection fluid.  This is a hack.
      switch lower(control.WCONINJE{i,2}(1))
        case 'w', compi = [1, 0, 0];  % Water, 1st phase
        case 'o', compi = [0, 1, 0];  % Oil, 2nd phase
        case 'g', compi = [0, 0, 1];  % Gas, 3rd phase
        otherwise
          dispif(opt.Verbose, ...
                ['Injection phase ''%s'' is unknown.  ', ...
                 'Well ignored.\n'], control.WCONINJE{i,2});
          continue
      end

      sizeW = numel(W);
      W = buildWell(W, G, rock, c2a, control, well_id(nm), ...
                    p, type, val, compi, opt.InnerProduct, 1, opt);

      if numel(W) > sizeW
         W(end).lims.rate = control.WCONINJE{i, 5};
         W(end).lims.bhp  = control.WCONINJE{i, 7};
         W(end).lims.thp  = control.WCONINJE{i, 8};
         W(end).status    = status;
      end
      W(end).vfp_index = control.WCONINJE{i, 9};
   end
end

%--------------------------------------------------------------------------

function W = process_wconinjh(W, control, G, rock, c2a, well_id, p, opt)
   for i = 1 : size(control.WCONINJH, 1)
      nm = control.WCONINJH{i,1};
      status = strcmp(control.WCONINJH{i,3}, 'OPEN');

      type = 'rate'; % History matched well are always rate controlled
      val = control.WCONINJH{i, 4};

      % Define injection fluid.  This is a hack.
      switch lower(control.WCONINJH{i,2}(1))
        case 'w', compi = [1, 0, 0];  % Water, 1st phase
        case 'o', compi = [0, 1, 0];  % Oil, 2nd phase
        case 'g', compi = [0, 0, 1];  % Gas, 3rd phase
        otherwise
          dispif(opt.Verbose, ...
                 ['Injection phase ''%s'' is unknown.  ', ...
                  'Well ignored.\n'], control.WCONINJH{i,2});
          continue
      end

      sizeW = numel(W);
      W = buildWell(W, G, rock, c2a, control, well_id(nm), ...
                    p, type, val, compi, opt.InnerProduct, 1, opt);

      if numel(W) > sizeW
         % lims.bhp is \approx 100e3 psia, as for default wconinje.

         W(end).lims.rate = control.WCONINJH{i, 4};
         W(end).lims.bhp  = 6895 * barsa;
         W(end).status    = status;
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wpolymer(W, control, varargin)
   npoly = size(control.WPOLYMER, 1);

   if npoly == 0, return, end

   if ~isfield(W, 'cp')
       [W.cp] = deal([]);
   end

   for i = 1:numel(W)
       % Add zeros for polymer
       W(i).cp = [W(i).cp, 0];
   end

   for i = 1 : npoly
      for j = 1:size(W,1)
         if strcmp(W(j).name, control.WPOLYMER{i,1})
            W(j).cp(end) = control.WPOLYMER{i,2};
         end
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wsurfact(W, control, varargin)
   nsurf = size(control.WSURFACT, 1);

   if nsurf == 0, return, end

   if ~isfield(W, 'cs')
       [W.cs] = deal([]);
   end

   for i = 1:numel(W)
       % Add zeros for surfactant
       W(i).cs = [W(i).cs, 0];
   end

   if ~isempty(W)
      Wn = { W.name };

      for i = 1 : nsurf
         j = find(strcmp(Wn, control.WSURFACT{i,1}));

         if ~isempty(j)
            [ W(j).cs(end) ] = deal(control.WSURFACT{i,2});
         end
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wsolvent(W, control, varargin)
   nsol = size(control.WSOLVENT, 1);

   if nsol == 0, return, end

   if ~isfield(W, 'solventFrac')
       [W.solventFrac] = deal([]);
   end

   for i = 1:numel(W)
       % Add zeros for solvent fraction
       W(i).solventFrac = 0;
   end

   if ~isempty(W)
      Wn = { W.name };

      for i = 1 : nsol
         j = find(strcmp(Wn, control.WSOLVENT{i,1}));

         if ~isempty(j)
            [ W(j).solventFrac ] = deal(control.WSOLVENT{i,2});
         end
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wtemp(W, control, varargin)
   nsurf = size(control.WTEMP, 1);

   if nsurf == 0, return, end

   if ~isfield(W, 'T')
       [W.T] = deal([]);
   end

   for i = 1:numel(W)
       % Add zeros for surfactant
       W(i).T = [W(i).T, 0];
   end

   if ~isempty(W)
      Wn = { W.name };

      for i = 1 : nsurf
         j = find(strcmp(Wn, control.WTEMP{i,1}));

         if ~isempty(j)
            [ W(j).T(end) ] = deal(control.WTEMP{i,2});
         end
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wconhist(W, control, G, rock, c2a, well_id, p, opt)
   for i = 1 : size(control.WCONHIST, 1)
      nm = control.WCONHIST{i,1};
      status = strcmp(control.WCONHIST{i,2}, 'OPEN');

      type = lower(control.WCONHIST{i,3});
      switch type
        case 'orat'
          val = - control.WCONHIST{i, 4};
          compi = [0, 1, 0];  % Oil, 2nd phase
        case 'wrat'
          val = - control.WCONHIST{i, 5};
          compi = [1, 0, 0];  % Water, 1st phase
        case 'grat'
          val = - control.WCONHIST{i, 6};
          compi = [0, 0, 1];  % Gas, 3rd phase
        case 'lrat'
          rates = - ([control.WCONHIST{i, 4:5}]);
          val   = sum(rates);
          if val ~= 0
            compi = [rates, 0]./val;
          else
            compi = [1 1 0]/2;
          end
        case 'resv'
          rates = - ([control.WCONHIST{i, 4:6}]);
          val   = sum(rates);
          if val ~= 0
              % Account for OWG ordering. MRST uses WOG.
              compi = rates([2, 1, 3])./val;
          else
              compi = [1 1 1]/3;
          end
          % Historical RESV control
          type = 'RESV_HISTORY';
        case 'bhp'
          val = control.WCONHIST{i, 10};

        otherwise
          dispif(opt.Verbose && status, ...
                 ['Control mode ''%s'' unsupported for producer ', ...
                  '''%s''.  Ignored.\n'], upper(type), nm);
          val = 0;
      end

      sizeW = numel(W);

      W = buildWell(W, G, rock, c2a, control, well_id(nm), ...
                    p, type, val, compi, opt.InnerProduct, -1, opt);

      if numel(W) > sizeW
          W(end).lims = struct('orat', -inf, ...
                               'wrat', -inf, ...
                               'grat', -inf, ...
                               'lrat', -inf, ...
                               'resv', -inf, ...
                               'bhp', 1*atm);
         switch type
           case 'orat'
             W(end).lims.orat = -control.WCONHIST{i, 4};
           case 'wrat'
             W(end).lims.wrat = -control.WCONHIST{i, 5};
           case 'grat'
             W(end).lims.grat = -control.WCONHIST{i, 6};
           case 'lrat'
             W(end).lims.lrat = -sum([control.WCONHIST{i, 4:5}]);
           case 'bhp'
             W(end).lims.bhp  = control.WCONHIST{i, 10};
         end
         W(end).status    = status;
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wconprod(W, control, G, rock, c2a, well_id, p, opt)
   for i = 1 : size(control.WCONPROD,1)
      nm = control.WCONPROD{i,1};
      status = strcmp(control.WCONPROD{i,2}, 'OPEN');

      type = lower(control.WCONPROD{i,3});
      switch type
        case 'orat'
          % 'val' is INJECTION rate.
          val = -control.WCONPROD{i, 4};
          compi = [ 0, 1, 0 ];

        case 'wrat'
          % 'val' is INJECTION rate.
          val = -control.WCONPROD{i, 5};
          compi = [ 1, 0, 0 ];

        case 'grat'
          % 'val' is INJECTION rate.
          val = -control.WCONPROD{i, 6};
          compi = [ 0, 0, 1 ];

        case 'lrat'
          % 'val' is INJECTION rate.
          val   = -control.WCONPROD{i, 7};
          compi = [ 1, 1, 0 ];  % LIQUID rate == water+oil rate (surface)

        case 'resv'
          % 'val' is INJECTION rate.
          val   = -control.WCONPROD{i, 8};
          compi = [ 1, 1, 1 ];

        case 'bhp'
          val   = control.WCONPROD{i, 9};
          compi = [1, 1, 1];  % Doesn't matter.

        case 'thp'
          val   = control.WCONPROD{i, 10};
          compi = [1, 1, 1];  % Doesn't matter.

        otherwise
          dispif(opt.Verbose && status, ...
                 ['Control mode ''%s'' unsupported for producer ', ...
                  '''%s''.  Ignored.\n'], upper(type), nm);
          val = 0;
          compi = [0 1 0]; 
      end

      sizeW = numel(W);

      W = buildWell(W, G, rock, c2a, control, well_id(nm), ...
                    p, type, val, compi, opt.InnerProduct, -1, opt);

      if numel(W) > sizeW
         W(end).lims.orat = -control.WCONPROD{i,  4};
         W(end).lims.wrat = -control.WCONPROD{i,  5};
         W(end).lims.grat = -control.WCONPROD{i,  6};
         W(end).lims.lrat = -control.WCONPROD{i,  7};
         W(end).lims.bhp  =  control.WCONPROD{i,  9};
         W(end).lims.thp  =  control.WCONPROD{i, 10};
         W(end).status    =  status;
      end
      W(end).vfp_index = control.WCONPROD{i, 11};
   end
end

%--------------------------------------------------------------------------

function wid = enumerateWells(control)
   [wn, i] = sort(control.WELSPECS(:,1));
   assert (numel(unique(wn)) == numel(wn));

   try
      % Java's O(1) hash table string search support.
      ht = java.util.Hashtable;

      for n = 1 : numel(wn)
         ht.put(wn{n}, i(n));
      end

      wid = @(s) ht.get(s);

   catch
      % Fall back to (probably) linear structure field name search if Java
      % is unavailable (e.g., -nojvm or a different interpreter).

      ht = struct();

      for n = 1 : numel(wn)
         ht.(regexprep(wn{n}, '\W', '_')) = i(n);
      end

      wid = @(s) ht.(regexprep(s, '\W', '_'));
   end
end

%--------------------------------------------------------------------------

function control = insertDefaults(control, id, nlayers)
   control.COMPDAT  = insertDefaultCOMPDAT (control, id, nlayers);
   control.WCONINJE = insertDefaultWCONINJE(control);
   control.WCONPROD = insertDefaultWCONPROD(control);
end

%--------------------------------------------------------------------------

function [control, p] = orderCompletions(control)
   if isempty(control.COMPDAT)
       p = [];
       return;
   end
   comp = [vertcat(control.COMPDAT{:,1}), ...
           (1 : size(control.COMPDAT,1)).'];
   comp = sortrows(comp);

   [w, p] = rlencode(comp(:,1), 1);
   assert (all(diff(w) == 1));

   p = cumsum([1; p]);
   control.COMPDAT = control.COMPDAT(comp(:,2), :);  % Reorder completions.
end

%--------------------------------------------------------------------------

function W = buildWell(W, G, rock, c2a, control, i, p, ...
                       type, val, compi, ip, sgn, opt)
   comp  = control.COMPDAT(p(i) : p(i + 1) - 1, :);
   name  = control.WELSPECS{i, 1};
   perf  = arrayfun(@(i) active_perf(G, comp(i,:), c2a), ...
                    (1 : size(comp,1)).', 'UniformOutput', false);

   nperf = reshape(cellfun('prodofsize', perf), [], 1);

   wdir  = rldecode([comp{:,13}], nperf, 2) .';
   Kh    = rldecode([comp{:,10}], nperf, 2) .';
   WI    = rldecode([comp{:, 8}], nperf, 2) .';
   Wdiam = rldecode([comp{:, 9}], nperf, 2) .';
   openShutFlag = rldecode(comp(:, 6), nperf);

   % Zero input value means to calculate quantities ourselves.
   % Set negative whence 'addWell' will DTRT[tm].
   %
   defaultedKh = Kh <= 0;
   defaultedWI = WI <= 0;
   Kh(defaultedKh) = -1;
   WI(defaultedWI) = -1;

   Wdiam(~ (Wdiam > 0)) = 1*ft;

   assert (all(Wdiam(defaultedKh & defaultedWI) > 0), ...
          ['Well bore diameter must be defined unless productivity ', ...
           'index is explicitly provided.']);

   RefDepth = control.WELSPECS{i,5};
   assert (isnumeric(RefDepth));

   defaultedDepth = isnan(RefDepth);
   if defaultedDepth
      RefDepth = min(G.cells.centroids(vertcat(perf{:}), 3));
   end

   % Remove cells from connection which are shutdown
   perf = vertcat(perf{:});
   [ia, ia] = unique(perf, 'last');                             %#ok<ASGLU>
   ia = sort(ia); % Don't sort well connections according to cell ID

   cstatus = ~(strcmp('SHUT', openShutFlag(ia)) | strcmp('HEAT', openShutFlag(ia)));
   %ia = ia(ia_open);

   perf  =  perf(ia);
   wdir  =  wdir(ia);
   Kh    =  Kh(ia);
   WI    =  WI(ia);
   Wdiam =  Wdiam(ia);

   if sum(nperf) > 0
      sizeW = numel(W);
      W = addWell(W, G, rock, perf,                                   ...
                  'Type', type, 'Val', val, 'Dir', wdir, 'Kh', Kh,    ...
                  'Sign', sgn, 'Wi', WI, 'Comp_i', compi,             ...
                  'Name', name, 'Radius', Wdiam./2, ...
                  'RefDepth', RefDepth, 'InnerProduct', ip,           ...
                  'cellDims', opt.cellDims, 'heel', [control.WELSPECS{i, 3:4}]);

      if numel(W) > sizeW
         W(end).lims = [];
         W(end).cstatus = cstatus;
      end
   end
   if opt.outputDefaulted
       W(end).defaulted = struct('Kh', defaultedKh, ...
                                 'WI', defaultedWI, ...
                                 'refDepth', defaultedDepth);
   end
end

%--------------------------------------------------------------------------

function compdat = insertDefaultCOMPDAT(control, id, nlayers)
   compdat = control.COMPDAT;
   present = false([size(control.WELSPECS,1), 1]);

   if ~isempty(compdat)
      compdat(:,1) = cellfun(id, compdat(:,1), 'UniformOutput', false);
      present(vertcat(compdat{:,1})) = true;
   end

   % Wells for which no completions have been explicitly set using the
   % 'COMPDAT' input keyword get vertical completions in all layers.
   %
   tmpl = {-1, -1, 1, nlayers, 'OPEN', 0, 0.0, 0.1, -1, 0.0, ...
           'Default', 'Z', -1};
   for w = reshape(find(~present), 1, [])
      compdat = [compdat; [w, tmpl]]; %#ok
   end
   if isempty(compdat)
       return;
   end
   % Fill in defaulted (I,J) locations (of well heel).
   i = vertcat(compdat{:,2}) < 1;
   compdat(i,2) = control.WELSPECS(vertcat(compdat{i,1}), 3);

   i = vertcat(compdat{:,3}) < 1;
   compdat(i,3) = control.WELSPECS(vertcat(compdat{i,1}), 4);
end

%--------------------------------------------------------------------------

function wconinje = insertDefaultWCONINJE(control)
   wconinje = control.WCONINJE;

   if numel(wconinje) > 0
      % Insert default pressure target (large value) where needed
      i = isnan(vertcat(wconinje{:,7}));
      wconinje(i,7) = { 6895 * barsa };  % \approx 1e5 psia
   end
end

%--------------------------------------------------------------------------

function wconprod = insertDefaultWCONPROD(control)
   wconprod = control.WCONPROD;

   if numel(wconprod) > 0
      % Insert default pressure value (atmospheric pressure) where needed.
      i = isnan(vertcat(wconprod{:,9}));
      wconprod(i,9) = { 1 * atm };
   end
end

%--------------------------------------------------------------------------

function perf = active_perf(G, comp, c2a)
   [ijk{1:3}] = ndgrid(comp{2}, comp{3}, comp{4} : comp{5});

   perf = c2a(sub2ind(G.cartDims, ijk{1}(:), ijk{2}(:), ijk{3}(:)));
   perf = perf(perf > 0);
end

%--------------------------------------------------------------------------

function W = createDefaultWell(G, rock)
% default well will be handled as injector in checkLims
   rock.perm = ones([G.cells.num, 1]);

   W = addWell([], G, rock, 1, 'Val', 0, 'Type', 'rate', ...
               'sign', 1, 'Comp_i', 1/3*[1, 1, 1], ...
               'refDepth', G.cells.centroids(1,3));

   W.lims.rate = -inf;
   W.lims.bhp  =  inf;
end
