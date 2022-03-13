function res = SampledProp2D(name, file, varargin)
% Create a structure with functions to interpolate a 2D sampled property and its
% derivatives, optionally with a discontinuity (phase boundary)
%
% SYNOPSIS:
%   function res = SampledProp2D(name, file, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   name     - name of property (e.g. 'density', 'viscosity', 'rho', ...)
%   file     - filename containing the sampled values as MATLAB objects.  The
%              objects should be:
%              * A 2D matrix named 'vals', representing the actual sampled
%                values
%              * Two structs, representing the two variables for which the
%                table has been sampled.  The structs contain the fields:
%                - num      - number of sample points for this variable
%                - span     - vector with 2 components, representing the range
%                             of the variable
%                - stepsize - the size of each step.  Should equal diff(span)/num.
%
%   varargin - optional parameters (as key/value pairs) are:
%                - 'assert_in_range':   If 'true', throw an error if user tried
%                                       to extrapolate outside valid range.
%                - 'nan_outside_range': If 'true' (default), return NaN
%                                       values outside valid range.
%                                       Otherwise, extrapolate as constant function.
%                - 'phase_boundary':    Cell array.  If nonempty, the first
%                                       cell contains the parameter values for
%                                       the critical point (start of the
%                                       discontinuity line).  The second cell
%                                       contains a function that describes the
%                                       dicontinuity line as a relation between
%                                       the two variables, on the form v2 = f(v1).
%                - 'const_derivatives': If 'true' (default), only first-order
%                                       partial derivatives are included.
%                                       Otherwise, second-order partial
%                                       derivatives are included as well.
%                                       For code based on automatic
%                                       differentiation, it is recommended to
%                                       avoid the second-order partial
%                                       derivatives.
%
% RETURNS:
%   res - Struct containing the generated functions.
%

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

   opt.assert_in_range   = true; % throw error if sampled outside range
   opt.nan_outside_range = true; % return NaN values if sampled outside range
   opt.phase_boundary    = [];   % {[critical point], boundary fun v1 = f(v2)}
   opt.const_derivatives = true; % consider constant derivative between each
                                 % sample point (precludes functions with
                                 % higher derivatives)
   opt = merge_options(opt, varargin{:});

   prop   = load(file);
   vn     = get_varnames(prop);
   tables = establish_data_tables(prop, vn, opt.phase_boundary, opt.const_derivatives);

   % setup the actual interpolation functions
   cfuns = setup_cfuns(tables, opt.const_derivatives);

   % Adding interpolating functions to the result object
   d1name = [name, 'D', vn{1}];
   d2name = [name, 'D', vn{2}];

   res.(name)   = @(v1, v2) calcMulti(v1, v2, cfuns.calcVal, cfuns.calcD1,  cfuns.calcD2);
   res.(d1name) = @(v1, v2) calcMulti(v1, v2, cfuns.calcD1,  cfuns.calcD11, cfuns.calcD12);
   res.(d2name) = @(v1, v2) calcMulti(v1, v2, cfuns.calcD2,  cfuns.calcD12, cfuns.calcD22);

   if ~opt.const_derivatives
      % Include higher derivative functions
      res.([name, 'D', vn{1}, vn{1}]) = @(v1, v2) calcMulti(v1, v2, cfuns.calcD11, cfuns.calcD111, cfuns.calcD112);
      res.([name, 'D', vn{2}, vn{2}]) = @(v1, v2) calcMulti(v1, v2, cfuns.calcD22, cfuns.calcD221, cfuns.calcD222);
      res.([name, 'D', vn{1}, vn{2}]) = @(v1, v2) calcMulti(v1, v2, cfuns.calcD12, cfuns.calcD112, cfuns.calcD122);
      res.([name, 'D', vn{1}, vn{1}, vn{1}]) = @(v1, v2) calcMulti(v1, v2, cfuns.calcD111, @noder, @noder);
      res.([name, 'D', vn{2}, vn{2}, vn{2}]) = @(v1, v2) calcMulti(v1, v2, cfuns.calcD222, @noder, @noder);
      res.([name, 'D', vn{1}, vn{1}, vn{2}]) = @(v1, v2) calcMulti(v1, v2, cfuns.calcD112, @noder, @noder);
      res.([name, 'D', vn{1}, vn{2}, vn{2}]) = @(v1, v2) calcMulti(v1, v2, cfuns.calcD122, @noder, @noder);
   end

   res.([name, 'Data']) = @() dataObjects(tables);
   res.phaseOf = @(v1, v2) phase_of(v1, v2, opt.phase_boundary);

   % ----------------------------------------------------------------------------

   function cfuns = setup_cfuns(t, const_derivatives)%#ok

   % Calculating the interpolated values
      cfuns.calcVal = @(v1, v2) extractValues(v1, v2, t.v1.span, t.v2.span, t.vals, t.vals, t.vals, false, false);

      if opt.const_derivatives
         %% We only need first partial derivatives, and do not shrink domain
         cfuns.calcD1   = @(v1, v2) extractValues(v1, v2, t.v1.span, t.v2.span, t.sup.d1, t.liq.d1, t.gas.d1, true, false);
         cfuns.calcD2   = @(v1, v2) extractValues(v1, v2, t.v1.span, t.v2.span, t.sup.d2, t.liq.d2, t.gas.d2, false, true);
         cfuns.calcD11  = @noder;
         cfuns.calcD22  = @noder;
         cfuns.calcD12  = @noder;
         cfuns.calcD111 = @noder;
         cfuns.calcD222 = @noder;
         cfuns.calcD112 = @noder;
         cfuns.calcD122 = @noder;
      else
         %% First partial derivative interpolating functions
         v1span = shrink_span(t.v1.span, t.v1.stepsize);
         v2span = shrink_span(t.v2.span, t.v2.stepsize);
         cfuns.calcD1 = @(v1, v2) extractValues(v1, v2, v1span, t.v2.span, t.sup.d1, t.liq.d1, t.gas.d1, false, false);
         cfuns.calcD2 = @(v1, v2) extractValues(v1, v2, t.v1.span, v2span, t.sup.d2, t.liq.d2, t.gas.d2, false, false);

         %% Second partial derivative interpolating functions
         v1span = shrink_span(t.v1.span, 2 * t.v1.stepsize);
         v2span = shrink_span(t.v2.span, 2 * t.v2.stepsize);
         cfuns.calcD11 = @(v1, v2) extractValues(v1, v2, v1span, t.v2.span, t.sup.d11, t.liq.d11, t.gas.d11, false, false);
         cfuns.calcD22 = @(v1, v2) extractValues(v1, v2, t.v1.span, v2span, t.sup.d22, t.liq.d22, t.gas.d22, false, false);

         %% Cross derivative
         v1span = shrink_span(t.v1.span, t.v1.stepsize);
         v2span = shrink_span(t.v2.span, t.v2.stepsize);
         cfuns.calcD12 = @(v1, v2) extractValues(v1, v2, v1span, v2span, t.sup.d12, t.liq.d12, t.gas.d12, false, false);

         %% Third order partial derivatives
         v1span = shrink_span(t.v1.span, 3 * t.v1.stepsize);
         v2span = shrink_span(t.v2.span, 3 * t.v2.stepsize);
         cfuns.calcD111 = @(v1, v2) extractValues(v1, v2, v1span, t.v2.span, t.sup.d111, t.liq.d111, t.gas.d111, false,false);
         cfuns.calcD222 = @(v1, v2) extractValues(v1, v2, t.v1.span, v2span, t.sup.d222, t.liq.d222, t.gas.d222, false,false);

         %% Third order mixed partial derivatives
         v1span = shrink_span(t.v1.span, 2 * t.v1.stepsize);
         v2span = shrink_span(t.v2.span, 1 * t.v2.stepsize);
         cfuns.calcD112 = @(v1, v2) extractValues(v1, v2, v1span, v2span, t.sup.d112, t.liq.d112, t.gas.d112, false, false);
         v1span = shrink_span(t.v1.span, 1 * t.v1.stepsize);
         v2span = shrink_span(t.v2.span, 2 * t.v2.stepsize);
         cfuns.calcD122 = @(v1, v2) extractValues(v1, v2, v1span, v2span, t.sup.d122, t.liq.d122, t.gas.d122, false, false);

      end
   end

   % ----------------------------------------------------------------------------
   function res = extractValues(v1, v2, span_v1, span_v2, supsamples, liqsamples, gassamples, const_v1, const_v2)

      if opt.assert_in_range
         assert_in_range(v1, span_v1);
         assert_in_range(v2, span_v2);
      else
         [v1, v2] = truncate_vectorized(v1, v2, span_v1, span_v2, opt.nan_outside_range);
      end

      phase    = phase_of(v1, v2, opt.phase_boundary);
      sph      = (phase == 0);
      lph      = (phase == 1);
      gph      = (phase == 2);

      res      = nan * ones(numel(v1), 1);
      res(sph) = extract_val_vectorized(supsamples, span_v1, span_v2, v1(sph), v2(sph), const_v1, const_v2);
      res(lph) = extract_val_vectorized(liqsamples, span_v1, span_v2, v1(lph), v2(lph), const_v1, const_v2);
      res(gph) = extract_val_vectorized(gassamples, span_v1, span_v2, v1(gph), v2(gph), const_v1, const_v2);

   end

end

% ============================================================================

function [v1, v2, vals] = dataObjects(m)
    v1 = m.v1;
    v2 = m.v2;
    vals = size(m.vals);
end
% ----------------------------------------------------------------------------
function phase = phase_of(v1, v2, phase_boundary)
% Flags are:
% * 0 designate "supercritical"
% * 1 designate "liquid"
% * 2 designate "gas"
    if isempty(phase_boundary)
        % We do not consider the sharp boundary, and can therefore make a
        % shortcut by considering everything as a single phase
        phase = zeros(numel(v1), 1); % everything considered "supercritical"
        return;
    end

    phase = nan * ones(numel(v1), 1);

    v1_c = phase_boundary{1}(1);
    v2_c = phase_boundary{1}(2);

    v1_super = (v1 >= v1_c);
    v2_super = (v2 >= v2_c);

    phase( v2_super &  v1_super) = 0; % hypercritical
    phase( v2_super & ~v1_super) = 2; % gas
    phase(~v2_super &  v1_super) = 1; % liquid

    % For the remaining entries, we are not beyond the critical value in either
    % direction.  We will have to do a more careful analysis to determine at
    % which side of the vapor-liquid boundary (v1, v2) is located.

    remains = find(isnan(phase) & ~isnan(v1) & ~isnan(v2));

    if numel(remains) > 0
        v2_bnd = v2(remains); % vector of temperature of the remaining cases
        v1_bnd = phase_boundary{2}(v2_bnd);
        phase_rem = ones(numel(remains), 1); % set to liquid by default
        phase_rem(v1(remains) < v1_bnd) = 2; % when pressure is lower than
                                             % vaporation pressure, the phase is gas
        phase(remains) = phase_rem;
    end
end

% ----------------------------------------------------------------------------

function phase = compute_phase(rows, cols, span_v1, span_v2, sharp_phase_boundary)

    v1_vals = linspace(span_v1(1), span_v1(2), rows);
    v2_vals = linspace(span_v2(1), span_v2(2), cols);
    [v2, v1] = meshgrid(v2_vals, v1_vals);
    phase = reshape(phase_of(v1(:), v2(:), sharp_phase_boundary), rows, cols);
end

% ----------------------------------------------------------------------------

function res = extract_val_vectorized(grid, span_v1, span_v2, v1, v2, const_v1, const_v2)
% Vectorized version of 'extract_val', where v1 and v2 are
% allowed to be vectors
   [p_ix, t_ix, p, t, nans] = ix_and_local_par(grid, span_v1, span_v2, v1, v2);
   res = NaN * ones(numel(v1), 1);
   lines = size(grid, 1);

   if const_v1  p = 0; end%#ok
   if const_v2  t = 0; end%#ok

   corners = [grid(p_ix   + lines * (t_ix - 1)), ...
              grid(p_ix+1 + lines * (t_ix - 1)), ...
              grid(p_ix   + lines * (t_ix    )), ...
              grid(p_ix+1 + lines * (t_ix    ))];

   weight = [(1-p).*(1-t), p.*(1-t), (1-p).*t, p.*t];

   res(~nans) = sum(corners .* weight, 2);
end

% ----------------------------------------------------------------------------

function [liq, gas] = prepare_by_phase(base, span_v1, span_v2, boundary)
% Prepare separate density matrices for the liquid and gas regions ('rliq',
%'rgas'). The returned matrices will have the same size as 'base',
% but will only contain numeric values for the samples corresponding to the
% liquid (or gas) phase, as well as a few extrapolated values beyond the
% discontinuity in order to faciltiate estimation of derivatives in this
% region.

% NB: This function still requires that the sharp boundary intersects the
% domain boundary at (v2 == v2_min) and (v1 > v1_min).  Fix when time. @@

    % This function should only be called if we use a sharp phase boundary,
    % hence the 'true' argument passed along to 'compute_phase' here.
    phase = compute_phase(size(base, 1), size(base, 2), span_v1, span_v2, boundary);

    liq = NaN * ones(size(base));
    gas = NaN * ones(size(base));

    lphase = (phase ~= 2); % everything that is not gas (i.e. including supercritical)
    gphase = (phase ~= 1); % everything that is not liquid (i.e. including supercritical)

    % Copying relevant values from the original sample matrix to the new,
    % phase-specific matrices
    liq(lphase(:)) = base(lphase(:));
    gas(gphase(:)) = base(gphase(:));

    % Determining largest column number below critical temperature (i.e. the
    % largest column for which extrapolation of values across discontinuity
    % is needed below)
    maxcol = max(find((cellfun(@(x)any(isnan(x)), num2cell(gas(end,:), 1)))));%#ok
    %dc = 1:maxcol+1;
    dc = 1:maxcol+3; %@ Strictly speaking, we should extrapolate along
                     %T-direction in order to get the additional values
                     %across boundary necessary for computing derivative.
                     %Currently, extrapolation is only in T.  Fix when time.

    % Determining row indices along discontinuity (column indices is given by (1:maxcol))
    dr = cellfun(@(x)min(find(~isnan(x)))-1, num2cell(liq(:, dc), 1));%#ok

    % Removing an extra row of values on each side of discontinuity, to
    % ensure that no value ended on the 'wrong side'
    drl = dr+1;
    drg = dr-1;
    liq(sub2ind(size(liq), drl, 1:maxcol+3)) = NaN;
    gas(sub2ind(size(gas), dr , 1:maxcol+3)) = NaN;

    % Extrapolating values along the boundary
    % @@ Perhaps extrapolation should also be based on temperature derivative?
    % @@ As of now, the extrapolation takes place along isotherms, using pressure derivatives
    % @@ only.
    ev = 5; % number of rows of extrapolated values
    for i = 0:ev-1
        % quadratic extrapolation of liquid density across boundary
        liq(sub2ind(size(liq),drl-i, dc)) = 3 * (liq(sub2ind(size(liq), drl-i+1, dc)) - ...
                                                liq(sub2ind(size(liq), drl-i+2, dc))) + ...
                                                liq(sub2ind(size(liq),drl-i+3, dc));

        % quadratic extrapolation of gas density across boundary
        gas(sub2ind(size(gas), drg+i+1, dc)) = 3 * (gas(sub2ind(size(gas), drg+i, dc)) - ...
                                                   gas(sub2ind(size(gas), drg+i-1, dc))) + ...
                                                   gas(sub2ind(size(gas), drg+i-2, dc));
    end
end

% ----------------------------------------------------------------------------

function [p_ix, t_ix, p, t, nans] = ix_and_local_par(grid, spanP, spanT, P, T)

    nans = isnan(P) | isnan(T);

    Psteps = size(grid, 1) - 1;
    Tsteps = size(grid, 2) - 1;

    dSpanP = diff(spanP);
    dSpanT = diff(spanT);
    Dp = dSpanP/Psteps;
    Dt = dSpanT/Tsteps;

    p_ix = 1 + floor(((P(~nans) - spanP(1)) * Psteps)/dSpanP);
    t_ix = 1 + floor(((T(~nans) - spanT(1)) * Tsteps)/dSpanT);

    p_ix = min(p_ix, Psteps);
    t_ix = min(t_ix, Tsteps);

    Psampled = spanP(1) + (p_ix - 1) * Dp;
    Tsampled = spanT(1) + (t_ix - 1) * Dt;

    p = (P(~nans) - Psampled) / Dp; % should be between 0 and 1
    t = (T(~nans) - Tsampled) / Dt; % should be between 0 and 1

end


% ----------------------------------------------------------------------------

function res = calcMulti(v1, v2, fun, fun_dv1, fun_dv2)
% function used as wrapper around the functions used to compute rho, enthalpy, and
% their derivatives.  The reason for this wrapper is to add support for the ADI
% framework when applicable.

    assert(isa(v2,'ADI') | isvector(v2));
    assert(isa(v1,'ADI') | isvector(v1));

    if(isa(v1,'ADI'))
        v1_val=v1.val;
    else
        v1_val=v1;
    end

    if (isa(v2,'ADI'))
        v2_val = v2.val;
    else
        v2_val = v2;
        if(numel(v2_val)==1)
           v2_val=ones(numel(v1_val),1)*v2;
        else
           v2_val = v2;
        end
    end

    res = fun(v1_val(:), v2_val(:));

    if(isa(v1,'ADI'))
        % Check if the function has an available derivative function
        dres_dp = fun_dv1(v1_val, v2_val);
        if(isa(v2,'ADI'))
            dres_dt = fun_dv2(v1_val, v2_val);
            assert(numel(v1.jac)==numel(v2.jac));

            %res = ADI(res, lMultDiag(dres_dp, v1.jac) + lMultDiag(dres_dt, v2.jac));
            tmp1 = lMultDiag(dres_dp, v1.jac);
            tmp2 = lMultDiag(dres_dt, v2.jac);
            res = ADI(res, cell_adder(tmp1, tmp2));%{tmp1{1} + tmp2{1}, tmp1{2} + tmp2{2}});
        else
            res = ADI(res, lMultDiag(dres_dp, v1.jac));
       end
    else
        if(isa(v2,'ADI'))
            dres_dt = fun_dv2(v1_val, v2_val);
            res = ADI(res, lMultDiag(dres_dt, v2.jac));
        end
    end
end



% ----------------------------------------------------------------------------
function noder(v1, v2)%#ok
% Function used as placeholder where no derivative function is available/implemented.
% Only used to throw an error
    error('The requested derivative is not implemented')
end

% ----------------------------------------------------------------------------

function tables = establish_data_tables(prop, vn, sharp_phase_boundary, const_derivatives)
% Based on the orignal table, generate additional, separate tables for each
% phase and for partial derivatives, if required

    v1_steplen = prop.(vn{1}).stepsize;
    v2_steplen = prop.(vn{2}).stepsize;

    tables.vals = prop.vals;
    tables.v1   = prop.(vn{1});
    tables.v2   = prop.(vn{2});

    % 'supercritical' table of derivatives (also the default table if we don't need separate
    % tables for each phase)
    tables.sup = deriv_tables(tables.vals, v1_steplen, v2_steplen, const_derivatives);

    if ~isempty(sharp_phase_boundary)
        % Establish separate tables for hypercritical, liquid and gas regions
        [liq, gas] = prepare_by_phase(tables.vals, tables.v1.span, tables.v2.span, sharp_phase_boundary);

        % compute table of derivatives
        tables.liq = deriv_tables(liq , v1_steplen, v2_steplen, const_derivatives);
        tables.gas = deriv_tables(gas , v1_steplen, v2_steplen, const_derivatives);
    else
        % This is the easy way.  A more comprehensive approach would ensure
        % that these tables were not needed at all if we do not enforce the
        % sharp boundary.
        tables.liq = tables.sup;
        tables.gas = tables.sup;
    end
end

% ----------------------------------------------------------------------------

function tables = deriv_tables(base, delta_v1, delta_v2, const_derivatives)

    % computing tables of first derivatives
    tables.d1 = diff(base, 1, 1) / delta_v1;
    tables.d2 = diff(base, 1, 2) / delta_v2;

    if const_derivatives
       % insert dummy row/col to keep derivative tables at same size as
       % original (the values in these rows will never be used)
       tables.d1 = [tables.d1; zeros(1, size(tables.d1, 2))];
       tables.d2 = [tables.d2, zeros(size(tables.d2, 1), 1)];
    else
       % derivatives not considered constant between each data point.
       % In this case, we also admit functions for higher derivatives.

       tables.d11  = diff(tables.d1,  1, 1) / delta_v1;
       tables.d111 = diff(tables.d11, 1, 1) / delta_v1;
       tables.d22  = diff(tables.d2,  1, 2) / delta_v2;
       tables.d222 = diff(tables.d22, 1, 2) / delta_v2;

       % cross derivatives
       tables.d12  = diff(tables.d1,  1, 2) / delta_v2;
       tables.d112 = diff(tables.d12, 1, 1) / delta_v1;
       tables.d122 = diff(tables.d12, 1, 2) / delta_v2;
    end
end


% ----------------------------------------------------------------------------

function res = get_varnames(prop)

   fnames = fieldnames(prop);

   % there should be three fields, one of which is named 'vals' (the sampled
   % table). The names of the other fields are the variable names)
   assert(numel(fnames) == 3);
   res = fnames(cellfun(@(x) ~strcmpi(x, 'vals'), fnames));

   assert(numel(res) == 2);
end

% ----------------------------------------------------------------------------
% @@ taken from ADI.m, where it is a private function
function J = lMultDiag(d, J1)
    n = numel(d);
    D = sparse((1:n)', (1:n)', d, n, n);
    J = cell(1, numel(J1));
    for k = 1:numel(J)
        J{k} = D*J1{k};
    end
end

% ----------------------------------------------------------------------------
function c = cell_adder(c1, c2)
    for i = 1:numel(c1)
        c{i} = c1{i} + c2{i};%#ok
    end
end

% ----------------------------------------------------------------------------

function assert_in_range(val, range)
   if (any(val < range(1)) || any (val > range(2)))
      error('Tried to extrapolate outside range of sampled table.');
   end
end

% ----------------------------------------------------------------------------

function [v1, v2] = truncate_vectorized(v1, v2, v1_span, v2_span, use_nan)

   if (use_nan)
      v1(v1 < v1_span(1) | v1 > v1_span(2)) = NaN;
      v2(v2 < v2_span(1) | v2 > v2_span(2)) = NaN;
   else
      v1_tr_1 = v1 < v1_span(1);
      v1_tr_2 = v1 > v1_span(2);
      v2_tr_1 = v2 < v2_span(1);
      v2_tr_2 = v2 > v2_span(2);
      if any(v1_tr_1 | v1_tr_2 | v2_tr_1 | v2_tr_2)
         % warning('There were values outside sampled table range.  Truncating.');
         v1(v1_tr_1) = v1_span(1);
         v1(v1_tr_2) = v1_span(2);
         v2(v2_tr_1) = v2_span(1);
         v2(v2_tr_2) = v2_span(2);
      end
   end
end
% ----------------------------------------------------------------------------
function span = shrink_span(span, ds)
    % shrink span by amount ds, centered
    span = [span(1) + ds/2, span(2) - ds/2];
end
