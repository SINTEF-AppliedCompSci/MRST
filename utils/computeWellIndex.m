function WI = computeWellIndex(G, rock, radius, cells, varargin)
nc = numel(cells);
opt = struct('Dir', 'z', ...
             'Subset', [], ...
             'cellDims', [], ...
             'InnerProduct', 'ip_tpf', ...
             'Skin', zeros(nc, 1), ...
             'Kh', repmat(-1, nc, 1));
opt = merge_options(opt, varargin{:});

if ~isempty(opt.cellDims)
    [dx, dy, dz] = deal(opt.cellDims(cells, 1), ...
                        opt.cellDims(cells, 2), ...
                        opt.cellDims(cells, 3));
else
    if(isfield(G,'nodes'))
       [dx, dy, dz] = cellDims(G, cells);
    else
       [dx, dy, dz] = cellDimsCG(G, cells);
    end
end
if G.griddim > 2
   k = permDiag3D(rock, cells);
else
   k = permDiag2D(rock, cells);
   kz = 1./(1./k(:, 1) + 1./k(:, 2));
   k = [k, kz];
end
welldir = lower(opt.Dir);

if numel(welldir) == 1, welldir = welldir(ones([size(k,1), 1])); end
if numel(radius)  == 1, radius  = radius (ones([size(k,1), 1])); end

assert (numel(welldir) == size(k,1));
assert (numel(radius)  == size(k,1));

[d1, d2, ell, k1, k2] = deal(zeros([size(k,1), 1]));

ci = welldir == 'x';
[d1(ci), d2(ci), ell(ci), k1(ci), k2(ci)] = ...
   deal(dy(ci), dz(ci), dx(ci), k(ci,2), k(ci,3));

ci = welldir == 'y';
[d1(ci), d2(ci), ell(ci), k1(ci), k2(ci)] = ...
   deal(dx(ci), dz(ci), dy(ci), k(ci,1), k(ci,3));

ci = welldir == 'z';
[d1(ci), d2(ci), ell(ci), k1(ci), k2(ci)] = ...
   deal(dx(ci), dy(ci), dz(ci), k(ci,1), k(ci,2));

% Table look-up (interpolation) for mimetic or 0.14 for tpf
wc  = wellConstant(d1, d2, opt.InnerProduct);

re1 = 2 * wc .* sqrt((d1.^2).*sqrt(k2 ./ k1) + ...
                     (d2.^2).*sqrt(k1 ./ k2));
re2 = (k2 ./ k1).^(1/4) + (k1 ./ k2).^(1/4);

re  = reshape(re1 ./ re2, [], 1);
ke  = sqrt(k1 .* k2);

Kh = reshape(opt.Kh, [], 1); i = Kh < 0;
if G.griddim > 2
   Kh(i) = ell(i) .* ke(i);
else
   Kh(i) =           ke(i);
end

WI = 2 * pi * Kh ./ (log(re ./ radius) + reshape(opt.Skin, [], 1));

if any(WI < 0)
   if any(re < radius)
      error(id('WellRadius'), ...
           ['Equivalent radius in well model smaller than well ', ...
            'radius causing negative well index'].');
   else
      error(id('SkinFactor'), ...
            'Large negative skin factor causing negative well index.');
   end
end

if ~isempty(opt.Subset)
    % Only return calculated WI for requested cells
    WI = WI(opt.Subset);
end
end
%--------------------------------------------------------------------------

function wellConst = wellConstant(d1, d2, innerProd)
% table= [ratio mixedWellConstant]
table = [ 1, 0.292; ...
          2, 0.278; ...
          3, 0.262; ...
          4, 0.252; ...
          5, 0.244; ...
          8, 0.231; ...
          9, 0.229; ...
         16, 0.220; ...
         17, 0.219; ...
         32, 0.213; ...
         33, 0.213; ...
         64, 0.210; ...
         65, 0.210];

switch innerProd
   case {'ip_tpf', 'ip_quasitpf'}
      wellConst = 0.14;
   case {'ip_rt', 'ip_simple', 'ip_quasirt'}
      ratio = max(round(d1./d2), round(d2./d1));
      wellConst = interp1(table(:,1), table(:,2), ratio, ...
                          'linear', 'extrap');
    otherwise
      error(id('InnerProduct:Unknown'), ...
            'Unknown inner product ''%s''.', innerProd);
end
end

function p = permDiag3D(rock, inx)
if isempty(rock)
   error(id('Rock:Empty'), ...
         'Empty input argument ''rock'' is not supported');
elseif ~isfield(rock, 'perm')
   error(id('Rock:NoPerm'), ...
         '''rock'' must include permeability data');
elseif size(rock.perm, 2) == 1
   p = rock.perm(inx, [1, 1, 1]);
elseif size(rock.perm, 2) == 3
   p = rock.perm(inx, :);
else
   p = rock.perm(inx, [1, 4, 6]);
end
end

%--------------------------------------------------------------------------

function p = permDiag2D(rock, inx)
if isempty(rock)
   error(id('Rock:Empty'), ...
         'Empty input argument ''rock'' is not supported');
elseif ~isfield(rock, 'perm')
   error(id('Rock:NoPerm'), ...
         '''rock'' must include permeability data');
elseif size(rock.perm, 2) == 1
   p = rock.perm(inx, [1, 1]);
elseif size(rock.perm, 2) == 2
   p = rock.perm(inx, :);
else
   p = rock.perm(inx, [1, 3]);
end
end


function s = id(s)
    s = ['computeWellIndex:', s];
end