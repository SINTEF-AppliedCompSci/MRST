function fluid = initSimpleADIFluid(varargin)
   opt = struct('mu', [1 1 1], 'rho', [1 1 1], 'n', [1 1 1]);
   opt = merge_options(opt, varargin{:});

   krW = @(sw, varargin) sw.^opt.n(1);
   krO = @(so, varargin) so.^opt.n(2);
   krG = @(sg, varargin) sg.^opt.n(3);
   relperms = {krW, krO, krG};

   fluid.relPerm = @(sw, sg, varargin) relPerm(krW, krO, krG, sw, sg, varargin{:});


   names = {'W', 'O', 'G'};
   for i = 1:numel(names)
       n = names{i};
       bf = @(p, varargin) constantUnitBfactor(p, varargin{:});

       fluid.(['rho', n]) = opt.rho(i);
       fluid.(['rho', n, 'S']) = opt.rho(i);
       fluid.(['b', n]) = bf;
       fluid.(['B', n]) = bf;
       fluid.(['mu', n]) = @(p, varargin) constantViscosity(opt.mu(i), p, varargin{:});
       fluid.(['kr', n]) = relperms{i};
       fluid.(['krO', n]) = relperms{i};
   end
   fluid.rsSat = @(varargin) varargin{1}*0;
end

function [krW, krO, krG] = relPerm(krW, krO, krG, sw, sg, varargin)
    krW = krW(sw, varargin{:});
    krO = krO(1 - sw - sg, varargin{:});
    krG = krG(sg, varargin{:});
end

function B = constantUnitBfactor(p, varargin)
    B = p*0 + 1;
end

function mu = constantViscosity(mu, p, varargin)
    mu = p*0 + mu;
end
