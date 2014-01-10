function test_mex_opmtransport
   clear all;
   !rm mex_opmtransport_implicit.mexa64

   gravity([0,1000,0])
   gravity on
   disp('init grid');
   G         = computeGeometry(cartGrid([150,150,1]));
   rock.perm = ones(G.cells.num, 1)*milli*darcy;
   rock.poro = ones(G.cells.num, 1);

   src    = [];
   src    = addSource(src,  1,           1e-3/day, 'sat', [1,0]);
   src    = addSource(src, G.cells.num, -1e-3/day, 'sat',[0,0]);
   disp('innerproducts')
   S      = computeMimeticIP(G, rock);
   resSol = initResSol(G, 0, 0);
   fluid  = initSimpleFluid('mu', [1,1], 'rho', [1000,10], 'n', [2,2]);
   disp('solving pressure');
   resSol = solveIncompFlow(resSol, G, S, fluid, 'src', src);

%[T, A] = computeTimeOfFlight(G, rock, resSol, 'src', src);
%full(A)

   sources = zeros(G.cells.num, 1);
   sources(src.cell)=src.rate;
   %%{
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   [mu, rho] = fluid.properties();
   g   = gravity() * (rho(1) - rho(2));
   dim = size(G.nodes.coords, 2);

   %harm          = zeros([G.faces.num, 1]);

   % nKg == n' * K * g *(rho1 - rho2) for all cellfaces.
   [K, r, c] = permTensor(rock, dim);

   assert (size(K,1) == G.cells.num, ...
          ['Permeability must be defined in active cells only.\n', ...
           'Got %d tensors, expected %d (== number of cells).'],   ...
          size(K,1), G.cells.num);

   nKg = sum(G.faces.normals(G.cells.faces(:,1), r) .* ...
      bsxfun(@times, K(cellNo,:), g(c)), 2);
   gflux = 2 ./ accumarray(G.cells.faces(:,1), 1 ./ nKg);
   gflux(any(G.faces.neighbors==0, 2))=0;
   %}
   mobo  = @(s) (1-s).*(1-s);
   dmobo = @(s) -2*(1-s);
   mobw  = @(s) s.*s;
   dmobw = @(s) 2*s;

   %mobo  = @(s) (1-s);
   %dmobo = @(s)-ones(size(s));
   %mobw  = @(s) s;
   %dmobw = @(s) ones(size(s));

   minSat = 0;     % Minimum (water) saturation
   maxSat = 1;     % Maximum (water) saturation
   capSat = @(sat) min(max(minSat, sat), maxSat);

   %g = gravity;
   s = zeros(G.cells.num, 1);
   disp('solving transport');
   dt = 150*year;
   for i=1:100,
      disp('explicit...');

      t=0;
      %      while t<dt,
      tic
      s0 = s;
      s  = s*0+0.5;
      res = inf;
      it = 0
      %while res>1e-6 & it<100,
      %   mob  = [ mobw(s),  mobo(s)];
      %   dmob = [dmobw(s), dmobo(s)];
         sigma = linspace(0,1,21)';tab = [mobw(sigma),mobo(sigma)];
         %ddt  = min(dt-t, estimate_dt(G, s, resSol.flux, sources, mob));
         [s, res]    = mex_opmtransport_implicit(G, s0, s, sigma(2)-sigma(1),0.0, int32(21), tab, resSol.flux, gflux, sources, dt);
      %   s = capSat(s);
      %   it = it+1;
      %   fprintf('norm %e\n', res);
      %end
      toc
      %t=t+dt;
      %return;
      %     end
      tic
      disp('implicit...');
      resSol = implicitTransport(resSol, G, dt, rock, fluid, 'src',src, 'verbose', true);
  toc
      subplot(1,2,1);
      cla;plotCellData(G, s);view(2);axis equal tight
      drawnow
      subplot(1,2,2);
      cla;plotCellData(G, resSol.s(:,1));view(2);axis equal tight
      drawnow
   end

end


function dt = estimate_dt(G, sat, flux, sources, mob)
   f = @(mob) mob(:,1)./sum(mob, 2);

   i = all(G.faces.neighbors > 0, 2);
   df    = f(mob(G.faces.neighbors(i, 2),:))-f(mob(G.faces.neighbors(i, 1),:));
   ds    = sat(G.faces.neighbors(i, 2)) - sat(G.faces.neighbors(i, 1));
   k     = ds~=0;
   df(k) = df(k)./ds(k);
   df(~k) = 0.0;
   df = df.*flux(i);

   fluxp = max(abs(df));

   s = sources > 0;
   df = sources(s).*(1.0 - f(mob(s, :)) )./(1.0 - sat(s));
   fluxp = max(fluxp, max(abs(df)));
   dt = 0.49*min(abs(G.cells.volumes/fluxp));

end

