function [state] = FlowNTPFA(G, bc, fluid, W, OSflux, u0, tol, maxiter, varargin)

opt = struct('MatrixOutput', false, ...
             'src', []);
opt = merge_options(opt, varargin{:});

dispif(mrstVerbose, 'FlowNTPFA\n');
[mu, rho] = fluid.properties();
mu = mu(1);
rho = rho(1);

% Fix pressure if there are only hom Neumann BC (i.e. no pressure bc)
% and no bhp wells
num_W_rate = getNoRateWells(W);
is_bhp_wells = numel(W) > num_W_rate;
is_pressure_bc = any(strcmpi(bc.type, 'pressure'));
well_posed = is_pressure_bc || is_bhp_wells;
dispif(mrstVerbose, ['FlowNTPFA problem is well posed ', num2str(well_posed), '\n']);

% Expand u0 if there are rate wells
u0 = [u0; ones(num_W_rate, 1)];
T = TransNTPFA(G, bc, mu, u0, OSflux);
[A, b] = AssemAb(G, bc, mu, rho, T, W, opt.src, well_posed);

if opt.MatrixOutput
    state.A = A;
    state.jac = jacobian(G, bc, mu, rho, u0, OSflux, W, opt, well_posed);
    %return;
end

iter = 0;
res = zeros(maxiter+1, 1);
bnorm = norm(b, inf);
res(1) = norm(A*u0-b, inf) / bnorm;
%res(1)=tol+1;
u = u0;
while (res(iter+1) > tol && iter < maxiter)
    dispif(mrstVerbose & mod(iter, 1)==0, ['iter=', num2str(iter), ' res=', num2str(res(iter+1)), '\n'])
    u = A \ b;
    T = TransNTPFA(G, bc, mu, u, OSflux);
    [A, b] = AssemAb(G, bc, mu, rho, T, W, opt.src, well_posed);
    iter = iter + 1;
    res(iter+1) = norm(A*u-b, inf) / bnorm;
end
dispif(mrstVerbose, ['iter=', num2str(iter), ' res=', num2str(res(iter+1)), '\n'])
if iter == maxiter
    warning(['Maximum number of iterations ', num2str(maxiter), ...
        ' reached. Residual=', num2str(res(end))]);
end

[flux, wellsol] = computeFlux(G, u, T, W, mu, rho);
state.pressure = u(1:G.cells.num);
state.flux = flux;
state.wellSol = wellsol;
state.iter = iter;
state.res = res(1:iter+1);

% if opt.MatrixOutput
%     state.A = A;
%     state.jac = jacobian(G, bc, mu, rho, u0, OSflux, opt, well_posed);
% end

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function jac = jacobian(G, bc, mu, rho, u0, OSflux, W, opt, well_posed)
dispif(mrstVerbose, 'jacobian\n');

u0 = initVariablesADI(u0);
T = TransNTPFA(G, bc, mu, u0, OSflux);

[A, b, I, J, V] = AssemAb(G, bc, mu, rho, T, W, opt.src, well_posed);

% Each entry i in Au is A(i,:)*u
Au = cell(numel(u0.val), 1);
for k = 1:numel(u0.val)
    ii = find(I == k);
    jj = J(ii);
    Avals = V(ii);
    uvals = u0(jj);
    Au{k} = sum(Avals.*uvals);
end

eqs = combineEquations(Au{:});
jac = eqs.jac{1};

% figure
% spy(A)
% title('A')
% figure
% spy(jac)
% title('jacobian')
% figure
% r = symrcm(jac);
% spy(jac(r,r))
% title('symrcm jacobian')

% keyboard

% solve
%x = -eqs.jac{1} \ eqs.val;
end

