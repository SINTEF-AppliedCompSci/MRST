function s = exactBL(model, state)

    G = model.G;
    fluid = model.fluid;

%     flux =  zeros(G.faces.num,1);
%     flux(model.operators.internalConn) = sum(state.flux,2);
    

%     flux  = sum(state.flux,2);
%     faces = G.cells.faces(:,1);
%     v = flux(faces).*(G.faces.normals(faces,1) > 0)./G.faces.areas(faces);
%     cells = rldecode((1:G.cells.num)', diff(G.cells.facePos),1);
%     v = accumarray(cells, v)/2;
%     v = v(2);

    vT = faceFlux2cellVelocity(G, sum(state.flux,2));
    v = vT(2,1);

    mobW = @(s) fluid.muW(1).*fluid.krW(s);
    mobN = @(s) fluid.muO(1).*fluid.krO(s);

    f = @(s) mobW(s)./(mobW(s) + mobN(1-s));
    
    s = linspace(0,1,1000);
    sAD = initVariablesADI(s);
    fAD = f(sAD);
    
    df = full(diag(fAD.jac{1}));
    df = @(sq) interp1(s, df, sq);
    
    ix  = s >= 0.5;
    dfinv = @(x) interp1(df(s(ix)), s(ix), x);

    fun = @(s) df(s) - (f(s) - f(0))./(s-0);
    s = linspace(0.05,1,1000)';
    ix = abs(fun(s)) == min(abs(fun(s)));
    sn = s(ix);
    
    sigma = (f(sn)-f(0))./(sn-0)*v;
    
    s = @(x,t) blsol(x, t, sigma, @(x) dfinv(x/v));
    
end
    
function s = blsol(x, t, sigma, dfinv)

    s =   1.*(x<=0) ...
        + dfinv(x./t).*(x > 0 & x <= t*sigma) ... 
        + 0.*(x>t*sigma);
    
    s(isnan(s)) = 0;
end
    
