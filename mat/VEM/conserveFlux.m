function flux = conserveFlux(state, G, rock, fluid,  bc, src)

f = G.cells.faces(:,1);
ncf = diff(G.cells.facePos);
fSgn = 1 - 2*(G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', ncf,1));

[ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
P = sparse(ii, jj, 1);
if ~isempty(src)
    rhs = sparse(src.cell, ones(numel(src.cell),1), src.rate, G.cells.num, 1);
else
    rhs = zeros(G.cells.num,1);
end
r = rhs-P*(state.flux(f).*fSgn);

if norm(r)/norm(state.flux) > 1e-15

    K = permTensor(rock, G.griddim)';
    
    [ii, jj] = blockDiagIndex(G.griddim*ones(G.cells.num,1));
    K = sparse(ii, jj, K(:));
    
    fn = bsxfun(@times,G.faces.normals(f,:), fSgn./G.faces.areas(f))';
    [ii, jj] = blockDiagIndex(G.griddim*ones(G.cells.num,1), ncf);
    fnMat = sparse(ii,jj, fn(:));
    delta = fnMat'*K*fnMat;
    delta = spdiags(delta, 0);

    
%     delta = (fn*K)';
%     [ii, jj] = blockDiagIndex(3*ones(numel(f),1), ones(numel(f),1));
%     delta = sparse(ii,jj,delta(:));
% 
%     fn = fn';
%     fn = sparse(ii,jj,fn(:));
%     delta = fn'*delta;
%     delta = delta*ones(size(delta,1),1);
% 
    ii = f;
    jj = (1:numel(f))';
    omega = sparse(ii, jj,1)*delta;

    for i = 1:G.faces.num
        d = delta(f == i);
        omega(i) = omega(i)/(numel(d)*prod(d));
    end
    
%     omega = ones(G.faces.num,1);
%     ff = G.cells.faces(G.cells.facePos(src.cell):G.cells.facePos(src.cell+1)-1);
%     omega(ff)= (1e-3)*omega(ff);

    c = G.faces.neighbors(f,:);
    c(fSgn == -1,:) = c(fSgn == -1, 2:-1:1);
    c = c(:,2);
    ii = rldecode((1:G.cells.num)', ncf,1);
    nz = c~=0; ii = ii(nz); c = c(nz);
    B = sparse(ii, c, -1./omega(f(nz)).*G.faces.areas(f(nz)).*fSgn(nz));
    
    neu = false(G.faces.num, 1);
    neu(bc.face(strcmp(bc.type, 'flux'))) = true;
    
    I = ones(numel(f),1);
    I(neu(f)) = 0;
    
    [ii, jj] = blockDiagIndex(ones(G.cells.num,1), ncf);
    ca = sparse(ii,jj,I)*(1./omega(f).*G.faces.areas(f).*fSgn);

    B = B + spdiags(-ca, 0, G.cells.num, G.cells.num);

    beta = B\r;
    beta = rldecode(beta, ncf,1);
    beta(neu(f)) = 0;
%     I = sparse(f, 1:numel(f), 1);
%     beta = I*beta.*G.faces.areas;

%     jj = 1:numel(f);
    
%   I = sparse(f(~neu(f)), jj(~neu(f)), 1, G.faces.num, numel(f));
    
    beta = sparse(f,1:numel(f),1)*beta.*G.faces.areas;
%     beta = I*beta.*G.faces.areas;
%     beta = beta(beta ~= 0).*G.faces.areas(~neu);
%     
    
%     isNeu = strcmp(bc.type, 'flux');
%     ii = true(G.faces.num, 1);
%     ii(bc.face(isNeu)) = false;
%     
%     flux = state.flux;
    
%     flux(~neu) = state.flux(~neu) - 1./omega(~neu).*beta(~neu);
    flux = state.flux - 1./omega.*beta;
%     flux(~neu) = state.flux(~neu) - 1./omega(~neu).*beta;
    
%     flux = state.flux - 1./omega.*beta;
    
    r = rhs-P*(flux(f).*fSgn);
    if norm(r)/norm(flux) > 1e-15;
        warning('Could not construct conservative flux field');
    end

else
    warning('Flux already conservative.');
    flux = state.flux;
end
    
end