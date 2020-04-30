function W = addTrajectories(W, G, np)
% add well trajectory for plotting, currently just np points along a single
% quadratic curve

for k = 1:numel(W)
    c  = G.cells.centroids(W(k).cells,:);
    if ~isempty(c)
    if ~strcmp(W(k).type, 'aquifer')
        nc = size(c,1);
        c0 = c(1,:);
        % take point furthest away as toe
        v  = c-ones(nc,1)*c0;
        d2 = sum(v.*v, 2);
        [~, ix] = max(d2);
        c2 = c(ix,:);
        
        s = linspace(0, 1, np)';
        [f0, f1, f2] = deal((1-s).^2, 2*(1-s).*s, s.^2);
        % x,y,z - coords take second point as mean shifted twice away from c0-cn
        mc = mean(c);
        %
        n   = (c2-c0)/norm(c2-c0);
        pmc = c0 - dot(c0-mc, n)*n;
        c1  =  2*mc - pmc;
        c1(3) = min(c2(3),c1(3));
        traj  = f0*c0 + f1*c1 + f2*c2;
        
        
        %s = (linspace(0, 1, np).^(1/3))';
        %[f0, f1, f2] = deal((1-s).^2, 2*(1-s).*s, s.^2);
        %traj(:, 3) = f0*c0(3) + f1*c1(3) + f2*c2(3);
        
        W(k).trajectory = traj;
    else
        W(k).trajectory = c(1,:);
    end
    end
end
end
