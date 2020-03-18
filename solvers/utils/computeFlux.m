function [flux, wellsol] = computeFlux(G, u, T, W, mu, rho)
    dispif(mrstVerbose, 'computeFlux\n');

    flux = 0 * T{1};

    ind = all(G.faces.neighbors ~= 0, 2);
    c1 = G.faces.neighbors(ind, 1);
    c2 = G.faces.neighbors(ind, 2);
    flux(ind) = T{1}(ind) .* u(c1) - T{2}(ind) .* u(c2);

    c1 = max(G.faces.neighbors(~ind, :), [], 2);
    flux(~ind) = T{1}(~ind) .* u(c1) - T{2}(~ind);

    ind = G.faces.neighbors(:, 1) == 0;
    flux(ind) = -flux(ind);

    g = gravity();
    wellsol = repmat(struct('pressure', [], 'flux', []), [numel(W), 1]);
    for i = 1:numel(W)
        if strcmpi(W(i).type, 'bhp')
            pbh = W(i).val;
            dZ = W(i).dZ;
            wellsol(i).pressure = pbh + rho * g(3) * dZ;
            wellsol(i).flux = W(i).WI ./ mu .* (wellsol(i).pressure - u(W(i).cells));
        elseif strcmpi(W(i).type, 'rate')
            rate = W(i).val;
            % dZ = W(i).dZ;
            wellsol(i).pressure = u(G.cells.num+i);
            wellsol(i).flux = W(i).WI ./ mu .* (wellsol(i).pressure - u(W(i).cells) + rate);
        else
            error('code under development!');
            % write code here babbabababaababababababababababababababababababab
        end
    end
end
