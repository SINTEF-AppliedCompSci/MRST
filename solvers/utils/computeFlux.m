function flux = computeFlux(G, u, T)
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
end
