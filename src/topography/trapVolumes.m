function v = trapVolumes(G, res, trap)
if nargin < 3
    trap = 1:max(res.traps);
end
% Find volume of a subset of traps
    v = zeros(1,numel(trap));
    for i = 1:numel(trap)
        ind = res.traps == trap(i);
        H=G.cells.H(ind);
        z = G.cells.z(ind);        
        fill_z = res.trap_z(trap(i));
        assert(all((fill_z - z)>=0));
        h_plume=min((fill_z - z),H);
        v(i) = sum(max(eps, G.cells.volumes(ind).*h_plume));
    end
end