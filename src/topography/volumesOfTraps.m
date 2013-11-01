function v = volumesOfTraps(G, res, varargin)
% Find volume of a subset of traps
if nargin == 2
    trap = 1:max(res.traps);
else
    trap = varargin{1};
end
v = zeros(1,numel(trap));
for i = 1:numel(trap)
    ind = res.traps == trap(i);
    z = G.cells.z(ind);
%     h = max(G.cells.z) - z;
    fill_z = res.trap_z(trap(i));
    H=G.cells.H(ind);
    assert(all((fill_z - z)>=0));
    h_plume=min((fill_z - z),H);    
    v(i) = sum(max(eps, G.cells.volumes(ind).*h_plume));
end
end
