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
    v(i) = sum(max(eps, G.cells.volumes(ind).*(fill_z - z)));
end
end
