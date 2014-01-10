function clusterlist = find_clusters(varargin)
    clusterlist = cell(nargin,1);
    for i = 1:nargin
        clusterlist{i} = components(varargin{i});
    end
end