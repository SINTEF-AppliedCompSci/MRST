function groups = addFacilityGroup(groups, members, varargin)
    

    group = struct('name' , sprintf('G%d', numel(groups) + 1), ...
                   'ftype', {{}}, ...
                   'type', 'none', ...
                   'val' , nan, ...
                   'sign', 0, ...
                   'status', true, ...
                   'lims'  , []);
    group = merge_options(group, varargin{:});
    group.members = members;
    
    group = checkInput(group);
    
    groups = [groups; group];

end

function group = checkInput(group)

    if ~iscell(group.members), group.members = {group.members}; end
    nm = numel(group.members);
    if isempty(group.ftype), group.ftype = 'well'; end
    if ~iscell(group.ftype), group.ftype = {group.ftype}; end
    if numel(group.ftype) == 1
        group.ftype = repmat(group.ftype, nm, 1);
    end
    assert(numel(group.ftype) == nm);

end