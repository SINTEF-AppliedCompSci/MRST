function groups = addFacilityGroup(groups, members, varargin)
    

    group = struct('name' , sprintf('G%d', numel(groups) + 1), ...
                   'types', {{}}, ...
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
    if isempty(group.types), group.types = repmat({'well'}, nm, 1); end
        
    assert(numel(group.types) == nm);

end