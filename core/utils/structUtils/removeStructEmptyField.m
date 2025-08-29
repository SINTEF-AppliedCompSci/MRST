function mstruct = removeStructEmptyField(mstruct)

    fds = fieldnames(mstruct);
    for ifd = 1 : numel(fds)
        fd = fds{ifd};
        if isempty(mstruct.(fd))
            mstruct = rmfield(mstruct, fd);
        elseif isstruct(mstruct.(fd))
            mstruct.(fd) = removeStructEmptyField(mstruct.(fd));
        else
            % do nothing
        end
    end
    
end

