function jsonstruct = removeJsonStructEmptyField(jsonstruct)

    fds = fieldnames(jsonstruct);
    for ifd = 1 : numel(fds)
        fd = fds{ifd};
        if isempty(jsonstruct.(fd))
            jsonstruct = rmfield(jsonstruct, fd);
        elseif isstruct(jsonstruct.(fd))
            jsonstruct.(fd) = removeJsonStructEmptyField(jsonstruct.(fd));
        else
            % do nothing
        end
    end
    
end

