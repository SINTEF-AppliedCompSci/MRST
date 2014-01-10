function deleteHandle(h)
% Delete one or more handles if they are handles. Helper function.
    arrayfun(@del, h);
end

function del(h)
    if ishandle(h)
        delete(h)
    end
end
