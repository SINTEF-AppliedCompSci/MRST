function ok = structPropEvaluated(s, name)
    ok = ~isempty(s.(name));
end