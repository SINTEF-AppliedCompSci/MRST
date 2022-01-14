function opt = mergeStructs(optA, optB)
    
    fnA  = fieldnames(optA);
    valA = struct2cell(optA);
    fnB  = fieldnames(optB);
    valB = struct2cell(optB);
    
    keep = ~ismember(fnB, fnA);
    fnB  = fnB(keep); valB = valB(keep);
    
    opt = cell2struct(vertcat(valA, valB), vertcat(fnA, fnB));
    
end