function z = response_function(model)
    
    z = zeros(11);
%     x = [0.326072076 0.99 0.99 4.437159657 2.004183679 0.037399768 1.022057961];
    iValues = 1:0.1:6;
    for idx = 1:numel(iValues)
        i = iValues(idx);
        jValues = 1:0.1:6;
        parfor (jdx = 1:numel(jValues), 6)
%         for jdx = 1:numel(jValues)
            j = jValues(jdx);
            y = [0.326072076 0.99 0.99 i j 0.037399768 1.022057961];
            z(idx,jdx) = objectivefun_sync(y,model);
        end
    end
%     surf(4:0.1:5,1.5:0.1:2.5,z)
end