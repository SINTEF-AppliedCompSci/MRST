
function num_W_rate = getNoRateWells(W)
    num_W_rate = 0;
    for i=1:numel(W)
        if(strcmpi(W(i).type,'rate'))
            num_W_rate = num_W_rate+1;
        end
    end
end
