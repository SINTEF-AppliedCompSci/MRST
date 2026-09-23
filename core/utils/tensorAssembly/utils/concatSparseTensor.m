function tensor = concatSparseTensor(tensor1, tensor2, varargin)

    opt = struct('checkUnique', true);
    opt = merge_options(opt, varargin{:});     

    tensor = tensor1;
    tensor.vals = [tensor1.vals; tensor2.vals];
    tensor.tbl = concatIndexArray(tensor1.tbl, tensor2.tbl, [], 'checkUnique', opt.checkUnique);
    
end

