function A = setupTensor(u, prod)
    A = SparseTensor();
    A = A.setFromTensorProd(u, prod);
end
