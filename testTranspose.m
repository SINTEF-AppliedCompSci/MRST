clear all

dim = 3;
coltbl.coldim = (1 : dim)';
coltbl.num = dim;
rowtbl = coltbl;
rowtbl = replacefield(rowtbl, {'coldim', 'rowdim'});

colrowtbl = crossTable(coltbl, rowtbl, {});

col2row2tbl.coldim2 = colrowtbl.coldim;
col2row2tbl.rowdim2 = colrowtbl.rowdim;
col2row2tbl.coldim1 = colrowtbl.rowdim;
col2row2tbl.rowdim1 = colrowtbl.coldim;
col2row2tbl.num = colrowtbl.num;

prod = TensorProd();
prod.tbl1 = col2row2tbl;
prod.tbl2 = colrowtbl;
prod.replacefds1 = {{'coldim1', 'coldim'}, ...
                    {'rowdim1', 'rowdim'}};
prod.replacefds2 = {{'coldim', 'coldim2'}, ...
                    {'rowdim', 'rowdim2'}};
prod.reducefds = {'coldim2', 'rowdim2'};
prod.prodtbl = colrowtbl;
prod = prod.setup();

trans_T = SparseTensor();
trans_T = trans_T.setFromTensorProd(ones(col2row2tbl.num, 1), prod);

A = (1 : (dim*dim))';
prod = TensorProd();
prod.tbl1 = colrowtbl;
prod.tbl2 = coltbl;
prod.reducefds = {'coldim'};
prod.prodtbl = rowtbl;
prod = prod.setup();

A = setupTensor(A, prod);
B = A.transpose;

printTensor(A);
printTensor(B);