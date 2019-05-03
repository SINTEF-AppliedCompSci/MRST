#ifndef OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED
#define OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED

namespace Opm{
  namespace Details{
    template <class DenseMatrix>
    DenseMatrix transposeDenseMatrix(const DenseMatrix& M)
    {
      DenseMatrix tmp;
      for (int i = 0; i < M.rows; ++i)
      for (int j = 0; j < M.cols; ++j)
	tmp[j][i] = M[i][j];
      
      return tmp;
    }
  }
  namespace Amg{
    template<class Matrix, class Vector>
    void getQuasiImpesWeights(const Matrix& matrix,const int pressureVarIndex, Vector& weights)
    {
      const Matrix& A = matrix;
      //Vector weights(rhs_->size());
      typedef typename Vector::block_type BlockVector;
      //typedef typename Matrix::MatrixBlock MatrixBlockType;
      typedef Dune::FieldMatrix<double, 3, 3> MatrixBlockType;
      BlockVector rhs(0.0);
      rhs[pressureVarIndex] = 1;
      const auto endi = A.end();
      for (auto i = A.begin(); i!=endi; ++i) {
	const auto endj = (*i).end();
	MatrixBlockType diag_block(0.0);
	for (auto j=(*i).begin(); j!=endj; ++j) {
	  if (i.index() == j.index()) {
	    diag_block = (*j);
	    break;
	  }
	}
	BlockVector bweights;
	auto diag_block_transpose = Opm::Details::transposeDenseMatrix(diag_block);
	diag_block_transpose.solve(bweights, rhs);
	double abs_max =
	  *std::max_element(
			    bweights.begin(),
			    bweights.end(), [](double a, double b){ return std::abs(a) < std::abs(b); } );
	bweights /= std::abs(abs_max);
	weights[i.index()] = bweights;
      }
      //return weights;
    }
  }
}
#endif
