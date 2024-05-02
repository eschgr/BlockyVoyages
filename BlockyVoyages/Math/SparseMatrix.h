#if !defined(__SPARSE_MATRIX_H__)
#define __SPARSE_MATRIX_H__

#include "MathMain.h"

#include <vector>

namespace BlockyVoyages {
namespace Math {

struct SparseMatrixCell {
    int32 col;
    float32 val;
};

typedef std::vector<SparseMatrixCell> MatrixRow;

// A sparse matrix class which can be used to solve large systems of equations
// which produce a spare matrix. Not a first class math object, so everything
// is not filled out.
class SparseMatrix {
public:
    SparseMatrix(int32 n);

    void SetRow(int32 ind, const MatrixRow& row);
    bool ConjagateGradientSolve(VectorNf r, VectorNf& x);
private:
    std::vector<MatrixRow> m_rows;
    void Mul(const VectorNf& x, VectorNf& out);
};

}  // namespace BlockyVoyages
}  // namespace Math

#endif