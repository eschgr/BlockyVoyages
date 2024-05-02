#include "stdafx.h"

#include "SparseMatrix.h"

namespace BlockyVoyages {
namespace Math {

SparseMatrix::SparseMatrix(int32 n)
    : m_rows(n)
{}

void SparseMatrix::SetRow(int32 ind, const MatrixRow& row) {
    m_rows[ind] = row;
}

bool SparseMatrix::ConjagateGradientSolve(VectorNf r, VectorNf& x) {
    // Calculate the Laplacian
    assert(m_rows.size() == r.n);

    int32 max_iters = m_rows.size() * 2;
    float32 tolerance = 1.0e-7f;

    int32 iter = 0;
    float32 rsqrd = r.MagnitudeSquared();
    float32 prev_rsqrd = rsqrd;
    VectorNf p(m_rows.size());

    // clear out x
    x.SetSize(r.n, 0.0f);
        
    while (rsqrd > tolerance && iter < max_iters) {
        iter++;
        if (iter == 1) {
            p = r;
        } else {
            float32 beta = rsqrd / prev_rsqrd;
            p *= beta;
            p += r;
        }
        VectorNf s;
        Mul(p, s);
        float32 alpha = rsqrd / (p.Dot(s));
        x += p * alpha;
        r -= s * alpha;
        prev_rsqrd = rsqrd;
        rsqrd = r.MagnitudeSquared();
    }
    return iter < max_iters;
}

void SparseMatrix::Mul(const VectorNf& x, VectorNf& out) {
    out.SetSize(x.n, 0.0f);
    for (uint32_t row = 0; row < m_rows.size(); ++row) {
        for (auto cell : m_rows[row]) {
            out[row] += cell.val * x[cell.col];
        }
    }
}

}  // namespace BlockyVoyages
}  // namespace Math