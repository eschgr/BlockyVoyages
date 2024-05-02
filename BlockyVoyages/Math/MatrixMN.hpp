#include <math.h>

namespace BlockyVoyages {
    // Constructors
    template<class T> inline MatrixMN<T>::MatrixMN(void) 
        : pntr(nullptr), m(0), n(0)
    {}

    template<class T> inline MatrixMN<T>::MatrixMN(int32 _m, int32 _n) 
        : pntr(nullptr)
    {
        SetSize(_m, _n);
    }

    template<class T> inline MatrixMN<T>::MatrixMN(T* v, int32 _m, int32 _n) 
        : pntr(nullptr)
    {
        Set(v, _m, _n);
    }

    template<class T> inline MatrixMN<T>::MatrixMN(T** v, int32 _m, int32 _n) 
        : pntr(nullptr)
    {
        Set(v, _m, _n);
    }

    template<class T> inline MatrixMN<T>::MatrixMN(const Matrix22<T>& other) 
        : pntr(nullptr)
    {
        Set(other);
    }

    template<class T> inline MatrixMN<T>::MatrixMN(const Matrix33<T>& other) 
        : pntr(nullptr)
    {
        Set(other);
    }

    template<class T> inline MatrixMN<T>::MatrixMN(const Matrix44<T>& other) 
        : pntr(nullptr)
    {
        Set(other);
    }

    template<class T> inline MatrixMN<T>::MatrixMN(const MatrixMN<T>& other) 
        : pntr(nullptr)
    {
        Set(other);
    }

    // Destructor
    template<class T> inline MatrixMN<T>::~MatrixMN(void) {
        delete[] pntr;
    }

    // Assignment methods
    template<class T> inline MatrixMN<T>& MatrixMN<T>::operator = (const MatrixMN<T>& other) {
        if (this != &other) {
            Set(other.pntr, other.m, other.n);
        }
        return *this;
    }

    template<class T> inline void MatrixMN<T>::SetSize(int32 _m, int32 _n) {
        delete[] pntr;

        if (_m == 0 || _n == 0) {
            throw "Invalid matrix size in function MatrixMN<T>::SetSize.";
        }

        m = _m;
        n = _n;
        pntr = new T[m * n];
    }

    template<class T> inline void MatrixMN<T>::Clear(void) {
        delete[] pntr;
        pntr = nullptr;
        m = 0;
        n = 0;
    }

    template<class T> inline void MatrixMN<T>::Set(T* v, int32 _m, int32 _n) {
        SetSize(_m, _n);
        for(int32 ind = 0; ind <m * n; ++ind) {
            pntr[ind] = v[ind];
        }
    }

    template<class T> inline void MatrixMN<T>::Set(T** v, int32 _m, int32 _n) {
        SetSize(_m, _n);
        for(int32 ind = 0; ind <m * n; ++ind) {
            pntr[ind] = v[ind];
        }
    }

    template<class T> inline void MatrixMN<T>::Set(const Matrix22<T>& other) {
        SetSize(2, 2);

        pntr[0] = other.m00;
        pntr[1] = other.m01;

        pntr[2] = other.m10;
        pntr[3] = other.m11;
    }

    template<class T> inline void MatrixMN<T>::Set(const Matrix33<T>& other) {
        SetSize(3, 3);

        pntr[0] = other.m00;
        pntr[1] = other.m01;
        pntr[2] = other.m02;
        pntr[3] = other.m10;
        pntr[4] = other.m11;
        pntr[5] = other.m12;
        pntr[6] = other.m20;
        pntr[7] = other.m21;
        pntr[8] = other.m22;
    }

    template<class T> inline void MatrixMN<T>::Set(const Matrix44<T>& other) {
        SetSize(4, 4);

        pntr[ 0] = other.m00;
        pntr[ 1] = other.m01;
        pntr[ 2] = other.m02;
        pntr[ 3] = other.m03;

        pntr[ 4] = other.m10;
        pntr[ 5] = other.m11;
        pntr[ 6] = other.m12;
        pntr[ 7] = other.m13;

        pntr[ 8] = other.m20;
        pntr[ 9] = other.m21;
        pntr[10] = other.m22;
        pntr[11] = other.m23;

        pntr[12] = other.m30;
        pntr[13] = other.m31;
        pntr[14] = other.m32;
        pntr[15] = other.m33;
    }

    template<class T> inline void MatrixMN<T>::Set(const MatrixMN<T>& other) {
        Set(other.pntr, other.m, other.n);
    }

    // General utility
    template<class T> inline void MatrixMN<T>::Clean(void) { // Cleans this matrix. Sets all close to zero values to zero
        for(int32 ind = 0; ind <m * n; ++ind) {
            if (IsCloseToZero(pntr[ind])) {
                pntr[ind] = 0;
            }
        }
    }

    template<class T> inline void MatrixMN<T>::Identity(void) {
        if (m != n) {
            throw "Invalid matrix size in function MatrixMN<T>::Identity.";
        }

        for(int32 row = 0, ind = 0; row <m; ++row) {
            for(int32 col = 0; col <row; ++col, ++ind) {
                pntr[ind] = 0;
            }
            pntr[ind++] = 1;
            for(int32 col = 0; col <n; ++col, ++ind) {
                pntr[ind] = 0;
            }
        }
    }

    template<class T> inline void MatrixMN<T>::Zero(void) {
        for(int32 ind = 0; ind <m * n; ++ind) {
            pntr[ind] = 0;
        }
    }

    // Gets the row specified.
    template<class T> inline const VectorN<T> MatrixMN<T>::GetRow(int32 row) const {
        return VectorN<T>(&pntr[row * n], n);
    }

    // Gets the column specified
    template<class T> inline const VectorN<T> MatrixMN<T>::GetColumn(int32 col) const {
        VectorN<T> colV(m);
        for(int32 matInd = col, vecInd = 0; vecInd < m; matInd += m, ++vecInd) {
            colV[vecInd] = pntr[matInd];
        }
        return colV;
    }

    template<class T> inline void MatrixMN<T>::SetRow(const VectorN<T>& v, int32 row) {
        if (v.n != n) {
            throw "Invalid vector size in function MatrixMN<T>::SetRow.";
        }

        for(int32 matInd = row * n, vecInd = 0; vecInd < n; ++matInd, ++vecInd) {
            pntr[matInd] = v[vecInd];
        }
    }

    template<class T> inline void MatrixMN<T>::SetColumn(const VectorN<T>& v, int32 col) {
        if (v.n != m) {
            throw "Invalid vector size in function MatrixMN<T>::SetRow.";
        }

        for(int32 matInd = col, vecInd = 0; vecInd <m; matInd += n, ++vecInd) {
            pntr[matInd] = v[vecInd];
        }
    }

    template<class T> inline bool MatrixMN<T>::InverseOf(const MatrixMN<T>& mat) {
        Set(mat);

        return Invert();
    }

    template<class T> inline bool MatrixMN<T>::Invert(void) {
        if (m != n) {
            throw "Invalid matrix dimensions in function MatrixMN<T>::Invert.";
        }

        int32* rowSwaps;         // which row have we swapped with the current one?
        rowSwaps = new int32[n];

        // do one pass for each diagonal element
        for(int32 pivot = 0; pivot <n; ++pivot) {
            int32 row, col;  // counters

            // find the largest magnitude element in the current column
            int32 maxRow = pivot;
            T maxElem = fabs(pntr[maxRow + n * pivot]);
            for (row = pivot + 1; row <n; ++row)
            {
                T elem = fabs(pntr[row * n + pivot]);
                if (elem> maxElem) {
                    maxElem = elem;
                    maxRow = row;
                }
            }

            // if max is zero, stop!
            if (IsZero(maxElem)) {
                delete[] rowSwaps;
                Clear();
                return false;
            }

            // if not in the current row, swap rows
            rowSwaps[pivot] = maxRow;
            if (maxRow != pivot) {
                // swap the row
                for(col = 0; col <n; ++col) {
                    T temp = pntr[maxRow * n + col];
                    pntr[maxRow * n + col] = pntr[pivot * n + col];
                    pntr[pivot * n + col] = temp;
                }
            }

            // multiply current row by 1/pivot to "set" pivot to 1
            T pivotRecip = 1.0f / pntr[n * pivot + pivot];
            for (col = 0; col <n; ++col) {
                pntr[pivot * n + col] *= pivotRecip;
            }

            // copy 1/pivot to pivot point (doing inverse in place)
            pntr[pivot * n + pivot] = pivotRecip;

            // now zero out pivot column in other rows 
            for (row = 0; row <n; ++row) {
                // don't subtract from pivot row
                if (row != pivot) {
                    // subtract multiple of pivot row from current row,
                    // such that pivot column element becomes 0
                    T factor = pntr[row * n + pivot];

                    // clear pivot column element (doing inverse in place)
                    // will end up setting this element to -factor*pivotInverse
                    pntr[row * n + pivot] = 0.0f;

                    // subtract multiple of row
                    for (col = 0; col <n; ++col) {
                        pntr[row * n + col] -= factor * pntr[pivot * n + col];    
                    }
                }
            }
        }

        // done, undo swaps in column direction, in reverse order
        int32 p = n;

        while (p> 0) {
            --p;
            // if row has been swapped
            if (rowSwaps[p] != p)
            {
                // swap the corresponding column
                for (int32 row = 0; row <n; ++row) {
                    T temp = pntr[row * n + rowSwaps[p]];
                    pntr[row * n + rowSwaps[p]] = pntr[row * n + p];
                    pntr[row * n + p] = temp;
                }
            }
        }

        delete[] rowSwaps;

        return true;
    }

    template<class T> inline void MatrixMN<T>::Transpose(void) {
        MatrixMN<T> trans;
        trans.TransposeOf(*this);
        Set(trans);
    }

    template<class T> inline void MatrixMN<T>::TransposeOf(const MatrixMN<T> &mat) {
        SetSize(mat.n, mat.m);
        for(int32 row = 0, matInd = 0; row <mat.m; ++row) {
            for(int32 col = 0; col <mat.n; ++col, ++matInd) {
                pntr[col * n + row] = mat(row, col);
            }
        }
    }

    template<class T> inline T MatrixMN<T>::Determinant(void) const {
        if (m != n) {
            return 0.0;
        }
        // Use LU Decomposition to decompose matrix
        MatrixMN<T> mat(*this);
        int32* rowPerms = new int32[n];
        T det;
        if (!mat.LUDecompose(rowPerms, det)) {
            delete[] rowPerms;
            return 0.0f;
        }

        // use the product of the diagonal elements of the LUDecomped matrix
        for(int32 i = 0; i <n; ++i) {
            det *= mat.pntr[i * n + i];
        }

        delete[] rowPerms;
        return det;
    }

    template<class T> inline T MatrixMN<T>::Trace(void) const {
        if (m != n || n == 0) {
            throw "Invalid matrix dimensions in function MatrixMN<T>::Trace.";
        }

        int32 sum = pntr[0];
        for(int32 ind = n + 1; ind <n * n; ind += n + 1) {
            sum += pntr[ind];
        }
        return sum;
    }

    template<class T> inline MatrixMN<T>& MatrixMN<T>::Tensor(const VectorN<T>& v1, const VectorN<T>& v2) {
        SetSize(v1.n, v2.n);
        for(int32 row = 0, ind = 0; row <m; ++row) {
            for(int32 col = 0; col <n; ++col, ++ind) {
                pntr[ind] = v1[row] * v2[col];
            }
        }
    }

    template<class T> inline const MatrixMN<T> MatrixMN<T>::operator *(const MatrixMN<T> &rhs) const {
        if (n != rhs.m) {
            throw "Invalid matrix dimensions in function MatrixMN<T>::operator *.";
        }

        MatrixMN<T> ans(m, rhs.n);
        for(int32 row = 0, ind = 0; row <m; ++row) {
            for(int32 col = 0; col <rhs.n; ++col, ++ind) {
                ans.pntr[ind] = pntr[row * n] * rhs.pntr[col];
                for(int32 val = 1; val <n; ++val) {
                    ans.pntr[ind] += pntr[row * n + val] * rhs.pntr[val * rhs.n + col];
                }
            }
        }

        return ans;
    }

    template<class T> inline const MatrixMN<T> MatrixMN<T>::operator *(	      T    rhs) const {
        return MatrixMN<T>(*this) *= rhs;
    }

    template<class T> inline const MatrixMN<T> MatrixMN<T>::operator +(const MatrixMN<T> &rhs) const {
        return MatrixMN<T>(*this) += rhs;
    }

    template<class T> inline const MatrixMN<T> MatrixMN<T>::operator -(const MatrixMN<T> &rhs) const {
        return MatrixMN<T>(*this) -= rhs;
    }

    template<class T> inline MatrixMN<T>& MatrixMN<T>::operator *=(const MatrixMN<T> &rhs) {
        return (*this) = (*this) * rhs;
    }

    template<class T> inline MatrixMN<T>& MatrixMN<T>::operator *=(      T    rhs) {
        for(int32 ind = 0; ind <m * n; ++ind) {
            pntr[ind] *= rhs;
        }

        return *this;
    }

    template<class T> inline MatrixMN<T>& MatrixMN<T>::operator +=(const MatrixMN<T> &rhs) {
        if (m != rhs.m && n != rhs.n) {
            throw "Invalid matrix dimensions in function MatrixMN<T>::operator +=.";
        }

        for(int32 ind = 0; ind <m * n; ++ind) {
            pntr[ind] += rhs.pntr[ind];
        }

        return *this;
    }

    template<class T> inline MatrixMN<T>& MatrixMN<T>::operator -=(const MatrixMN<T> &rhs) {
        if (m != rhs.m && n != rhs.n) {
            throw "Invalid matrix dimensions in function MatrixMN<T>::operator +=.";
        }

        for(int32 ind = 0; ind <m * n; ++ind) {
            pntr[ind] -= rhs.pntr[ind];
        }

        return *this;
    }

    // Methods used to replace the = sign. Use as ans = m1 op m2
    template<class T> inline MatrixMN<T>& MatrixMN<T>::Add(const MatrixMN<T>& m1, const MatrixMN<T>& m2) {
        if (m1.m != m2.m && m1.n != m2.n) {
            throw "Invalid matrix dimensions in function MatrixMN<T>::operator +=.";
        }

        SetSize(m1.m, m1.n);

        for(int32 ind = 0; ind <m * n; ++ind) {
            pntr[ind] = m1.pntr[ind] + m2.pntr[ind];
        }

        return *this;
    }

    template<class T> inline MatrixMN<T>& MatrixMN<T>::Sub(const MatrixMN<T>& m1, const MatrixMN<T>& m2) {
        if (m1.m != m2.m && m1.n != m2.n) {
            throw "Invalid matrix dimensions in function MatrixMN<T>::operator +=.";
        }

        SetSize(m1.m, m1.n);

        for(int32 ind = 0; ind <m * n; ++ind) {
            pntr[ind] = m1.pntr[ind] - m2.pntr[ind];
        }

        return *this;
    }

    template<class T> inline MatrixMN<T>& MatrixMN<T>::Mul(const MatrixMN<T>& m1, const MatrixMN<T>& m2) {
        if (m1.n != m2.m) {
            throw "Invalid matrix dimensions in function MatrixMN<T>::operator *.";
        }

        SetSize(m1.m, m2.n);
        for(int32 row = 0, ind = 0; row <m1.m; ++row) {
            for(int32 col = 0; col <m2.n; ++col, ++ind) {
                pntr[ind] = m1(row, 0) * m2(0, col);
                for(int32 val = 1; val <n; ++val) {
                    pntr[ind] += m1(row, val) * m2(val, col);
                }
            }
        }

        return *this;
    }

    template<class T> inline MatrixMN<T>& MatrixMN<T>::Mul(const MatrixMN<T>& m1,         T    f)  {
        SetSize(m1.m, m1.n);

        for(int32 ind = 0; ind <m * n; ++ind) {
            pntr[ind] = m1.pntr[ind] * f;
        }

        return *this;
    }

    template<class T> inline const VectorN<T> MatrixMN<T>::operator * (const VectorN<T> &rhs) const {
        if (n != rhs.n) {
            throw "Invalid vector dimensions in function MatrixMN<T>::operator *.";
        }

        VectorN<T> ans(m);
        SetSize(m1.m, m2.n);
        for(int32 row = 0, ind = 0; row <m1.m; ++row) {
            ans[row] = pntr[ind] * rhs[0];
            for(int32 col = 1; col <rhs.n; ++col, ++ind) {
                ans[row] += pntr[ind] * rhs[col];
            }
        }

        return ans;
    }

    template<class T> inline bool MatrixMN<T>::operator==(const MatrixMN<T>& other) {
        if (this == &other) {
            return true;
        }

        for(int32 ind = 0; ind <m * n; ++ind) {
            if (pntr[ind] != other.pntr[ind]) {
                return false;
            }
        }

        return true;
    }

    template<class T> inline bool MatrixMN<T>::operator!=(const MatrixMN<T>& other) {
        return !((*this) == other);
    }

    template<class T> inline bool MatrixMN<T>::GetSymEigenInfo(MatrixMN<T>& ev, VectorN<T>& ew) const {
        ev = *this;
        return ev.SymEigenDecompose(ew);
    }

    template<class T> inline bool MatrixMN<T>::SymEigenDecompose(VectorN<T>& ew) {
        throw "Can't use method MatrixMN<T>::SymEigenDecompose with the requested type.";
        return false;
    }

    template<> inline bool MatrixMN<float32>::SymEigenDecompose(VectorN<float32>& ew) {
        VectorN<float32> subd(n);
        ew.SetSize(n);
        Tridiagonalize(ew, subd);
        if (QLFactorize(ew, subd)) {
            // sort them
            SortEigenInfo(ew);
            return true;
        }
        return false;
    }

    template<> inline bool MatrixMN<float64>::SymEigenDecompose(VectorN<float64>& ew) {
        VectorN<float64> subd(n);
        ew.SetSize(n);
        Tridiagonalize(ew, subd);
        if (QLFactorize(ew, subd)) {
            // sort them
            SortEigenInfo(ew);
            return true;
        }
        return false;
    }
    // places the LU decomposiiton of mat into this
    // rowPerms: a vector used to how the rows are swapped during the
    //  decomposition
    // dir: +/-1 when succesful; +1 for even number of row swaps,
    //  -1 for odd number of row swaps
    template<class T> inline bool MatrixMN<T>::LUDecompositionOf(const MatrixMN<T>& mat, int32* rowPerms, T& dir) {
        // Copy over information into this to be able to use
        // another function. Slower, but easier to read and 
        // probably won't make a huge difference in efficiency
        Set(mat);

        // Perform the decomposition and check success
        return LUDecompose(rowPerms, dir);
    }

    // makes this matrix into its LU Decomposition
    // rowPerms: a vector used to how the rows are swapped during the
    //  decomposition
    // dir: +/-1 when succesful; +1 for even number of row swaps,
    //  -1 for odd number of row swaps
    template<class T> inline bool MatrixMN<T>::LUDecompose(int32* rowPerms, T& dir) {
        if (m != n) {
            throw "Invalid matrix dimensions in function MatrixMN<T>::InverseOf.";
            return false;
        }

        int32 i, iMax, j, k;
        T big, dum, sum, temp;
        VectorN<T> vv(n);

        dir = 1.0;

        for(i = 0; i <n; ++i) {
            big = 0.0;
            for(j = 0; j <n; ++j) {
                if ((temp = fabs(pntr[i * n + j]))> big) {
                    big = temp;
                }
            }
            if (big == 0.0) {
                Clear();
                return false;
            }
            vv[i] = 1.0f / big;
        }

        for(j = 0; j <n; ++j) {
            for(i = 0; i <j; ++i) {
                sum = pntr[i * n + j];
                for(k = 0; k <i; ++k) {
                    sum -= pntr[i * n + k] * pntr[k * n + j];
                }
                pntr[i * n + j] = sum;
            }

            big = 0.0;

            for (i = j; i <n; ++i) {
                sum = pntr[i * n + j];
                for(k = 0; k <j; ++k)
                    sum -= pntr[i * n + k] * pntr[k * n + j];
                pntr[i * n + j] = sum;

                if ((dum = vv[i] * fabs(sum))>= big) {
                    big = dum;
                    iMax = i;
                }
            }

            if (j != iMax) {
                for(k = 0; k <n; ++k) {
                    dum = pntr[iMax * n + k];
                    pntr[iMax * n + k] = pntr[j * n + k];
                    pntr[j * n + k] = dum;
                }
                dir *= -1;
                vv[iMax] = vv[j];
            }

            rowPerms[j] = iMax;
            if (pntr[j * n + j] == 0.0) {
                Clear();
                return false;
            }

            if (j != n - 1) {
                dum = 1.0f / pntr[j * n + j];
                for(i = j + 1; i <n; ++i) {
                    pntr[i * n + j] *= dum;
                }
            }
        }

        return true;
    }

    // Uses LU Backsubstitution to solve for a linear set of equations
    // Assumes this is a LU Decompoosed matrix
    // b: vector on the right side of the equations
    // x: vector to store the solutions
    // rowPerms is the row swap vector from the LU Decomposition
    // x is reset if no solution exists
    template<class T> inline bool MatrixMN<T>::LUBacksubstitute(const VectorN<T>& b, VectorN<T>& x, const int32* rowPerms) const {
        if (m != n) {
            x.Clear();
            throw "Invalid matrix dimensions in function MatrixMN<T>::InverseOf.";
            return false;
        }

        int32 i, ii = -1, ip, j;
        T sum;

        x = b;

        for(i = 0; i <n; ++i) {
            ip = rowPerms[i];
            sum = x[ip];
            x[ip] = x[i];
            if (ii>= 0) {
                for(j = ii; j <= i - 1; ++j) {
                    sum -= pntr[i * n + j] * x[j];
                }
            }
            else if (IsNotZero(sum)) {
                ii = i;
            }
            x[i] = sum;
        }

        for(i = n - 1; i>= 0; --i) {
            sum = x[i];
            for(j = i + 1; j <n; ++j) {
                sum -= pntr[i * n + j] * x[j];
            }
            x[i] = sum / pntr[i * n + i];
        }

        return true;
    }

    // Solves a linear set of equations using LU Decomposition
    // Assums nothing about this matrix
    // b: vector on the right side of the equations
    // x: vector to store the solution if one exists
    // x is reset if no solution exists
    template<class T> inline bool MatrixMN<T>::LULinearSolver(const VectorN<T>& b, VectorN<T>& x) const {
        if (m != n) {
            x.Clear();
            throw "Invalid matrix dimensions in function MatrixMN<T>::InverseOf.";
            return false;
        }

        MatrixMN luMat(*this);
        T dir;
        int32* rowPerms = new int32[n];
        if (!luMat.LUDecompose(rowPerms, dir)) {
            x.Clear();
            delete[] rowPerms;
            return false;
        }

        luMat.LUBacksubstitute(b, x, rowPerms);

        delete[] rowPerms;
        return true;
    }

    template<class T> inline bool MatrixMN<T>::SVDecompose(VectorN<T>& sigma, MatrixMN<T>& v) {
        int flag, i, its, j, jj, k, l = 0, nm = 0;
        T anorm, c, f, g, h, s, scale, x, y, z;
        T* rv1 = new T[n];
        g = scale = anorm = 0.0;

        sigma.SetSize(n);
        v.SetSize(n, n);
        
        for(i = 0; i <n; ++i) {
            l = i + 1;
            rv1[i] = scale * g;
            g = s = scale = 0.0;
            if (i <m) {
                for(k = i; k <m; ++k) {
                    scale += fabs(pntr[k * n + i]);
                }

                if (IsNotZero(scale)) {
                    for(k = i; k <m; ++k) {
                        int32 ind = k * n + i;
                        pntr[ind] /= scale;
                        s += pntr[ind] * pntr[ind];
                    }
                    f = pntr[i * (1 + n)];
                    g = -sign(sqrt(s), f);
                    h = f * g - s;
                    pntr[i * (1 + n)]=f-g;
                    for (j = l; j <n; ++j) {
                        for(s = 0.0, k = i; k <m; ++k) {
                            s += pntr[k * n + i] * pntr[k * n + j];
                        }
                        f = s / h;
                        for(k = i; k <m; ++k) {
                            pntr[k * n + j] += f * pntr[k * n + i];
                        }
                    }
                    for(k = i; k <m; ++k) {
                        pntr[k * n + i] *= scale;
                    }
                }
            }
            sigma[i] = scale * g;
            g = s = scale = 0.0;
            if (i <m && i != n - 1) {
                for(k = l; k <n; ++k) {
                    scale += fabs(pntr[i * n + k]);
                }

                if (IsNotZero(scale)) {
                    for(k = l; k <n; ++k) {
                        pntr[i * n + k] /= scale;
                        s += pntr[i * n + k] * pntr[i * n + k];
                    }
                    f = pntr[i * n + l];
                    g = -sign(sqrt(s), f);
                    h = f * g - s;
                    pntr[i * n + l] = f - g;
                    for(k = l; k <n; ++k) 
                        rv1[k] = pntr[i * n + k] / h;
                    for (j = l; j <m; ++j) {
                        s = 0.0;
                        for(k = l; k <n; ++k) {
                            s += pntr[j * n + k] * pntr[i * n + k];
                        }

                        for(k = l; k <n; ++k) {
                            pntr[j * n + k] += s * rv1[k];
                        }
                    }
                    for(k = l; k <n; ++k) {
                        pntr[i * n + k] *= scale;
                    }
                }
            }
            anorm = std::max(anorm, fabs(sigma[i]) + fabs(rv1[i]));
        }

        for(i = n - 1; i>= 0; --i) {
            if (i <n - 1) {
                if (IsNotZero(g)) {
                    T invG = 1.0f / g;
                    for(j = l; j <n; ++j) {
                        v.pntr[j * n + i] = (pntr[i * n + j] / pntr[i * n + l]) * invG;
                    }

                    for(j = l; j <n; ++j) {
                        s = 0.0f;
                        for(k = l; k <n; ++k) {
                            s += pntr[i * n + k] * v.pntr[k * n + j];
                        }

                        for(k = l; k <n; ++k) {
                            v.pntr[k * n + j] += s * v.pntr[k * n + i];
                        }
                    }
                }
                for(j = l; j <n; ++j) {
                    v.pntr[i * n + j] = v.pntr[j * n + i]=0.0;
                }
            }

            v.pntr[i * (n + 1)] = 1.0;
            g = rv1[i];
            l = i;
        }

        for(i = std::min(m, n) - 1; i>= 0; --i) {
            l = i + 1;
            g = sigma[i];
            for(j = l; j <n; ++j) {
                pntr[i * n + j] = 0.0;
            }

            if (IsNotZero(g)) {
                g = 1.0f / g;
                for(j = l; j <n; ++j) {
                    s = 0.0f;
                    for(k = l; k <m; ++k) {
                        s += pntr[k * n + i] * pntr[k * n + j];
                    }
                    f = (s / pntr[i * n + i]) * g;
                    for(k = i; k <m; ++k) 
                        pntr[k * n + j] += f * pntr[k * n + i];
                }
                for(j = i; j <m; ++j) 
                    pntr[j * n + i] *= g;

            } 
            else {
                for(j = i; j <m; ++j) {
                    pntr[j * n + i]=0.0;
                }
            }
            ++pntr[i * n + i];
        }
        for(k = n - 1; k>= 0; --k) {
            for(its = 0; its <30; ++its) {
                flag = 1;
                for(l = k; l>= 0; --l) {
                    nm = l - 1;
                    if (IsZero(rv1[l])) {
                        flag=0;
                        break;
                    }
                    if (IsZero(sigma[nm])) {
                        break;
                    }
                }
                if (flag) {
                    c = 0.0f;
                    s = 1.0f;
                    for(i = l; i <= k; ++i) {
                        f=s*rv1[i];
                        rv1[i]=c*rv1[i];
                        if (IsZero(f)) {
                            break;
                        }

                        g = sigma[i];
                        h = Pythag(f, g);
                        sigma[i] = h;
                        h = 1.0f / h;
                        c = g * h;
                        s = -f * h;
                        for(j = 0; j <m; ++j) {
                            y=pntr[j * n + nm];
                            z=pntr[j * n + i];
                            pntr[j * n + nm] = y * c + z * s;
                            pntr[j * n + i] = z * c - y * s;
                        }
                    }
                }
                z=sigma[k];
                if (l == k) {
                    if (z <0.0) {
                        sigma[k] = -z;
                        for(j = 0; j <n; ++j) {
                            v.pntr[j * n + k] = -v.pntr[j * n + k];
                        }
                    }
                    break;
                }
                if (its == 30) {
                    return false;
                }

                x = sigma[l];
                nm = k - 1;
                y = sigma[nm];
                g = rv1[nm];
                h = rv1[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);
                g = Pythag(f, 1.0f);
                f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x;
                c = s = 1.0f;
                for(j = l; j <= nm; ++j) {
                    i = j + 1;
                    g = rv1[i];
                    y = sigma[i];
                    h = s * g;
                    g = c * g;
                    z = Pythag(f, h);
                    rv1[j] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y *= c;
                    for(jj = 0; jj <n; ++jj) {
                        x = v.pntr[jj * n + j];
                        z = v.pntr[jj * n + i];
                        v.pntr[jj * n + j] = x * c + z * s;
                        v.pntr[jj * n + i] = z * c - x * s;
                    }
                    z = Pythag(f, h);
                    sigma[j] = z;
                    if (IsNotZero(z)) {
                        z = 1.0f / z;
                        c = f * z;
                        s = h * z;
                    }
                    f = c * g + s * y;
                    x = c * y - s * g;
                    for(jj = 0; jj <m; ++jj) {
                        y = pntr[jj * n + j];
                        z = pntr[jj * n + i];
                        pntr[jj * n + j]= y * c + z * s;
                        pntr[jj * n + i]= z * c - y * s;
                    }
                }
                rv1[l] = 0.0f;
                rv1[k] = f;
                sigma[k] = x;
            }
        }
        delete[] rv1;
        return true;
    }

    template<class T> inline void MatrixMN<T>::SVBacksubstitute(const VectorN<T>& sigma, const MatrixMN<T>& v, const VectorN<T>& b, VectorN<T>& x) const {
        // Solves A*X = B for a vector X, where A is specifed by the arrays u[1..m][1..n], w[1..n],
        // v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
        // square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
        // No input quantities are destroyed, so the routine may be called sequentially with different b's.

        if (b.n != m) {
            throw "Vector and matrix dimensions do not match in MatrixMN<T>::SVBacksubstitute.";
        }

        x.SetSize(n);

        int jj, j, i;
        T s, *tmp;
        tmp = new T[n];
        for (j = 0; j <n; ++j) { //Calculate UTB.
            s = 0.0;
            if (sigma[j]) { //Nonzero result only if wj is nonzero.
                for(i = 0; i <m; ++i) {
                    s += pntr[i * n + j]*b[i];
                }
                s /= sigma[j]; //This is the divide by wj .
            }
            tmp[j] = s;
        }
        for(j = 0; j <n; ++j) { //Matrix multiply by V to get answer.
            s=0.0;
            for(jj = 0; jj <n; ++jj) {
                s += v.pntr[j * v.n + jj] * tmp[jj];
            }

            x[j]=s;
        }
        delete[] tmp;
    }

    template<class T> inline bool MatrixMN<T>::SVDLinearSolver(const VectorN<T>& b, VectorN<T>& x) const {
        MatrixMN<T> u(*this);
        VectorN<T> sigma(n);
        MatrixMN<T> v(n, n);

        if (!u.SVDecompose(sigma, v)) {
            x.Clear();
            return false;
        }

        // clean up the singular values we got out
        float64 wmax = 0.0; //Will be the maximum singular value obtained.
        float64 wmin = 0.0;

        for(int32 j = 0; j <3; ++j) 
            if (sigma[j]> wmax)
                wmax = sigma[j];

        /*This is where we set the threshold for singular values allowed to be nonzero. The constant
        is typical, but not universal. You have to experiment with your own application.*/
        wmin = wmax * 1.0e-6;
        for(j=0; j <3; ++j) 
            if (sigma[j] <wmin) 
                sigma[j] = 0.0;

        u.SVBacksubstitute(sigma, v, b, x); //Now we can backsubstitute.
        return true;
    }

    template<class T> inline void MatrixMN<T>::Tridiagonalize(VectorN<T>& diag, VectorN<T>& subd) {
        int32 i0, i1, i2, i3;

        T h, scale;
        for(i0 = n - 1, i3 = n - 2; i0>= 2; --i0, --i3)
        {
            h = 0.0f;
            scale = 0.0f;
            for(i2 = 0; i2 <= i3; ++i2)
                scale += fabs(pntr[i0 * n + i2]);

            if (IsZero(scale))
                subd[i0] = pntr[i0 * n + i3];
            else
            {
                T invScale = 1.0f / scale;
                for (i2 = 0; i2 <= i3; ++i2)
                {
                    pntr[i0 * n + i2] *= invScale;
                    h += pntr[i0 * n + i2] * pntr[i0 * n + i2];
                }

                T f = pntr[i0 * n + i3];
                T g = sqrt(h);
                if (f> 0.0f)
                    g = -g;

                subd[i0] = scale * g;
                h -= f * g;
                pntr[i0 * n + i3] = f - g;
                f = 0.0f;
                T invH = 1.0f / h;
                for (i1 = 0; i1 <= i3; ++i1)
                {
                    pntr[i1 * n + i0] = pntr[i0 * n + i1] * invH;
                    g = 0.0f;
                    for (i2 = 0; i2 <= i1; ++i2)
                        g += pntr[i1 * n + i2] * pntr[i0 * n + i2];

                    for (i2 = i1+1; i2 <= i3; ++i2)
                        g += pntr[i2 * n + i1] * pntr[i0 * n + i2];

                    subd[i1] = g * invH;
                    f += subd[i1] * pntr[i0 * n + i1];
                }

                T halfdivH = 0.5f * f * invH;
                for (i1 = 0; i1 <= i3; i1++)
                {
                    f = pntr[i0 * n + i1];
                    g = subd[i1] - halfdivH*f;
                    subd[i1] = g;
                    for (i2 = 0; i2 <= i1; ++i2)
                        pntr[i1 * n + i2] -= f * subd[i2] + g * pntr[i0 * n + i2];
                }
            }

            diag[i0] = h;
        }

        subd[1] = pntr[n];
        diag[1] = 0.0f;

        diag[0] = subd[0] = 0;
        for (i0 = 0, i3 = 0; i0 <n; ++i0, ++i3) {
            if (IsNotZero(diag[i0])) {
                for (i1 = 0; i1 <i3; ++i1) {
                    T sum = 0.0f;
                    for (i2 = 0; i2 <i3; ++i2)
                        sum += pntr[i0 * n + i2] * pntr[i2 * n + i1];

                    for (i2 = 0; i2 <i3; ++i2)
                        pntr[i2 * n + i1] -= sum * pntr[i2 * n + i0];
                }
            }

            diag[i0] = pntr[i0 * n + i0];
            pntr[i0 * n + i0] = 1;
            for (i1 = 0; i1 <i3; ++i1)
                pntr[i1 * n + i0] = pntr[i0 * n + i1] = 0.0f;
        }

        // re-ordering if MgcEigen::QLAlgorithm is used subsequently
        for (i0 = 1, i3 = 0; i0 <n; i0++, i3++)
            subd[i3] = subd[i0];

        subd[n-1] = 0;
    }

    template<class T> inline bool MatrixMN<T>::QLFactorize(VectorN<T>& diag, VectorN<T>& subd) {
        const int32 iMaxIter = 32;

		int32 i0, i1, i2, i3, i4;
		for (i0 = 0; i0 <n; i0++)
		{
			for (i1 = 0; i1 <iMaxIter; i1++)
			{
				i2 = i0;
				while(i2 <n-1 && IsNotZero(subd[i2]))
					++i2;

				if (i2 == i0)
					break;

				T g = (diag[i0 + 1] - diag[i0]) / (2.0f * subd[i0]);
				T r = sqrt(g * g + 1.0f);

				if (g <0.0f)
					g = diag[i2] - diag[i0] + subd[i0] / (g - r);
				else
					g = diag[i2] - diag[i0] + subd[i0] / (g + r);
				
				T fSin = 1.0f, fCos = 1.0f, p = 0.0f;
				for (i3 = i2; i3> i0; --i3)
				{
					T f = fSin * subd[i3 - 1];
					T b = fCos * subd[i3 - 1];
					if (fabs(f)>= fabs(g))
					{
						fCos = g / f;
						r = sqrt(fCos * fCos + 1.0f);
						subd[i3] = f * r;
						fSin = 1.0f / r;
						fCos *= fSin;
					}
					else
					{
						fSin = f / g;
						r = sqrt(fSin * fSin + 1.0f);
						subd[i3] = g * r;
						fCos = 1.0f / r;
						fSin *= fCos;
					}

					g = diag[i3] - p;
					r = (diag[i3 - 1] - g) * fSin + 2.0f * b * fCos;
					p = fSin * r;
					diag[i3] = g + p;
					g = fCos * r - b;

					for (i4 = 0; i4 <n; i4++)
					{
						f = pntr[i4 * n + i3];
						pntr[i4 * n + i3] = fSin * pntr[i4 * n + i3 - 1] + fCos * f;
						pntr[i4 * n + i3 - 1] = fCos * pntr[i4 * n + i3 - 1] - fSin * f;
					}
				}
				diag[i0] -= p;
				subd[i0] = g;
				subd[i2] = 0.0f;
			}
			if (i1 == iMaxIter)
				return false;
		}
		return true;
    }

    template<class T> inline void MatrixMN<T>::SortEigenInfo(VectorN<T>& diag) {
        int32 i0, i1, i2;
		for(i0 = 0; i0 <n - 1; ++i0) {
			// locate maximum eigenvalue
			i1 = i0;
			T max = fabs(diag[i1]);
			for(i2 = i0 + 1; i2 <n; ++i2) {
				if (fabs(diag[i2])> fabs(max)) {
					i1 = i2;
					max = diag[i1];
				}
			}

			if (i1 != i0) {
				// swap eigenvalues
				diag[i1] = diag[i0];
				diag[i0] = max;

				// swap eigenvectors
                for(i2 = 0; i2 <n; ++i2) {
                    std::swap(pntr[i2 * n + i0], pntr[i2 * n + i1]);
                }
			}
		}
    }

    template<class T> inline T MatrixMN<T>::Pythag(T a, T b) const {
        T absa,absb;
        absa=fabs(a);
        absb=fabs(b);
        if (absa> absb) 
            return absa * sqrt(1.0f + sqr(absb / absa));
        else if (absb != 0.0) {
            return absb * sqrt(1.0f +  sqr(absa / absb));
        }
        else {
            return 0.0;
        }
    }
}