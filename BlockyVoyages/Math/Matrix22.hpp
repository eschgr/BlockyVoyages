#include <math.h>

namespace BlockyVoyages {
    // constructors
    template<class T> inline Matrix22<T>::Matrix22(void) {}
    template<class T> inline Matrix22<T>::Matrix22(const Matrix22<T>& other)
        : m00(other.m00), m10(other.m10),
          m01(other.m01), m11(other.m11)
    {}

    template<class T> inline Matrix22<T>::Matrix22(const Matrix33<T>& other) 
        : m00(other.m00), m10(other.m10),
          m01(other.m01), m11(other.m11)
    {}

    template<class T> inline Matrix22<T>::Matrix22(const Matrix44<T>& other) 
        : m00(other.m00), m10(other.m10),
          m01(other.m01), m11(other.m11)
    {}

    template<class T> inline Matrix22<T>::Matrix22(const T _v[4]) 
        : m00(_v[0]), m10(_v[1]), 
          m01(_v[2]), m11(_v[3])
    { // Set the matrix to the values in pntr. pntr must be in column major form
    }

    template<class T> inline Matrix22<T>::Matrix22(T _00, T _01,
                                                                         T _10, T _11)
        : m00(_00), m10(_10),
          m01(_01), m11(_11)
    {}

    template<class T> inline Matrix22<T>::~Matrix22(void) {}

    template<class T> inline Matrix22<T>& Matrix22<T>::operator = (const Matrix22<T>& other) {
        if (this != &other) {
            m00 = other.m00;
            m10 = other.m10;
            m01 = other.m01;
            m11 = other.m11;
        }

        return *this;
    }

    // set the matrix to a type other than Matrix22<T>
    template<class T> inline void Matrix22<T>::Set(const Matrix22<T>& other) {
        m00 = other.m00; m10 = other.m10;
        m01 = other.m01; m11 = other.m11;
    }

    template<class T> inline void Matrix22<T>::Set(const Matrix33<T>& other) {
        m00 = other.m00; m10 = other.m10;
        m01 = other.m01; m11 = other.m11;
    }

    template<class T> inline void Matrix22<T>::Set(const Matrix44<T>& other) {
        m00 = other.m00; m10 = other.m10;
        m01 = other.m01; m11 = other.m11;
    }

    template<class T> inline void Matrix22<T>::Set(const T _v[4]) {
        m00 = _v[0]; m10 = _v[1];
        m01 = _v[2]; m11 = _v[3];
    }

    template<class T> inline void Matrix22<T>::Set(T _00, T _01,
                                                                         T _10, T _11) {
        m00 = _00; m10 = _10;
        m01 = _01; m11 = _11;
    }

    template<class T> inline void Matrix22<T>::SetRow(const Vector2<T>& vals, int32 row) {
        pntr[row     ] = vals.pntr[0];
        pntr[row + 2] = vals.pntr[1];
    }

    template<class T> inline void Matrix22<T>::SetRows(const Vector2<T>& row0, const Vector2<T>& row1) {
        m00 = row0.pntr[0]; m01 = row0.pntr[1];
        m10 = row1.pntr[0]; m11 = row1.pntr[1];
    }

    template<class T> inline void Matrix22<T>::SetColumn(const Vector2<T>& vals, int32 col) {
        pntr[(col << 1)    ] = vals.pntr[0];
        pntr[(col << 1) + 1] = vals.pntr[1];
    }

    template<class T> inline void Matrix22<T>::SetColumns(const Vector2<T>& col0, const Vector2<T>& col1) {
        m00 = col0.pntr[0]; m10 = col0.pntr[1];
        m10 = col1.pntr[0]; m11 = col1.pntr[1];
    }

    template<class T> inline Vector2<T> Matrix22<T>::GetRow(int32 row) const {
        return Vector2<T>(pntr[row], pntr[row + 2]);
    }

    template<class T> inline void Matrix22<T>::GetRows(Vector2<T>& row0, Vector2<T>& row1) const {
        row0.pntr[0] = m00; row0.pntr[1] = m01;
        row1.pntr[0] = m10; row1.pntr[1] = m11;
    }

    template<class T> inline Vector2<T> Matrix22<T>::GetColumn(int32 col) const {
        return Vector2<T>(&pntr[2 * col]);
    }

    template<class T> inline void Matrix22<T>::GetColumns(Vector2<T>& col0, Vector2<T>& col1) const {
        col0.pntr[0] = m00; col1.pntr[0] = m01;
        col0.pntr[1] = m10; col1.pntr[1] = m11;
    }

    template<class T> inline void Matrix22<T>::Clean(void) {
        if (abs(m00) <EPSILON) {
            m00 = 0;
        }

        if (abs(m10) <EPSILON) {
            m10 = 0;
        }

        if (abs(m01) <EPSILON) {
            m01 = 0;
        }

        if (abs(m11) <EPSILON) {
            m11 = 0;
        }
    }

    template<class T> inline void Matrix22<T>::Zero(void) {
        m00 = m10 = m01 = m11 = 0;
    }

    template<class T> inline void Matrix22<T>::Identity(void)  {
        m00 = m11 = 1;
        m01 = m10 = 0;
    }

    // Creates a rotation matrix which rotates about the z axis
    template<class T> inline void Matrix22<T>::Rotation(T ang) {
        m00 = m11 = cos(ang);
        m10 = sin(ang);
        m01 = -m10;
    }

    template<class T> inline void Matrix22<T>::Reflection(T ang) {
        m00 = cos(ang);
        m01 = m10 = sin(ang);
        m11 = -m00;
    }

    template<class T> inline void Matrix22<T>::TransposeOf(const Matrix22<T>& mat) {
        m00 = mat.m00;
        m10 = mat.m01;
        m01 = mat.m10;
        m11 = mat.m11;
    }

    template<class T> inline void Matrix22<T>::Transpose(void) {
        T temp = m10;
        m10 = m01;
        m01 = temp;
    }

    template<class T> inline bool Matrix22<T>::InverseOf(const Matrix22<T>& mat) {
        T det = mat.Determinant();

        if (det == 0) {
            return false;
        }

        // create adjoint matrix and multiply by 1/det to get inverse
        m00 =  mat.m11 / det;
        m10 = -mat.m10 / det;
        m01 = -mat.m01 / det;
        m11 =  mat.m00 / det;

        return true;
    }

    template<> inline bool Matrix22<float32>::InverseOf(const Matrix22<float32>& mat) {
        float32 det = mat.Determinant();

        if (abs(det) <TINY) {
            return false;
        }

        // create adjoint matrix and multiply by 1/det to get inverse
        det = 1.0f / det;
        m00 =  mat.m11 * det;
        m10 = -mat.m10 * det;
        m01 = -mat.m01 * det;
        m11 =  mat.m00 * det;

        return true;
    }

    template<> inline bool Matrix22<float64>::InverseOf(const Matrix22<float64>& mat) {
        float64 det = mat.Determinant();

        if (abs(det) <TINY) {
            return false;
        }

        det = 1.0 / det;
        // create adjoint matrix and multiply by 1/det to get inverse
        m00 =  mat.m11 * det;
        m10 = -mat.m10 * det;
        m01 = -mat.m01 * det;
        m11 =  mat.m00 * det;

        return true;
    }

    template<class T> inline bool Matrix22<T>::Invert(void) {
        T det = Determinant();

        if (det == 0) {
            return false;
        }

        // create adjoint matrix and multiply by 1/det to get inverse
        m10 /= -det;
        m01 /= -det;

        T t00 = m11 / det;
        m11 = m00 / det;
        m00 = t00;

        return true;
    }

    template<> inline bool Matrix22<float32>::Invert(void) {
        float32 det = Determinant();

        if (abs(det) <TINY) {
            return false;
        }

        // create adjoint matrix and multiply by 1/det to get inverse
        det = 1.0f / det;
        m10 *= -det;
        m01 *= -det;

        float32 t00 = det * m11;
        m11 = det * m00;
        m00 = t00;

        return true;
    }

    template<> inline bool Matrix22<float64>::Invert(void) {
        float64 det = Determinant();

        if (abs(det) <TINY) {
            return false;
        }

        // create adjoint matrix and multiply by 1/det to get inverse
        det = 1.0 / det;
        m10 *= -det;
        m01 *= -det;

        float64 t00 = det * m11;
        m11 = det * m00;
        m00 = t00;

        return true;
    }

    template<class T> inline T Matrix22<T>::Determinant(void) const {
        return m00 * m11 - m01 * m10;
    }

    template<class T> inline T Matrix22<T>::Trace(void) const {
        return m00 + m11;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::Tensor(const Vector2<T>& v1, const Vector2<T>& v2) {
        m00 = v1.x * v2.x;
        m01 = v1.x * v2.y;
        m10 = v1.y * v2.x;
        m11 = v1.y * v2.y;

        return *this;
    }

    template<class T> inline const Matrix22<T> Matrix22<T>::operator *(const Matrix22<T> &rhs) const {
        return Matrix22<T>(*this) *= rhs;
    }

    template<class T> inline const Matrix22<T> Matrix22<T>::operator *(        T    rhs) const {
        return Matrix22<T>(*this) *= rhs;
    }

    template<class T> inline const Matrix22<T> Matrix22<T>::operator /(        T    rhs) const {
        return Matrix22<T>(*this) /= rhs;
    }

    template<class T> inline const Matrix22<T> Matrix22<T>::operator +(const Matrix22<T> &rhs) const {
        return Matrix22<T>(*this) += rhs;
    }

    template<class T> inline const Matrix22<T> Matrix22<T>::operator -(const Matrix22<T> &rhs) const {
        return Matrix22<T>(*this) -= rhs;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::operator *=(const Matrix22<T> &rhs) {
        T tmp;
        tmp = m00 * rhs.m00 + m01 * rhs.m10;
        m01 = m00 * rhs.m01 + m01 * rhs.m11;
        m00 = tmp;

        tmp = m10 * rhs.m00 + m11 * rhs.m10;
        m11 = m10 * rhs.m01 + m11 * rhs.m11;
        m10 = tmp;

        return *this;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::operator *=(        T    rhs) {
        m00 *= rhs;
        m01 *= rhs;
        m10 *= rhs;
        m11 *= rhs;

        return *this;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::operator /=(        T    rhs) {
        m00 /= rhs;
        m01 /= rhs;
        m10 /= rhs;
        m11 /= rhs;

        return *this;
    }

    template<> inline Matrix22<float32>& Matrix22<float32>::operator /=(float32 rhs) {
        rhs = 1.0f / rhs;
        m00 *= rhs;
        m01 *= rhs;
        m10 *= rhs;
        m11 *= rhs;

        return *this;
    }

    template<> inline Matrix22<float64>& Matrix22<float64>::operator /=(float64 rhs) {
        rhs = 1.0 / rhs;
        m00 *= rhs;
        m01 *= rhs;
        m10 *= rhs;
        m11 *= rhs;

        return *this;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::operator +=(const Matrix22<T> &rhs) {
        m00 += rhs.m00;
        m01 += rhs.m01;
        m10 += rhs.m10;
        m11 += rhs.m11;

        return *this;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::operator -=(const Matrix22<T> &rhs) {
        m00 -= rhs.m00;
        m01 -= rhs.m01;
        m10 -= rhs.m10;
        m11 -= rhs.m11;

        return *this;
    }

    // Methods used to replace the = sign. Use as ans = m1 op m2 -> ans.op(m1, m2)
    template<class T> inline Matrix22<T>& Matrix22<T>::Add(const Matrix22<T>& m1, const Matrix22<T>& m2) {
        m00 = m1.m00 + m2.m00;
        m01 = m1.m01 + m2.m01;
        m10 = m1.m10 + m2.m10;
        m11 = m1.m11 + m2.m11;

        return *this;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::Sub(const Matrix22<T>& m1, const Matrix22<T>& m2) {
        m00 = m1.m00 + m2.m00;
        m01 = m1.m01 + m2.m01;
        m10 = m1.m10 + m2.m10;
        m11 = m1.m11 + m2.m11;

        return *this;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::Mul(const Matrix22<T>& m1, const Matrix22<T>& m2) {
        m00 = m1.m00 * m2.m00 + m1.m01 * m2.m10;
        m01 = m1.m00 * m2.m01 + m1.m01 * m2.m11;
        m10 = m1.m10 * m2.m00 + m1.m11 * m2.m10;
        m11 = m1.m10 * m2.m01 + m1.m11 * m2.m11;
        
        return *this;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::Mul(const Matrix22<T>& m1,         T f) {
        return (*this = m1) *= m2;
        return *this;
    }

    template<class T> inline Matrix22<T>& Matrix22<T>::Div(const Matrix22<T>& m1,         T f) {
        return (*this = m1) /= m2;
    }

    template<class T> inline Vector2<T>  Matrix22<T>::operator *(const Vector2<T> &rhs) const {
        return Vector2<T>(m00 * rhs.x + m01 * rhs.y,
                                    m10 * rhs.x + m11 * rhs.y);
    }

    template<class T> inline bool Matrix22<T>::operator==(const Matrix22<T>& other) const {
        return (this == &other) || (m00 == other.m00 && m10 == other.m10 && 
                                              m01 == other.m01 && m11 == other.m11);
    }

    template<class T> inline bool Matrix22<T>::operator!=(const Matrix22<T>& other) const {
        return !((*this) == other);
    }

    template<class T> inline bool Matrix22<T>::GetSymEigenInfo(Matrix22<T>& ev, Vector2<T>& ew) const {
        ev = *this;
        return ev.SymEigenDecompose(ew);
    }

    template<class T> inline bool Matrix22<T>::SymEigenDecompose(Vector2<T>& ew) {
        throw "Can't use method Matrix22<T>::SymEigenDecompose with the requested type.";
        return false;
    }

    template<> inline bool Matrix22<float32>::SymEigenDecompose(Vector2<float32>& ew) {
        Vector2<float32> subd;
        Tridiagonalize(ew, subd);
        if (QLFactorize(ew, subd)) {
            // sort them
            SortEigenInfo(ew);
            return true;
        }
        return false;
    }

    template<> inline bool Matrix22<float64>::SymEigenDecompose(Vector2<float64>& ew) {
        Vector2<float64> subd;
        Tridiagonalize(ew, subd);
        if (QLFactorize(ew, subd)) {
            // sort them
            SortEigenInfo(ew);
            return true;
        }
        return false;
    }

    template<class T> inline void Matrix22<T>::Tridiagonalize(Vector2<T>& diag, Vector2<T>& subd) {
        diag.x = m00;
		diag.y = m11;
		
        subd.x = m01;
		subd.y = 0.0f;

		Identity();
    }

    template<class T> inline bool Matrix22<T>::QLFactorize(Vector2<T>& diag, Vector2<T>& subd) {
        const int32 iMaxIter = 32;

		int32 iter, i4;
		for(iter = 0; iter <iMaxIter; ++iter) {
			if (subd[0] == 0) {
				break;
            }

			T g = (diag[1] - diag[0]) / (static_cast<T>(2.0) * subd[0]);
			T r = sqrt(g * g + static_cast<T>(1.0));

			if (g < 0.0f)
				g = diag[1] - diag[0] + subd[0] / (g - r);
			else
				g = diag[1] - diag[0] + subd[0] / (g + r);
			
			T fSin = static_cast<T>(1.0), fCos = static_cast<T>(1.0), p = static_cast<T>(0.0);
			T f = fSin * subd[0];
			T b = fCos * subd[0];
			if (fabs(f) >= fabs(g)) {
				fCos = g / f;
				r = sqrt(fCos * fCos + 1);
				fSin = static_cast<T>(1.0) / r;
				fCos *= fSin;
			}
			else {
				fSin = f / g;
				r = sqrt(fSin * fSin + 1);
				fCos = static_cast<T>(1.0) / r;
				fSin *= fCos;
			}

			g = diag[1] - p;
			r = (diag[0] - g) * fSin + static_cast<T>(2.0) * b * fCos;
			p = fSin * r;
			diag[1] = g + p;
			g = fCos * r - b;

			for (i4 = 0; i4 <2; i4++)
			{
				f = pntr[i4 + 2]; //mat(i4,1);
				//mat(i4, 1) = fSin * mat(i4, 0) + fCos * f;
                //mat(i4, 0) = fCos * mat(i4, 0) - fSin * f;

                pntr[i4 + 2] = fSin * pntr[i4] + fCos * f;
                pntr[i4     ] = fCos * pntr[i4] - fSin * f;
			}

			diag[0] -= p;
			subd[0] = g;
		}

        if (iter == iMaxIter)
			return false;
		
        return true;
    }

    template<class T> inline void Matrix22<T>::SortEigenInfo(Vector2<T>& diag) {
        if (fabs(diag[1])> fabs(diag[0])) {
			// swap eigenvalues
            T temp = diag[1];
			diag[1] = diag[0];
			diag[0] = temp;

			// swap eigenvectors
            for(int32 i2 = 0; i2 <2; ++i2) {
                std::swap(pntr[i2], pntr[i2 + 2]);
            }
		}
    }
}