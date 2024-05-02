#include <math.h>

namespace BlockyVoyages {
    template<class T> inline Matrix33<T>::Matrix33() {}

    template<class T> inline Matrix33<T>::Matrix33(const Matrix33<T>& other) 
        : m00(other.m00), m10(other.m10), m20(other.m20), 
          m01(other.m01), m11(other.m11), m21(other.m21),
          m02(other.m02), m12(other.m12), m22(other.m22)
    {
        memcpy(pntr, other.pntr, sizeof(T) * 4);
    }

    // Set the matrix to the values in v. v must be in column major form
    template<class T> inline Matrix33<T>::Matrix33(const T v[9]) 
        : m00(v[0]), m10(v[1]), m20(v[2]),
          m01(v[3]), m11(v[4]), m21(v[5]),
          m02(v[6]), m12(v[7]), m22(v[8]),
    {}

    template<class T> inline Matrix33<T>::Matrix33(const Matrix22<T>& other) 
        : m00(other.m00), m10(other.m10), m20(0),
          m01(other.m01), m11(other.m11), m21(0),
          m02(0),            m12(0),            m22(0)
    {}

    template<class T> inline Matrix33<T>::Matrix33(T _00, T _01, T _02,
                                                                         T _10, T _11, T _12,
                                                                         T _20, T _21, T _22)
        : m00(_00), m10(_10), m20(_20),
          m01(_01), m11(_11), m21(_21), 
          m02(_02), m12(_12), m22(_22)
    {}

    template<class T> inline Matrix33<T>::Matrix33(const Matrix44<T>& other) 
        : m00(other.m00), m10(other.m10), m20(other.m20), 
          m01(other.m01), m11(other.m11), m21(other.m21),
          m02(other.m02), m12(other.m12), m22(other.m22)
    {}

    template<class T> inline Matrix33<T>::Matrix33(const Quaternion<T>& q) {
        Set(q);
    }

    template<class T> inline Matrix33<T>::~Matrix33() {}

    template<class T> inline Matrix33<T>& Matrix33<T>::operator = (const Matrix33<T>& other) {
        m00 = other.m00; m10 = other.m10; m20 = other.m20;
        m01 = other.m01; m11 = other.m11; m21 = other.m21;
        m02 = other.m02; m12 = other.m12; m22 = other.m22;

        return *this;
    }

    // Set the matrix to the values in v. v must be in column major form
    template<class T> inline void Matrix33<T>::Set(const T v[9]) {
        memcpy(pntr, v, sizeof(T) * 9);
    }

    template<class T> inline void Matrix33<T>::Set(T _00, T _01, T _02,
                                                                         T _10, T _11, T _12, 
                                                                         T _20, T _21, T _22) {
        m00 = _00; m10 = _10; m20 = _20;
        m01 = _01; m11 = _11; m21 = _21;
        m02 = _02; m12 = _12; m22 = _22;
    }

    template<class T> inline void Matrix33<T>::Set(const Matrix33<T>& other) {
        m00 = other.m00; m10 = other.m10; m20 = other.m20;
        m01 = other.m01; m11 = other.m11; m21 = other.m21;
        m02 = other.m02; m12 = other.m12; m22 = other.m22;
    }

    template<class T> inline void Matrix33<T>::Set(const Matrix44<T>& other) {
        m00 = other.m00; m10 = other.m10; m20 = other.m20;
        m01 = other.m01; m11 = other.m11; m21 = other.m21;
        m02 = other.m02; m12 = other.m12; m22 = other.m22;
    }

    template<class T> inline void Matrix33<T>::Set(const Quaternion<T>& q) {
        T s, xs, ys, zs, wx, wy, wz, xx, xy, xz, yy, yz, zz;

		s = 2.0f / (q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
		xs = s * q.x;	ys = s * q.y;	zs = s * q.z;
		wx = q.w * xs;	wy = q.w * ys;	wz = q.w * zs;
		xx = q.x * xs;	xy = q.x * ys;	xz = q.x * zs;
		yy = q.y * ys;	yz = q.y * zs;	zz = q.z * zs;
          
		m00 = 1.0f - (yy + zz); m10 = xy + wz;             m20 = xz - wy;
        m01 = xy - wz;             m11 = 1.0f - (xx + zz); m21 = yz + wx;
        m02 = xz + wy;             m12 = yz - wx;             m22 = 1.0f - (xx + yy);
    }

    template<class T> inline void Matrix33<T>::SetRow(const Vector3<T>& vals, int32 row) {
        pntr[row     ] = vals.x;
        pntr[row + 3] = vals.y;
        pntr[row + 6] = vals.z;
    }

    template<class T> inline void Matrix33<T>::SetRows(const Vector3<T>& row0, const Vector3<T>& row1, const Vector3<T>& row2) {
        m00 = row0.x;
        m10 = row1.x;
        m20 = row2.x;

        m01 = row0.y;
        m11 = row1.y;
        m21 = row2.y;
        
        m02 = row0.z;
        m12 = row1.z;
        m22 = row2.z;
    }

    template<class T> inline void Matrix33<T>::SetColumn(const Vector3<T>& vals, int32 col) {
        col = col * 3;
        pntr[col     ] = vals.x;
        pntr[col + 1] = vals.y;
        pntr[col + 2] = vals.z;
    }

    template<class T> inline void Matrix33<T>::SetColumns(const Vector3<T>& col0, const Vector3<T>& col1, const Vector3<T>& col2) {
        m00 = col0.x;
        m10 = col0.y;
        m20 = col0.z;

        m01 = col1.x;
        m11 = col1.y;
        m21 = col1.z;

        m02 = col2.x;
        m12 = col2.y;
        m22 = col2.z;
    }

    template<class T> inline Vector3<T> Matrix33<T>::GetRow(int32 row) const {
        return Vector3<T>(pntr[row], pntr[row + 3], pntr[row + 6]);
    }

    template<class T> inline void Matrix33<T>::GetRows(Vector3<T>& row0, Vector3<T>& row1, Vector3<T>& row2) const {
        row0.Set(m00, m01, m02);
        row1.Set(m10, m11, m12);
        row2.Set(m20, m21, m22);
    }

    template<class T> inline Vector3<T> Matrix33<T>::GetColumn(int32 col) const {
        return Vector3<T>(&pntr[col * 3]);
    }

    template<class T> inline void Matrix33<T>::GetColumns(Vector3<T>& col0, Vector3<T>& col1, Vector3<T>& col2) const {
        col0.Set(pntr        );
        col1.Set(&pntr[3]);
        col2.Set(&pntr[6]);
    }

    template<class T> inline void Matrix33<T>::Clean(void) {
        for(int32 ind = 0; ind <9; ++ind) {
            if (fabs(pntr[ind]) <EPSILON) {
                pntr[ind] = 0.0f;
            }
        }
    }

    template<class T> inline void Matrix33<T>::Zero(void) {
        for(int32 ind = 0; ind <9; ++ind) {
            pntr[ind] = 0;
        }
    }

    template<class T> inline void Matrix33<T>::Identity(void) {
        m00 = m11 = m22 = 1;
        m10 = m20 = m01 =
        m21 = m02 = m12 = 0;
    }

    template<class T> inline void Matrix33<T>::Rotation(T angle, T x, T y, T z) {
        throw "Can't use method Matrix22<T>::Rotation with the requested type.";
    }

    template<> inline void Matrix33<float32>::Rotation(float32 angle, float32 x, float32 y, float32 z) {
        float32 fCos = cos(angle);
		float32 fSin = sin(angle);
		float32 fSum = 1.0f - fCos;

		float32 vectLen = x * x + y * y + z * z;

		if (vectLen> TINY) {
			vectLen = 1 / sqrt(vectLen);
			x *= vectLen;
			y *= vectLen;
			z *= vectLen;
		}

		m00 = (x * x) * fSum + fCos;
		m01 = (x * y) * fSum - (z * fSin);
		m02 = (x * z) * fSum + (y * fSin);

		m10 = (y * x) * fSum + (z * fSin);
		m11 = (y * y) * fSum + fCos;
		m12 = (y * z) * fSum - (x * fSin);

		m20 = (z * x) * fSum - (y * fSin);
		m21 = (z * y) * fSum + (x * fSin);
		m22 = (z * z) * fSum + fCos;
    }

    template<> inline void Matrix33<float64>::Rotation(float64 angle, float64 x, float64 y, float64 z) {
        float64 fCos = cos(angle);
		float64 fSin = sin(angle);
		float64 fSum = 1.0f - fCos;

		float64 vectLen = x * x + y * y + z * z;

		if (vectLen> TINY) {
			vectLen = 1 / sqrt(vectLen);
			x *= vectLen;
			y *= vectLen;
			z *= vectLen;
		}

		m00 = (x * x) * fSum + fCos;
		m01 = (x * y) * fSum - (z * fSin);
		m02 = (x * z) * fSum + (y * fSin);

		m10 = (y * x) * fSum + (z * fSin);
		m11 = (y * y) * fSum + fCos;
		m12 = (y * z) * fSum - (x * fSin);

		m20 = (z * x) * fSum - (y * fSin);
		m21 = (z * y) * fSum + (x * fSin);
		m22 = (z * z) * fSum + fCos;
    }

    

    template<class T> inline void Matrix33<T>::Scale(T x, T y, T z) {
        m01 = m02 = m10 = m12 = m20 = m21 = 0.0f;
        m00 = x;
        m11 = y;
        m22 = z;
    }

    template<class T> inline void Matrix33<T>::Shear(const Vector3<T>& normal, const Vector3<T>& shear) {
        m00 = 1 + normal.x * shear.x;
		m01 =      normal.y * shear.x;
		m02 =      normal.z * shear.x;

		m10 =      normal.x * shear.y;
		m11 = 1 + normal.y * shear.y;
		m12 =      normal.z * shear.y;

		m20 =      normal.x * shear.z;
		m21 =      normal.y * shear.z;
		m22 = 1 + normal.z * shear.z;
    }

    template<class T> inline void Matrix33<T>::TransposeOf(const Matrix33<T>& mat) {
        m00 = mat.m00; m01 = mat.m10; m02 = mat.m20;
        m10 = mat.m01; m11 = mat.m11; m12 = mat.m21;
        m20 = mat.m02; m21 = mat.m12; m22 = mat.m22;
    }

    template<class T> inline void Matrix33<T>::Transpose(void) {
        std::swap(m01, m10);
        std::swap(m02, m20);
        std::swap(m12, m21);
    }

    template<class T> inline bool Matrix33<T>::InverseOf(const Matrix33<T>& mat) {
        // compute determinant
        T cofactor0 = m11 * m22 - m21 * m12;
		T cofactor3 = m20 * m12 - m10 * m22;
		T cofactor6 = m10 * m21 - m20 * m11;

		T det = m00 * cofactor0 + m01 * cofactor3 + m02 * cofactor6;
		if (IsZero(det)) {
			return false;
		}

		// create adjoint matrix and multiply by 1/det to get inverse
		m00 = cofactor0 / det;
		m10 = cofactor3 / det;
		m20 = cofactor6 / det;
	    
		m01 = (mat.m21 * mat.m02 - mat.m01 * mat.m22) / det;
		m11 = (mat.m00 * mat.m22 - mat.m20 * mat.m02) / det;
		m21 = (mat.m20 * mat.m01 - mat.m00 * mat.m21) / det;

		m02 = (mat.m01 * mat.m12 - mat.m11 * mat.m02) / det;
		m12 = (mat.m10 * mat.m02 - mat.m00 * mat.m12) / det;
		m22 = (mat.m00 * mat.m11 - mat.m10 * mat.m01) / det;

		return true;
    }

    template<> inline bool Matrix33<float32>::InverseOf(const Matrix33<float32>& mat) {
        // compute determinant
        float32 cofactor0 = m11 * m22 - m21 * m12;
		float32 cofactor3 = m20 * m12 - m10 * m22;
		float32 cofactor6 = m10 * m21 - m20 * m11;

		float32 det = m00 * cofactor0 + m01 * cofactor3 + m02 * cofactor6;;
		if (IsZero(det)) {
			return false;
		}

		// create adjoint matrix and multiply by 1/det to get inverse
		det = 1.0f / det;
		m00 = det * cofactor0;
		m10 = det * cofactor3;
		m20 = det * cofactor6;
	    
		m01 = det * (mat.m21 * mat.m02 - mat.m01 * mat.m22);
		m11 = det * (mat.m00 * mat.m22 - mat.m20 * mat.m02);
		m21 = det * (mat.m20 * mat.m01 - mat.m00 * mat.m21);

		m02 = det * (mat.m01 * mat.m12 - mat.m11 * mat.m02);
		m12 = det * (mat.m10 * mat.m02 - mat.m00 * mat.m12);
		m22 = det * (mat.m00 * mat.m11 - mat.m10 * mat.m01);

		return true;
    }

    template<> inline bool Matrix33<float64>::InverseOf(const Matrix33<float64>& mat) {
        // compute determinant
        float64 cofactor0 = m11 * m22 - m21 * m12;
		float64 cofactor3 = m20 * m12 - m10 * m22;
		float64 cofactor6 = m10 * m21 - m20 * m11;

		float64 det = m00 * cofactor0 + m01 * cofactor3 + m02 * cofactor6;;
		if (IsZero(det)) {
			return false;
		}

		// create adjoint matrix and multiply by 1/det to get inverse
		det = 1.0f / det;
		m00 = det * cofactor0;
		m10 = det * cofactor3;
		m20 = det * cofactor6;
	    
		m01 = det * (mat.m21 * mat.m02 - mat.m01 * mat.m22);
		m11 = det * (mat.m00 * mat.m22 - mat.m20 * mat.m02);
		m21 = det * (mat.m20 * mat.m01 - mat.m00 * mat.m21);

		m02 = det * (mat.m01 * mat.m12 - mat.m11 * mat.m02);
		m12 = det * (mat.m10 * mat.m02 - mat.m00 * mat.m12);
		m22 = det * (mat.m00 * mat.m11 - mat.m10 * mat.m01);

		return true;
    }

    template<class T> inline bool Matrix33<T>::Invert(void) {
        Matrix33<T> mat;
        if (mat.InverseOf(*this)) {
            Set(mat);
            return true;
        }
        return false;
    }

    template<class T> inline T Matrix33<T>::Determinant(void) const {
        T cofactor0 = m11 * m22 - m21 * m12;
		T cofactor3 = m20 * m12 - m10 * m22;
		T cofactor6 = m10 * m21 - m20 * m11;
        return m00 * cofactor0 + m01 * cofactor3 + m02 * cofactor6;
    }

    template<class T> inline T Matrix33<T>::Trace(void) const {
        return m00 + m11 + m22;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::Tensor(const Vector3<T>& v1, const Vector3<T>& v2) {
        m00 = v1.x * v2.x;  m01 = v1.x * v2.y; m02 = v1.x * v2.z;
        m10 = v1.y * v2.x;  m11 = v1.y * v2.y; m12 = v1.y * v2.z;
        m20 = v1.z * v2.x;  m21 = v1.z * v2.y; m22 = v1.z * v2.z;
    }

    template<class T> inline const Matrix33<T> Matrix33<T>::operator * (const Matrix33<T> &rhs) const {
        return Matrix33<T>(m00 * rhs.m00 + m01 * rhs.m10 + m02 * rhs.m20, 
                                     m00 * rhs.m01 + m01 * rhs.m11 + m02 * rhs.m21, 
                                     m00 * rhs.m02 + m01 * rhs.m12 + m02 * rhs.m22,

                                     m10 * rhs.m00 + m11 * rhs.m10 + m12 * rhs.m20, 
                                     m10 * rhs.m01 + m11 * rhs.m11 + m12 * rhs.m21, 
                                     m10 * rhs.m02 + m11 * rhs.m12 + m12 * rhs.m22,

                                     m20 * rhs.m00 + m21 * rhs.m10 + m22 * rhs.m20, 
                                     m20 * rhs.m01 + m21 * rhs.m11 + m22 * rhs.m21, 
                                     m20 * rhs.m02 + m21 * rhs.m12 + m22 * rhs.m22);
    }

    template<class T> inline const Matrix33<T> Matrix33<T>::operator * (        T    rhs) const {
        return Matrix33<T>(*this) *= rhs;
    }

    template<class T> inline const Matrix33<T> Matrix33<T>::operator / (        T    rhs) const {
        return Matrix33<T>(*this) /= rhs;
    }

    template<class T> inline const Matrix33<T> Matrix33<T>::operator + (const Matrix33<T> &rhs) const {
        return Matrix33<T>(*this) += rhs;
    }

    template<class T> inline const Matrix33<T> Matrix33<T>::operator - (const Matrix33<T> &rhs) const {
        return Matrix33<T>(*this) -= rhs;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::operator *=(const Matrix33<T> &rhs) {
        Matrix33<T> temp(m00 * rhs.m00 + m01 * rhs.m10 + m02 * rhs.m20, 
                                  m00 * rhs.m01 + m01 * rhs.m11 + m02 * rhs.m21, 
                                  m00 * rhs.m02 + m01 * rhs.m12 + m02 * rhs.m22,

                                  m10 * rhs.m00 + m11 * rhs.m10 + m12 * rhs.m20, 
                                  m10 * rhs.m01 + m11 * rhs.m11 + m12 * rhs.m21, 
                                  m10 * rhs.m02 + m11 * rhs.m12 + m12 * rhs.m22,

                                  m20 * rhs.m00 + m21 * rhs.m10 + m22 * rhs.m20, 
                                  m20 * rhs.m01 + m21 * rhs.m11 + m22 * rhs.m21, 
                                  m20 * rhs.m02 + m21 * rhs.m12 + m22 * rhs.m22);

        (*this) = temp;
        return *this;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::operator *=(        T    rhs) {
        m00 *= rhs;
        m10 *= rhs;
        m20 *= rhs;

        m01 *= rhs;
        m11 *= rhs;
        m21 *= rhs;

        m02 *= rhs;
        m12 *= rhs;
        m22 *= rhs;

        return *this;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::operator /=(        T    rhs) {
        m00 /= rhs;
        m10 /= rhs;
        m20 /= rhs;

        m01 /= rhs;
        m11 /= rhs;
        m21 /= rhs;

        m02 /= rhs;
        m12 /= rhs;
        m22 /= rhs;

        return *this;
    }

    template<> inline Matrix33<float32>& Matrix33<float32>::operator /=(float32 rhs) {
        rhs = 1.0f / rhs;
        m00 *= rhs;
        m10 *= rhs;
        m20 *= rhs;

        m01 *= rhs;
        m11 *= rhs;
        m21 *= rhs;

        m02 *= rhs;
        m12 *= rhs;
        m22 *= rhs;

        return *this;
    }

    template<> inline Matrix33<float64>& Matrix33<float64>::operator /=(float64 rhs) {
        rhs = 1.0 / rhs;
        m00 *= rhs;
        m10 *= rhs;
        m20 *= rhs;

        m01 *= rhs;
        m11 *= rhs;
        m21 *= rhs;

        m02 *= rhs;
        m12 *= rhs;
        m22 *= rhs;

        return *this;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::operator +=(const Matrix33<T> &rhs) {
        m00 += rhs.m00;
        m01 += rhs.m01;
        m02 += rhs.m02;

        m10 += rhs.m10;
        m11 += rhs.m11;
        m12 += rhs.m12;

        m20 += rhs.m20;
        m21 += rhs.m21;
        m22 += rhs.m22;

        return *this;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::operator -=(const Matrix33<T> &rhs) {
        m00 -= rhs.m00;
        m10 -= rhs.m10;
        m20 -= rhs.m20;

        m01 -= rhs.m01;
        m11 -= rhs.m11;
        m21 -= rhs.m21;

        m02 -= rhs.m02;
        m12 -= rhs.m12;
        m22 -= rhs.m22;

        return *this;
    }

    // Methods used to replace the = sign. Use as ans = m1 op m2
    template<class T> inline Matrix33<T>& Matrix33<T>::Add(const Matrix33<T>& m1, const Matrix33<T>& m2) {
        return (*this = m1) += m2;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::Sub(const Matrix33<T>& m1, const Matrix33<T>& m2) {
        return (*this = m1) -= m2;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::Mul(const Matrix33<T>& m1, const Matrix33<T>& m2) {
        m00 = m1.m00 * m2.m00 + m1.m01 * m2.m10 + m1.m02 * m2.m20; 
        m01 = m1.m00 * m2.m01 + m1.m01 * m2.m11 + m1.m02 * m2.m21; 
        m02 = m1.m00 * m2.m02 + m1.m01 * m2.m12 + m1.m02 * m2.m22;

        m10 = m1.m10 * m2.m00 + m1.m11 * m2.m10 + m1.m12 * m2.m20; 
        m11 = m1.m10 * m2.m01 + m1.m11 * m2.m11 + m1.m12 * m2.m21; 
        m12 = m1.m10 * m2.m02 + m1.m11 * m2.m12 + m1.m12 * m2.m22;

        m20 = m1.m20 * m2.m00 + m1.m21 * m2.m10 + m1.m22 * m2.m20; 
        m21 = m1.m20 * m2.m01 + m1.m21 * m2.m11 + m1.m22 * m2.m21; 
        m22 = m1.m20 * m2.m02 + m1.m21 * m2.m12 + m1.m22 * m2.m22;

        return *this;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::Mul(const Matrix33<T>& m1, T    f) {
        return (*this = m1) *= f;
    }

    template<class T> inline Matrix33<T>& Matrix33<T>::Div(const Matrix33<T>& m1, T    f) {
        return (*this = m1) /= f;
    }

    template<class T> inline const Vector3<T> Matrix33<T>::operator * (const Vector3<T>  &rhs) const {
        return Vector3<T>(m00 * rhs.x + m01 * rhs.y + m02 * rhs.z,
                                    m10 * rhs.x + m11 * rhs.y + m12 * rhs.z,
                                    m20 * rhs.x + m21 * rhs.y + m22 * rhs.z);
    }

    template<class T> inline bool Matrix33<T>::operator==(const Matrix33<T>& other) const {
        if (this == &other) {
            return true;
        }

        for(int32 ind = 0; ind <9; ++ind) {
            if (pntr[ind] != other.pntr[ind]) {
                return false;
            }
        }

        return true;
    }

    template<class T> inline bool Matrix33<T>::operator!=(const Matrix33<T>& other) const {
        return !((*this) == other);
    }

    template<class T> inline bool Matrix33<T>::GetSymEigenInfo(Matrix33<T>& ev, Vector3<T>& ew) const {
        ev = *this;
        return ev.SymEigenDecompose(ew);
    }

    template<class T> inline bool Matrix33<T>::SymEigenDecompose(Vector3<T>& ew) {
        throw "Can't use method Matrix33<T>::SymEigenDecompose with the requested type.";
        return false;
    }

    template<> inline bool Matrix33<float32>::SymEigenDecompose(Vector3<float32>& ew) {
        Vector3<float32> subd;
        Tridiagonalize(ew, subd);
        if (QLFactorize(ew, subd)) {
            // sort them
            SortEigenInfo(ew);
            return true;
        }
        return false;
    }

    template<> inline bool Matrix33<float64>::SymEigenDecompose(Vector3<float64>& ew) {
        Vector3<float64> subd;
        Tridiagonalize(ew, subd);
        if (QLFactorize(ew, subd)) {
            // sort them
            SortEigenInfo(ew);
            return true;
        }
        return false;
    }

    template<class T> inline void Matrix33<T>::Tridiagonalize(Vector3<T>& diag, Vector3<T>& subd) {
        int32 i0, i1, i2;

        T scale;
        scale = fabs(pntr[2]) + fabs(pntr[5]);

        if (scale == 0) {
            subd[2] = pntr[5];

            diag[2] = 0;
        }
        else
        {
            T invScale = 1.0f / scale;

            pntr[2] *= invScale;
            pntr[5] *= invScale;

            T h = pntr[2] * pntr[2] + pntr[5] * pntr[5];

            T f = pntr[5];
            T g = sqrt(h);
            if (f> 0.0f)
                g = -g;

            subd[2] = scale * g;
            h -= f * g;
            pntr[5] = f - g;

            T invH = 1.0f / h;

            pntr[6] = pntr[2] * invH;
            pntr[7] = pntr[5] * invH;

            g = pntr[0] * pntr[2] + pntr[1] * pntr[5];
            subd[0] = g * invH;

            g = pntr[1] * pntr[2] + pntr[4] * pntr[5];
            subd[1] = g * invH;

            f = subd[0] * pntr[2] + subd[1] * pntr[5];

            T halfdivH = 0.5f * f * invH;

            f = pntr[2];
            g = subd[0] - halfdivH*f;
            subd[0] = g;
            pntr[0] -= f * subd[0] + g * pntr[2];

            f = pntr[5];
            g = subd[1] - halfdivH*f;
            subd[1] = g;
            pntr[1] -= f * subd[0] + g * pntr[2];
            pntr[4] -= f * subd[1] + g * pntr[5];

            diag[2] = h;
        }

        subd[1] = pntr[1];
        diag[1] = 0.0f;

        diag[0] = subd[0] = 0;
        for (i0 = 0; i0 <3; ++i0)
        {
            if (diag[i0] != 0)
            {
                for (i1 = 0; i1 <i0; ++i1)
                {
                    T sum = 0.0f;
                    for (i2 = 0; i2 <i0; ++i2)
                        sum += pntr[i0 + i2 * 3] * pntr[i2 + i1 * 3];

                    for (i2 = 0; i2 <i0; ++i2)
                        pntr[i2 + i1 * 3] -= sum * pntr[i2 + i0 * 3];
                }
            }

            diag[i0] = pntr[i0 * 4];
            pntr[i0 * 4] = 1;
            for (i1 = 0; i1 <i0; ++i1)
                pntr[i1 + i0 * 3] = pntr[i0 + i1 * 3] = 0.0f;
        }

        // shift the sub diagonal up one
        subd[0] = subd[1];
        subd[1] = subd[2];
        subd[2] = 0;
        //*/
    }

    template<class T> inline bool Matrix33<T>::QLFactorize(Vector3<T>& diag, Vector3<T>& subd) {
        const int32 iMaxIter = 32;

		int32 i0, i1, i2, i3, i4;
		for(i0 = 0; i0 <3; i0++) {
			for (i1 = 0; i1 <iMaxIter; i1++)
			{
				i2 = i0;
				while(i2 <2 && subd[i2] != 0)
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

					for (i4 = 0; i4 <3; i4++)
					{
						f = pntr[i4 + i3 * 3];
						pntr[i4 +    i3         * 3] = fSin * pntr[i4 + (i3 - 1) * 3] + fCos * f;
						pntr[i4 + (i3 - 1) * 3] = fCos * pntr[i4 + (i3 - 1) * 3] - fSin * f;
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

    template<class T> inline void Matrix33<T>::SortEigenInfo(Vector3<T>& diag) {
        int32 i0, i1, i2;
		for(i0 = 0; i0 <2; ++i0) {
			// locate maximum eigenvalue
			i1 = i0;
			T max = fabs(diag[i1]);
			for(i2 = i0 + 1; i2 <3; ++i2) {
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
                for(i2 = 0; i2 <3; ++i2) {
                    std::swap(pntr[i2 + i0 * 3], pntr[i2 + i1 * 3]);
                }
			}
		}
    }
}