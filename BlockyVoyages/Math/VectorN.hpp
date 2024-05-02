#include <assert.h>
#include <math.h>

namespace BlockyVoyages {
    template<class T> inline VectorN<T>::VectorN(void) 
        : n(0), pntr(nullptr)
    {}

    template<class T> inline VectorN<T>::VectorN(int32 _n, T value)
        : n(_n), pntr(new T[_n])
    {
        for (int32 i = 0; i < n; ++i) {
            pntr[i] = value;
        }
    }

    template<class T> inline VectorN<T>::VectorN(const VectorN<T>& other)
        : n(other.n), pntr(new T[other.n])
    {
        memcpy(pntr, other.pntr, sizeof(T) * n);
    }

    template<class T> inline VectorN<T>::VectorN(const Vector2<T>& other) 
        : n(2), pntr(new T[2])
    {
        pntr[0] = other.x;
        pntr[1] = other.y;
    }

    template<class T> inline VectorN<T>::VectorN(const Vector3<T>& other)
        : n(3), pntr(new T[3])
    {
        pntr[0] = other.x;
        pntr[1] = other.y;
        pntr[2] = other.z;
    }

    template<class T> inline VectorN<T>::VectorN(const Vector4<T>& other)
        : n(4), pntr(new T[4])
    {
        pntr[0] = other.x;
        pntr[1] = other.y;
        pntr[2] = other.z;
        pntr[3] = other.w;
    }

    // set from floats
    template<class T> inline VectorN<T>::VectorN(const T* _v, int32 _n)
        : n(_n), pntr(new T[_n])
    {
        memcpy(pntr, _v, sizeof(T) * n);
    }

    template<class T> inline VectorN<T>::VectorN(int32 _n) 
        : n(_n), pntr(new T[_n])
    {}

    template<class T> inline VectorN<T>::~VectorN(void) {
        delete[] pntr;
    }

     // sets the size of the array. Also initializes the vector to value
    template<class T> inline void VectorN<T>::SetSize(int32 nSlots, T value) {
        if (nSlots != n) {
            delete[] pntr;
            n = nSlots;
            pntr = new T[nSlots];
        }
        for (int32 i = 0; i < n; ++i) {
            pntr[i] = value;
        }
    }

    template<class T> inline void VectorN<T>::SetSize(int32 nSlots) {
        if (nSlots != n) {
            delete[] pntr;

            n = nSlots;
            pntr = new T[nSlots];
        }
    }

    template<class T> inline void VectorN<T>::Clear(void) {
        delete[] pntr;
        pntr = nullptr;

        n = 0;
    }

    // setting methods
    template<class T> inline VectorN<T>& VectorN<T>::operator=(const VectorN<T>& other) {
        Set(other);
        return *this;
    }

    template<class T> inline void VectorN<T>::Set(const Vector2<T>& other) {
        SetSize(2);
        pntr[0] = other.x;
        pntr[1] = other.y;
    }

    template<class T> inline void VectorN<T>::Set(const Vector3<T>& other) {
        SetSize(3);
        pntr[0] = other.x;
        pntr[1] = other.y;
        pntr[2] = other.z;
    }

    template<class T> inline void VectorN<T>::Set(const Vector4<T>& other) {
        SetSize(4);
        pntr[0] = other.x;
        pntr[1] = other.y;
        pntr[2] = other.z;
        pntr[3] = other.w;
    }

    template<class T> inline void VectorN<T>::Set(const VectorN<T>& other) {
        SetSize(other.n);
        memcpy(pntr, other.pntr, sizeof(T) * n);
    }

    // set from floats
    template<class T> inline void VectorN<T>::Set(const T* _v, int32 _n) {
        SetSize(_n);
        memcpy(pntr, _v, sizeof(T) * n);
    }

    // clean up and fast implementation
    template<class T> inline void VectorN<T>::Zero(void) {
        for(int32 ind = 0; ind <n; ++ind) {
            pntr[ind] = 0.0f;
        }
    }

    template<class T> inline void VectorN<T>::Clean(void) {
        for(int32 ind = 0; ind <n; ++ind) {
            if (abs(pntr[ind]) <EPSILON) {
                pntr[ind] = 0.0f;
            }
        }
    }

    template<class T> inline VectorN<T>& VectorN<T>::Add(const VectorN<T>&  lhs, const VectorN<T>& rhs) {
        if (lhs.n == rhs.n) {
            SetSize(lhs.n);
            for(int32 ind = 0; ind <n; ++ind) {
                pntr[ind] = lhs.pntr[ind] + rhs.pntr[ind];
            }
        }
        else {
            Clear();
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::Sub(const VectorN<T>&  lhs, const VectorN<T>& rhs) {
        if (lhs.n == rhs.n) {
            SetSize(lhs.n);
            for(int32 ind = 0; ind <n; ++ind) {
                pntr[ind] = lhs.pntr[ind] - rhs.pntr[ind];
            }
        }
        else {
            Clear();
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::Mul(const VectorN<T>&  lhs, T rhs) {
        SetSize(lhs.n);
        for(int32 ind = 0; ind <n; ++ind) {
            pntr[ind] = lhs.pntr[ind] * rhs;
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::Mul(const VectorN<T>&  lhs, const MatrixMN<T>& rhs) {
        if (rhs.GetM() == lhs.n) {
            // row 0
            for(int32 col = 0; col <rhs.GetN(); ++col) {
                pntr[col] = rhs(0, col) * lhs.pntr[0];
            }

            // the rest of the matrix
            for(int32 row = 1; row <rhs.GetM(); ++row) {
                for(int32 col = 0; col <rhs.GetN(); ++col) {
                    pntr[col] += lhs.pntr[row] * rhs(row, col);
                }
            }
        }
        else {
            Clear();
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::Mul(const MatrixMN<T>& lhs, const VectorN<T>& rhs) {
        if (lhs.GetN() == rhs.n) {
            for(int32 row = 0; row <lhs.GetM(); ++row) {
                pntr[row] = rhs.pntr[0] * lhs(row, 0);
                for(int32 col = 1; col <lhs.GetN(); ++col) {
                    pntr[row] += rhs.pntr[0] * lhs(row, col);
                }
            }
        }
        else {
            Clear();
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::Div(const VectorN<T>&  lhs, T rhs) {
        SetSize(lhs.n);
        for(int32 ind = 0; ind <n; ++ind) {
            pntr[ind] = lhs.pntr[ind] / rhs;
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::Build(const VectorN<T>& init, const VectorN<T>& term) {
        return Sub(term, init);
    }

    template<class T> inline T VectorN<T>::Dot(const VectorN<T>& rhs) const {                 // computes the dot product between this and v2
        if (rhs.n == n && n> 0) {
            T sum = rhs.pntr[0] * pntr[0];
            for(int32 ind = 1; ind <n; ++ind) {
                sum += rhs.pntr[ind] * pntr[ind];
            }
            return sum;
        }
        return 0.0f;
    }

    template<class T> inline const VectorN<T> VectorN<T>::operator + (const VectorN<T> &rhs) const {        // Adds two vectors. use Add() or += instead if possible
        if (rhs.n == n) {
            VectorN<T> ans(n);
            for(int32 ind = 0; ind <n; ++ind) {
                ans.pntr[ind] = pntr[ind] + rhs.pntr[ind];
            }
            return ans;
        }
        return VectorN<T>();
    }

    template<class T> inline const VectorN<T> VectorN<T>::operator - (const VectorN<T> &rhs) const {        // subtracts two vectors. use Sub() or -= instead if possible
        if (rhs.n == n) {
            VectorN<T> ans(n);
            for(int32 ind = 0; ind <n; ++ind) {
                ans.pntr[ind] = pntr[ind] - rhs.pntr[ind];
            }
            return ans;
        }
        return VectorN<T>();
    }

    template<class T> inline T VectorN<T>::operator * (const VectorN<T> &rhs) const {        // computes the dot product between this and rhs
        if (rhs.n == n && n> 0) {
            T sum = pntr[0] * rhs.pntr[0];
            for(int32 ind = 1; ind <n; ++ind) {
                sum += pntr[ind] * rhs.pntr[ind];
            }
            return sum;
        }
        return 0.0;
    }

    template<class T> inline const VectorN<T> VectorN<T>::operator * (T rhs) const {             // multiplies this by a scalar. use Mul() or *= instead if possible
        VectorN<T> ans(n);
        for(int32 ind = 0; ind < n; ++ind) {
            ans.pntr[ind] = pntr[ind] * rhs;
        }
        return ans;
    }

    template<class T> inline const VectorN<T> VectorN<T>::operator * (const MatrixMN<T>& rhs) const {
        if (rhs.GetM() == n) {
            // row 0
            VectorN<T> ans(n);
            for(int32 col = 0; col <rhs.GetN(); ++col) {
                ans.pntr[col] = rhs(0, col) * pntr[0];
            }

            // the rest of the matrix
            for(int32 row = 1; row <rhs.GetM(); ++row) {
                for(int32 col = 0; col <rhs.GetN(); ++col) {
                    ans.pntr[col] += pntr[row] * rhs(row, col);
                }
            }

            return ans;
        }
        else {
            return VectorN<T>();
        }
    }

    template<class T> inline const VectorN<T> VectorN<T>::operator / (T rhs) const {             // divides this by a scalar. use Div() or /= instead if possible
        VectorN<T> ans(n);
        for(int32 ind = 0; ind < n; ++ind) {
            ans.pntr[ind] = pntr[ind] / rhs;
        }
        return ans;
    }

    template<class T> inline VectorN<T>& VectorN<T>::operator += (const VectorN<T> &rhs) {     // does an in-place addition with rhs
        assert(rhs.n == n);
        for(int32 ind = 0; ind < n; ++ind) {
            pntr[ind] += rhs.pntr[ind];
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::operator -= (const VectorN<T> &rhs) {     // does an in-place subtraction with rhs
        assert(rhs.n == n);
        for(int32 ind = 0; ind < n; ++ind) {
            pntr[ind] -= rhs.pntr[ind];
        }
        
        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::operator *= (T rhs) {          // does an in-place scalar multiplication
        for(int32 ind = 0; ind < n; ++ind) {
            pntr[ind] *= rhs;
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::operator *= (const MatrixMN<T>& mat) {          // does an in-place matrix multiplication
        if (mat.GetM() == n) {
            // row 0
            VectorN<T> ans(n);
            for(int32 col = 0; col <mat.GetN(); ++col) {
                ans.pntr[col] = pntr[0] * mat(0, col);
            }

            // the rest of the matrix
            for(int32 row = 1; row <mat.GetM(); ++row) {
                for(int32 col = 0; col <mat.GetN(); ++col) {
                    ans.pntr[col] += pntr[row] * mat(row, col);
                }
            }

            for(int32 i = 0; i <n; ++i) {
                pntr[i] = ans.pntr[i];
            }
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::operator /= (T rhs) {          // does an in-place scalar division
        for(int32 ind = 0; ind <n; ++ind) {
            pntr[ind] /= rhs;
        }

        return *this;
    }

    template<class T> inline T VectorN<T>::MagnitudeSquared(void) const { // returns the squared 2-norm of this vector
        return Dot(*this);
    }

    template<class T> inline T VectorN<T>::Magnitude(void) const {
        return sqrt(MagnitudeSquared());
    }

    template<class T> inline T VectorN<T>::MaxElement(void) const {
        if (n> 0) {
            T max = pntr[0];
            for(int32 ind = 1; ind <n; ++ind) {
                if (pntr[ind]> max) {
                    max = pntr[ind];
                }
            }
            return max;
        }
        return 0.0;
    }

    template<class T> inline int32 VectorN<T>::MaxIndex(void) const { // returns the index of the maximum element
        if (n> 0) {
            T max = pntr[0];
            int32 maxInd = 0;
            for(int32 ind = 1; ind <n; ++ind) {
                if (pntr[ind]> max) {
                    max = pntr[ind];
                    maxInd = ind;
                }
            }
            return maxInd;
        }
        return 0;
    }

    template<class T> inline VectorN<T>& VectorN<T>::NegationOf(const VectorN<T>& other) {
        SetSize(other.n);
        for(int32 ind = 0; ind <n; ++ind) {
            pntr[ind] = -other.pntr[ind];
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::Negate(void) {
        for(int32 ind = 0; ind <n; ++ind) {
            pntr[ind] *= -1.0f;
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::UnitOf(const VectorN<T>& other) {
        SetSize(other.n);
        T mag = other.MagnitudeSquared();
        if (mag> TINY) {
            mag = 1.0f / sqrt(mag);
            for(int32 ind = 0; ind <n; ++ind) {
                pntr[ind] = other.pntr[ind] * mag;
            }
        }
        else {
            Zero();
        }

        return *this;
    }

    template<class T> inline VectorN<T>& VectorN<T>::Normalize(void) {
        T mag = MagnitudeSquared();
        if (mag> TINY) {
            mag = 1.0f / sqrt(mag);
            for(int32 ind = 0; ind <n; ++ind) {
                pntr[ind] *= mag;
            }
        }
        else {
            Zero();
        }

        return *this;
    }

    template<class T> inline const MatrixMN<T> VectorN<T>::Tensor(const VectorN<T> &rhs) const {  // computes the outer product of this and rhs
        MatrixMN<T> ans(n, rhs.n);
        for(int32 row = 0; row <n; ++row) {
            for(int32 col = 0; col <rhs.n; ++col) {
                ans(row, col) = pntr[row] * rhs.pntr[col];
            }
        }

        return ans;
    }
    
    template<class T> inline T VectorN<T>::CosAngleBetween(const VectorN<T>& vect) const { // computes the cosine of the angle between this vector and vect
        T mag1 = MagnitudeSquared();
        T mag2 = vect.MagnitudeSquared();

        if (mag1> TINY && mag2> TINY) {
            return Dot(vect) / sqrt(mag1 * mag2);
        }
        return 0.0f;
    }

    template<class T> inline T VectorN<T>::AngleBetween(const VectorN<T>& vect) const {     // computes the angle between this vector and vect
        return acos(CosAngleBetween(vect));
    }

    template<class T> inline T VectorN<T>::DistanceSquared(const VectorN<T>& v2) const { // finds the square of the distance between this and v2
        if (n == v2.n && n> 0) {
            T sum = pntr[0] - v2.pntr[0];
            sum *= sum;
            for(int32 ind = 1; ind <n; ++ind) {
                T diff = pntr[ind] - v2.pntr[ind];
                sum = diff * diff;
            }

            return sum;
        }
        return 0.0f;
    }

    template<class T> inline T VectorN<T>::Distance(const VectorN<T>& v2) const {          // finds the distance between this and v2
        return sqrt(DistanceSquared(v2));
    }

    template<class T> inline bool VectorN<T>::operator == (const VectorN<T>& rhs) const {          // compares this to rhs and returns if they are equal
        if (this == &rhs) {
            return true;
        }

        if (n != rhs.n) {
            return false;
        }

        for(int32 ind = 0; ind <n; ++ind) {
            if (pntr[ind] != rhs.pntr[ind]) {
                return false;
            }
        }

        return true;
    }

    template<class T> inline bool VectorN<T>::operator != (const VectorN<T>& rhs) const {          // compares this to rhs and returns if they are not equal
        return !((*this) == rhs);
    }
}