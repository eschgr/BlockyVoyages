#include "stdafx.h"

#include "SimplexNoise.h"
#include "MathMain.h"

namespace BlockyVoyages {
namespace Math {

const Vector2f SimplexNoise::grad2[] = {
    Vector2f( 1.0f,  1.0f), Vector2f(-1.0f,  1.0f),
    Vector2f( 1.0f, -1.0f), Vector2f(-1.0f, -1.0f),
    Vector2f( 1.0f,  0.0f), Vector2f(-1.0f,  0.0f),
    Vector2f( 1.0f,  0.0f), Vector2f(-1.0f,  0.0f),
    Vector2f( 0.0f,  1.0f), Vector2f( 0.0f, -1.0f),
    Vector2f( 0.0f,  1.0f), Vector2f( 0.0f, -1.0f)
};

const Vector3f SimplexNoise::grad3[] = {
    Vector3f( 1.0f,  1.0f,  0.0f), Vector3f(-1.0f,  1.0f,  0.0f),
    Vector3f( 1.0f, -1.0f,  0.0f), Vector3f(-1.0f, -1.0f,  0.0f),
    Vector3f( 1.0f,  0.0f,  1.0f), Vector3f(-1.0f,  0.0f,  1.0f),
    Vector3f( 1.0f,  0.0f, -1.0f), Vector3f(-1.0f,  0.0f, -1.0f),
    Vector3f( 0.0f,  1.0f,  1.0f), Vector3f( 0.0f, -1.0f,  1.0f),
    Vector3f( 0.0f,  1.0f, -1.0f), Vector3f( 0.0f, -1.0f, -1.0f)
};

// this list is float32d
const uint8 SimplexNoise::perm[] = {
    151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7,
    225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190,
    6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117,
    35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136,
    171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146,
    158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41,
    55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80,
    73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116,
    188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226,
    250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207,
    206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170,
    213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167,
    43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232,
    178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238,
    210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239,
    107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115,
    121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67,
    29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,

    151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7,
    225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190,
    6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117,
    35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136,
    171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146,
    158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41,
    55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80,
    73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116,
    188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226,
    250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207,
    206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170,
    213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167,
    43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232,
    178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238,
    210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239,
    107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115,
    121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67,
    29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
};

// same list as above but everything mod 12
const uint8 SimplexNoise::permMod12[] = {
    7, 4, 5, 7, 6, 3, 11, 1, 9, 11, 0, 5, 2, 5, 7, 9, 8, 0, 7, 6, 9, 10, 8, 3,
    1, 0, 9, 10, 11, 10, 6, 4, 7, 0, 6, 3, 0, 2, 5, 2, 10, 0, 3, 11, 9, 11, 11,
    8, 9, 9, 9, 4, 9, 5, 8, 3, 6, 8, 5, 4, 3, 0, 8, 7, 2, 9, 11, 2, 7, 0, 3,
    10, 5, 2, 2, 3, 11, 3, 1, 2, 0, 7, 1, 2, 4, 9, 8, 5, 7, 10, 5, 4, 4, 6, 11,
    6, 5, 1, 3, 5, 1, 0, 8, 1, 5, 4, 0, 7, 4, 5, 6, 1, 8, 4, 3, 10, 8, 8, 3, 2,
    8, 4, 1, 6, 5, 6, 3, 4, 4, 1, 10, 10, 4, 3, 5, 10, 2, 3, 10, 6, 3, 10, 1,
    8, 3, 2, 11, 11, 11, 4, 10, 5, 2, 9, 4, 6, 7, 3, 2, 9, 11, 8, 8, 2, 8, 10,
    7, 10, 5, 9, 5, 11, 11, 7, 4, 9, 9, 10, 3, 1, 7, 2, 0, 2, 7, 5, 8, 4, 10,
    5, 4, 8, 2, 6, 1, 0, 11, 10, 2, 1, 10, 6, 0, 0, 11, 11, 6, 1, 9, 3, 1, 7,
    9, 2, 11, 11, 1, 0, 10, 7, 1, 7, 10, 1, 4, 0, 0, 8, 7, 1, 2, 9, 7, 4, 6, 2,
    6, 8, 1, 9, 6, 6, 7, 5, 0, 0, 3, 9, 8, 3, 6, 6, 11, 1, 0, 0,

    7, 4, 5, 7, 6, 3, 11, 1, 9, 11, 0, 5, 2, 5, 7, 9, 8, 0, 7, 6, 9, 10, 8, 3,
    1, 0, 9, 10, 11, 10, 6, 4, 7, 0, 6, 3, 0, 2, 5, 2, 10, 0, 3, 11, 9, 11, 11,
    8, 9, 9, 9, 4, 9, 5, 8, 3, 6, 8, 5, 4, 3, 0, 8, 7, 2, 9, 11, 2, 7, 0, 3,
    10, 5, 2, 2, 3, 11, 3, 1, 2, 0, 7, 1, 2, 4, 9, 8, 5, 7, 10, 5, 4, 4, 6, 11,
    6, 5, 1, 3, 5, 1, 0, 8, 1, 5, 4, 0, 7, 4, 5, 6, 1, 8, 4, 3, 10, 8, 8, 3, 2,
    8, 4, 1, 6, 5, 6, 3, 4, 4, 1, 10, 10, 4, 3, 5, 10, 2, 3, 10, 6, 3, 10, 1,
    8, 3, 2, 11, 11, 11, 4, 10, 5, 2, 9, 4, 6, 7, 3, 2, 9, 11, 8, 8, 2, 8, 10,
    7, 10, 5, 9, 5, 11, 11, 7, 4, 9, 9, 10, 3, 1, 7, 2, 0, 2, 7, 5, 8, 4, 10,
    5, 4, 8, 2, 6, 1, 0, 11, 10, 2, 1, 10, 6, 0, 0, 11, 11, 6, 1, 9, 3, 1, 7,
    9, 2, 11, 11, 1, 0, 10, 7, 1, 7, 10, 1, 4, 0, 0, 8, 7, 1, 2, 9, 7, 4, 6, 2,
    6, 8, 1, 9, 6, 6, 7, 5, 0, 0, 3, 9, 8, 3, 6, 6, 11, 1, 0, 0
};

// Skewing and unskewing factors for 2 and 3
const float32 SimplexNoise::F2 = 0.5f * (sqrt(3.0f) - 1.0f);
const float32 SimplexNoise::G2 = (3.0f - sqrt(3.0f)) / 6.0f;
const float32 SimplexNoise::F3 = 1.0f / 3.0f;
const float32 SimplexNoise::G3 = 1.0f / 6.0f;

float32 dot(Vector3f g, float32 x, float32 y) {
    return g.x * x + g.y * y;
}

float32 dot(Vector3f g, float32 x, float32 y, float32 z) {
    return g.x * x + g.y * y + g.z * z;
}

// 2D simplex noise
float32 SimplexNoise::noise2D(float32 xin, float32 yin) {
    float32 n0, n1, n2; // Noise contributions from the three corners
    // Skew the input space to determine which simplex cell we're in
    float32 s = (xin + yin) * F2; // Hairy factor for 2D
    int32 i = static_cast<int32>(fastfloor(xin + s));
    int32 j = static_cast<int32>(fastfloor(yin + s));
    float32 t = (i + j) * G2;
    float32 X0 = i - t; // Unskew the cell origin back to (x,y) space
    float32 Y0 = j - t;
    float32 x0 = xin - X0; // The x,y distances from the cell origin
    float32 y0 = yin - Y0;
    // For the 2D case, the simplex shape is an equilateral triangle.
    // Determine which simplex we are in.
    int i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
    if(x0 > y0) {
        // lower triangle, XY order: (0,0)->(1,0)->(1,1)
        i1 = 1;
        j1 = 0;
    }
    else {
        // upper triangle, YX order: (0,0)->(0,1)->(1,1)
        i1 = 0;
        j1 = 1;
    }
    // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
    // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
    // c = (3-sqrt(3))/6
    float32 x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
    float32 y1 = y0 - j1 + G2;
    float32 x2 = x0 - 1.0f + 2.0f * G2; // Offsets for last corner in (x,y) unskewed coords
    float32 y2 = y0 - 1.0f + 2.0f * G2;
    // Work out the hashed gradient indices of the three simplex corners
    int ii = i & 255;
    int jj = j & 255;
    int gi0 = permMod12[ii + perm[jj]];
    int gi1 = permMod12[ii + i1 + perm[jj + j1]];
    int gi2 = permMod12[ii + 1 + perm[jj + 1]];
    // Calculate the contribution from the three corners
    float32 t0 = 0.5f - x0 * x0 - y0 * y0;
    if (t0 < 0) {
        n0 = 0.0f;
    } else {
        t0 *= t0;
        n0 = t0 * t0 * grad2[gi0].Dot(Vector2f(x0, y0));  // (x,y) of grad used for 2D gradient
    }
    float32 t1 = 0.5f - x1 * x1 - y1 * y1;
    if (t1 < 0) {
        n1 = 0.0f;
    } else {
        t1 *= t1;
        n1 = t1 * t1 * grad2[gi1].Dot(Vector2f(x1, y1));
    }
    float32 t2 = 0.5f - x2 * x2 - y2 * y2;
    if (t2 < 0) {
        n2 = 0.0f;
    } else {
        t2 *= t2;
        n2 = t2 * t2 * grad2[gi2].Dot(Vector2f(x2, y2));
    }
    // Add contributions from each corner to get the final noise value.
    // The result is scaled to return values in the interval [-1,1].
    return 70.0f * (n0 + n1 + n2);
}


// 3D simplex noise
float32 SimplexNoise::noise3D(float32 xin, float32 yin, float32 zin) {
    float32 n0, n1, n2, n3; // Noise contributions from the four corners
    // Skew the input space to determine which simplex cell we're in
    float32 s = (xin + yin + zin) * F3; // Very nice and simple skew factor for 3D
    int32 i = static_cast<int32>(fastfloor(xin + s));
    int32 j = static_cast<int32>(fastfloor(yin + s));
    int32 k = static_cast<int32>(fastfloor(zin + s));
    float32 t = (i + j + k) * G3;
    float32 X0 = i - t; // Unskew the cell origin back to (x,y,z) space
    float32 Y0 = j - t;
    float32 Z0 = k - t;
    float32 x0 = xin - X0; // The x,y,z distances from the cell origin
    float32 y0 = yin - Y0;
    float32 z0 = zin - Z0;
    // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
    // Determine which simplex we are in.
    int i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
    int i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
    if(x0 >= y0) {
        if(y0 >= z0) {
            // X Y Z order
            i1 = 1; j1 = 0; k1 = 0;
            i2 = 1; j2 = 1; k2 = 0;
        } else if (x0 >= z0) {
            // X Z Y order
            i1 = 1; j1 = 0; k1 = 0;
            i2 = 1; j2 = 0; k2 = 1;
        } else {
            // Z X Y order
            i1 = 0; j1 = 0; k1 = 1;
            i2 = 1; j2 = 0; k2 = 1;
        }
    } else {
        // x0 < y0
        if (y0 < z0) {
            // Z Y X order
            i1 = 0; j1 = 0; k1 = 1;
            i2 = 0; j2 = 1; k2 = 1;
        } else if (x0 < z0) {
            // Y Z X order
            i1 = 0; j1 = 1; k1 = 0;
            i2 = 0; j2 = 1; k2 = 1;
        } else {
            // Y X Z order
            i1 = 0; j1 = 1; k1 = 0;
            i2 = 1; j2 = 1; k2 = 0;
        }
    }
    // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
    // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
    // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
    // c = 1/6.
    float32 x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
    float32 y1 = y0 - j1 + G3;
    float32 z1 = z0 - k1 + G3;
    float32 x2 = x0 - i2 + 2.0f * G3; // Offsets for third corner in (x,y,z) coords
    float32 y2 = y0 - j2 + 2.0f * G3;
    float32 z2 = z0 - k2 + 2.0f * G3;
    float32 x3 = x0 - 1.0f + 3.0f * G3; // Offsets for last corner in (x,y,z) coords
    float32 y3 = y0 - 1.0f + 3.0f * G3;
    float32 z3 = z0 - 1.0f + 3.0f * G3;
    // Work out the hashed gradient indices of the four simplex corners
    int ii = i & 255;
    int jj = j & 255;
    int kk = k & 255;
    int gi0 = permMod12[ii + perm[jj + perm[kk]]];
    int gi1 = permMod12[ii + i1 + perm[jj + j1 + perm[kk + k1]]];
    int gi2 = permMod12[ii + i2 + perm[jj + j2 + perm[kk + k2]]];
    int gi3 = permMod12[ii + 1 + perm[jj + 1 + perm[kk + 1]]];
    // Calculate the contribution from the four corners
    float32 t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
    if(t0 < 0) {
        n0 = 0.0f;
    } else {
        t0 *= t0;
        n0 = t0 * t0 * grad3[gi0].Dot(Vector3f(x0, y0, z0));
    }
    float32 t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
    if(t1 < 0) {
        n1 = 0.0f;
    } else {
        t1 *= t1;
        n1 = t1 * t1 * grad3[gi1].Dot(Vector3f(x1, y1, z1));
    }
    float32 t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
    if(t2 < 0) {
        n2 = 0.0f;
    } else {
        t2 *= t2;
        n2 = t2 * t2 * grad3[gi2].Dot(Vector3f(x2, y2, z2));
    }
    float32 t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
    if(t3 < 0) {
        n3 = 0.0f;
    } else {
        t3 *= t3;
        n3 = t3 * t3 * grad3[gi3].Dot(Vector3f(x3, y3, z3));
    }
    // Add contributions from each corner to get the final noise value.
    // The result is scaled to stay just inside [-1,1]
    return 32.0f * (n0 + n1 + n2 + n3);
}

float32 SimplexNoise::noise2DOctaves(float32 xin, float32 yin, float32 scale, int octaves) {
    const float32 gain = 0.5f;
    const float32 lacunarity = 2.0f;

    float32 total = 0.0f;
    float32 frequency = 1.0f / scale;
    float32 amplitude = 0.5f;

    for (int32 i = 0; i < octaves; ++i) {
        total += noise2D(xin * frequency, yin * frequency) * amplitude;         
        frequency *= lacunarity;
        amplitude *= gain;
    }
    return total;
}

float32 SimplexNoise::noise3DOctaves(float32 xin, float32 yin, float32 zin, float32 scale, int octaves) {
    const float32 gain = 0.5f;
    const float32 lacunarity = 2.0f;

    float32 total = 0.0f;
    float32 frequency = 1.0f / scale;
    float32 amplitude = 0.5f;

    for (int32 i = 0; i < octaves; ++i) {
        total += noise3D(xin * frequency, yin * frequency, zin * frequency) * amplitude;         
        frequency *= lacunarity;
        amplitude *= gain;
    }
    return total;
}

}
}