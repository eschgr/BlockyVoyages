#pragma once

#include "MathMain.h"

namespace BlockyVoyages {
namespace Math {

/**
  * A helper class used to calculate random noise using the simplex noise
  * function. This class is entirely stateless
  */
class SimplexNoise {
public:
    static float32 noise2D(float32 xin, float32 yin);
    static float32 noise3D(float32 xin, float32 yin, float32 zin);
    static float32 noise2DOctaves(float32 xin, float32 yin, float32 scale, int octaves);
    static float32 noise3DOctaves(float32 xin, float32 yin, float32 zin, float32 scale, int octaves);
private:
    static const Vector2f grad2[];
    static const Vector3f grad3[];
    static const uint8 perm[];
    static const uint8 permMod12[];
    static const float32 F2, F3;
    static const float32 G2, G3;
};

}
}