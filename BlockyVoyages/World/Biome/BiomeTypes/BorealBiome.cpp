#include "stdafx.h"

#include "BorealBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

BorealBiome::~BorealBiome() {}

float32 BorealBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f BorealBiome::GetMapColor() const {
    return Vector3f(0.15f, 0.6f, 0.3f);
}

}
}
}