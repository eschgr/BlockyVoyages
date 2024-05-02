#include "stdafx.h"

#include "SnowBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

SnowBiome::~SnowBiome() {}

float32 SnowBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f SnowBiome::GetMapColor() const {
    return Vector3f(1.0f, 1.0f, 1.0f);
}

}
}
}