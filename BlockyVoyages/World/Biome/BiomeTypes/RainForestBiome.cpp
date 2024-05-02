#include "stdafx.h"

#include "RainForestBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

RainForestBiome::~RainForestBiome() {}

float32 RainForestBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f RainForestBiome::GetMapColor() const {
    return Vector3f(0.2f, 0.5f, 0.1f);
}

}
}
}