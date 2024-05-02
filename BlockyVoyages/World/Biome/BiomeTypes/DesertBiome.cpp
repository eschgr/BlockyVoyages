#include "stdafx.h"

#include "DesertBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

DesertBiome::~DesertBiome() {}

float32 DesertBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f DesertBiome::GetMapColor() const {
    return Vector3f(1.0f, 0.95f, 0.5f);
}

}
}
}