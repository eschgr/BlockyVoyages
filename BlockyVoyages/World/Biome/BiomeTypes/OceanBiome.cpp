#include "stdafx.h"

#include "OceanBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

OceanBiome::~OceanBiome() {}

float32 OceanBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f OceanBiome::GetMapColor() const {
    return Vector3f(0.0f, 0.5f, 0.75f);
}

}
}
}