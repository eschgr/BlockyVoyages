#include "stdafx.h"

#include "DeepOceanBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

DeepOceanBiome::~DeepOceanBiome() {}

float32 DeepOceanBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f DeepOceanBiome::GetMapColor() const {
    return Vector3f(0.1f, 0.0f, 0.5f);
}

}
}
}