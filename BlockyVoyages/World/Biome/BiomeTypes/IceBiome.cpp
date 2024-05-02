#include "stdafx.h"

#include "IceBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

IceBiome::~IceBiome() {}

float32 IceBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f IceBiome::GetMapColor() const {
    return Vector3f(0.9f, 0.9f, 1.0f);
}

}
}
}