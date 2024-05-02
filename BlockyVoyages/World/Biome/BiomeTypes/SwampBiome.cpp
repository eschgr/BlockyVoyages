#include "stdafx.h"

#include "SwampBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

SwampBiome::~SwampBiome() {}

float32 SwampBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f SwampBiome::GetMapColor() const {
    return Vector3f(0.3f, 0.6f, 0.4f);
}

}
}
}