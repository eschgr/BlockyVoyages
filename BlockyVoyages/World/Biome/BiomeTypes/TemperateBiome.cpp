#include "stdafx.h"

#include "TemperateBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

TemperateBiome::~TemperateBiome() {}

float32 TemperateBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f TemperateBiome::GetMapColor() const {
    return Vector3f(0.25f, 0.65f, 0.1f);
}

}
}
}