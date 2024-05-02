#include "stdafx.h"

#include "GrasslandBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

GrasslandBiome::~GrasslandBiome() {}

float32 GrasslandBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f GrasslandBiome::GetMapColor() const {
    return Vector3f(0.55f, 0.75f, 0.0f);
}

}
}
}