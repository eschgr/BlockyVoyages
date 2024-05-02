#include "stdafx.h"

#include "ScorchedBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

ScorchedBiome::~ScorchedBiome() {}

float32 ScorchedBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f ScorchedBiome::GetMapColor() const {
    return Vector3f(0.4f, 0.35f, 0.35f);
}

}
}
}