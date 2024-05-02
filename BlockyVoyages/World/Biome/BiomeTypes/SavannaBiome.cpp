#include "stdafx.h"

#include "SavannaBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

SavannaBiome::~SavannaBiome() {}

float32 SavannaBiome::GetHeight(const Vector2f& pos) const {
    return 0.0f;
}

Vector3f SavannaBiome::GetMapColor() const {
    return Vector3f(0.7f, 0.8f, 0.4f);
}

}
}
}