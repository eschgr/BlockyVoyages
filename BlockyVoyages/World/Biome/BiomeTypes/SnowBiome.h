#if !defined(__SNOW_BIOME_H__)
#define __SNOW_BIOME_H__

#include "../Biome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

class TreeGenerator;

class SnowBiome : public Biome {
public:
    virtual ~SnowBiome();

    float32 GetHeight(const Vector2f& pos) const;
    // GetSurfaceType: Get the surface point at a specific spot. May change
    // across the biome
    BlockType GetSurfaceType(const Vector2f& pos) const { return kSnow; }
    
    // GetPrimarySurfaceType: Gets the main surface type for the biome. This
    // is the type that would be expected to be seen at a distance.
    BlockType GetPrimarySurfaceType() const { return kSnow; }
    Vector3f GetMapColor() const;
};

}  // namespace BlockyVoyages
}  // namespace World
}  // namespace Biome

#endif  // __SNOW_BIOME_H__