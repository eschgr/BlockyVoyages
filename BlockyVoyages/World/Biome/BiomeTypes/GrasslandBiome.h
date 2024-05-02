#if !defined(__GRASSLAND_BIOME_H__)
#define __GRASSLAND_BIOME_H__

#include "../Biome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

class TreeGenerator;

class GrasslandBiome : public Biome {
public:
    virtual ~GrasslandBiome();

    float32 GetHeight(const Vector2f& pos) const;
    // GetSurfaceType: Get the surface point at a specific spot. May change
    // across the biome
    BlockType GetSurfaceType(const Vector2f& pos) const { return kGrass; }
    
    // GetPrimarySurfaceType: Gets the main surface type for the biome. This
    // is the type that would be expected to be seen at a distance.
    BlockType GetPrimarySurfaceType() const { return kGrass; }

    Vector3f GetMapColor() const;
};

}  // namespace BlockyVoyages
}  // namespace World
}  // namespace Biome

#endif  // __GRASSLAND_BIOME_H__