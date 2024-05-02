#if !defined(__DEEP_OCEAN_BIOME_H__)
#define __DEEP_OCEAN_BIOME_H__

#include "../Biome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

class TreeGenerator;

class DeepOceanBiome : public Biome {
public:
    virtual ~DeepOceanBiome();

    float32 GetHeight(const Vector2f& pos) const;
    
    // GetSurfaceType: Get the surface point at a specific spot. May change
    // across the biome
    BlockType GetSurfaceType(const Vector2f& pos) const { return kDirt; }
    
    // GetPrimarySurfaceType: Gets the main surface type for the biome. This
    // is the type that would be expected to be seen at a distance.
    BlockType GetPrimarySurfaceType() const { return kDirt; }

    Vector3f GetMapColor() const;
};

}  // namespace BlockyVoyages
}  // namespace World
}  // namespace Biome

#endif  // __DEEP_OCEAN_BIOME_H__