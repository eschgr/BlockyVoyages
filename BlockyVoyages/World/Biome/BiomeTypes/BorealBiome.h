#if !defined(__BOREAL_BIOME_H__)
#define __BOREAL_BIOME_H__

#include "../Biome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

class TreeGenerator;

class BorealBiome : public Biome {
public:
    virtual ~BorealBiome();

    float32 GetHeight(const Vector2f& pos) const;
    
    // GetSurfaceType: Get the surface point at a specific spot. May change
    // across the biome
    BlockType GetSurfaceType(const Vector2f& pos) const { return kSnowyGrass; }
    
    // GetPrimarySurfaceType: Gets the main surface type for the biome. This
    // is the type that would be expected to be seen at a distance.
    BlockType GetPrimarySurfaceType() const { return kSnowyGrass; }

    Vector3f GetMapColor() const;
};

}  // namespace BlockyVoyages
}  // namespace World
}  // namespace Biome

#endif  // __BOREAL_BIOME_H__