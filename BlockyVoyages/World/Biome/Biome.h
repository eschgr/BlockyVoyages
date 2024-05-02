#if !defined(__BIOME_H__)
#define __BIOME_H__

#include "../../types.h"
#include "../../Math/MathMain.h"

#include "../BlockTypes.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

class TreeGenerator;

class Biome {
public:
    virtual ~Biome() = 0 {}

    // GetHeight returns the offset from the base height. This must be added to
    // the region height to get the world height.
    virtual float32 GetHeight(const Vector2f& pos) const = 0;
    virtual TreeGenerator* GetTreeGenerator() { return nullptr; }
    
    // GetSurfaceType: Get the surface point at a specific spot. May change
    // across the biome
    virtual BlockType GetSurfaceType(const Vector2f& pos) const = 0;
    
    // GetPrimarySurfaceType: Gets the main surface type for the biome. This
    // is the type that would be expected to be seen at a distance.
    virtual BlockType GetPrimarySurfaceType() const = 0;
    virtual Vector3f GetMapColor() const = 0;
};

}  // namespace BlockyVoyages
}  // namespace World
}  // namespace Biome

#endif  // !defined(__BIOME_H__)