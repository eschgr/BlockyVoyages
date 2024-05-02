#if !defined(__BIOME_REGION_H__)
#define __BIOME_REGION_H

#include "../../Graph/Cell.h"

#include <map>

namespace BlockyVoyages {
namespace Graphics {
class VertexArrayObject;
class Shader;
}
namespace World {
namespace Biome {

// see my google drive for the mappings
/*enum BiomeType {
    DESERT,
    ICE,
    SNOW,
    GRASSLAND,
    BOREAL,
    SAVANNA,
    TEMPERATE,
    RAINFOREST,
    SCORCHED,
    SWAMP,
    OCEAN,
    LAKE
};*/

class BiomeFactory;
class Biome;

// std::map<std::pair<int32, int32>, BiomeType> TempMoistureBiomeMap;

class BiomeRegion : public Graph::Cell {
public:
    BiomeRegion(const Cell* other);
    ~BiomeRegion();

    bool IsWater() const { return m_is_water; }
    void SetIsWater() { m_is_water = true; }

    void CalculateCenterStats();
    void SetBiome(const BiomeFactory& biome_factory);
    const Biome* GetBiome() const { return m_biome; }

    float32 GetHeightAt(const Vector2f& point);

    void Draw(Graphics::Shader* shader);
private:
    bool m_is_water;
    float32 m_height;
    float32 m_temperature;
    float32 m_moisture;

    Biome* m_biome;

    Graphics::VertexArrayObject* m_vertex_array;
};

}
}
}

#endif