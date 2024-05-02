#if !defined(__BIOME_FACTORY_H__)
#define __BIOME_FACTORY_H__

#include "../../types.h"
#include "../../Math/MathMain.h"

#include <map>

namespace BlockyVoyages {
namespace World {
namespace Biome {

class Biome;

class BiomeFactory {
public:
    BiomeFactory();
    Biome* GetBiome(float32 moisture, float32 temperature) const;
    Biome* GetSeaBiome(bool deep_ocean) const;

    void SetTemperatureBounds(const Vector2f& bounds) {
        m_temperature_bounds = bounds;
        m_temperature_scale = bounds.y - bounds.x;
    }
    void SetMoistureBounds(const Vector2f& bounds) {
        m_moisture_bounds = bounds;
        m_moisture_scale = bounds.y - bounds.x;
    }
private:
    typedef std::pair<int32, int32> BiomePair;
    std::map<std::pair<int32, int32>, Biome* (*)()> m_biome_map;

    Vector2f m_temperature_bounds;
    float32 m_temperature_scale;
    Vector2f m_moisture_bounds;
    float32 m_moisture_scale;
};

}  // namespace Biome
}  // namespace World
}  // namespace BlockyVoyages

#endif  // __BIOME_FACTORY_H__