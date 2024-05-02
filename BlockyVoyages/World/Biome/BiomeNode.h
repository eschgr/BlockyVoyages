#if !defined(__BIOME_NODE_H__)
#define __BIOME_NODE_H__

#include "../../Graph/Node.h"


namespace BlockyVoyages {
namespace World {
namespace Biome {

class BiomeNode : public Graph::Node {
public:
    BiomeNode(const Graph::Node* other)
        : Graph::Node(other),
          m_height(0.0f),
          m_dist_to_shore(-1.0f),
          m_temperature(0.0f),
          m_moisture(0.0f)
    {}

    float32 GetShoreDistance() const { return m_dist_to_shore; }
    float32 GetHeight() const { return m_height; }
    float32 GetTemperature() const { return m_temperature; }
    float32 GetMoisture() const { return m_moisture; }

    void SetShoreDistance(float32 distance) { m_dist_to_shore = distance; }
    void SetHeight(float32 height) { m_height = height; }
    void SetTemperature(float32 temp) { m_temperature = temp; }
    void SetMoisture(float32 moisture) { m_moisture = moisture; }

    // returns true if the moisture of this node changed
    bool CalculateMoisture(const Vector2f& wind_dir);

    bool IsShore() const;
    bool IsLand() const;
    bool IsOcean() const;
private:
    float32 m_height;
    float32 m_dist_to_shore;
    float32 m_temperature;
    float32 m_moisture;
};

}
}
}

#endif