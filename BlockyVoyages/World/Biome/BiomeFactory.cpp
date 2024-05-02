#include "stdafx.h"

#include "BiomeFactory.h"
#include "BiomeTypes/BorealBiome.h"
#include "BiomeTypes/DeepOceanBiome.h"
#include "BiomeTypes/DesertBiome.h"
#include "BiomeTypes/GrasslandBiome.h"
#include "BiomeTypes/IceBiome.h"
#include "BiomeTypes/OceanBiome.h"
#include "BiomeTypes/RainForestBiome.h"
#include "BiomeTypes/SavannaBiome.h"
#include "BiomeTypes/ScorchedBiome.h"
#include "BiomeTypes/SnowBiome.h"
#include "BiomeTypes/SwampBiome.h"
#include "BiomeTypes/TemperateBiome.h"

namespace BlockyVoyages {
namespace World {
namespace Biome {

Biome* GetDeepOceanBiome() {
    return new DeepOceanBiome();
}

Biome* GetOceanBiome() {
    return new OceanBiome();
}

Biome* GetDesertBiome() {
    return new DesertBiome();
}

Biome* GetIceBiome() {
    return new IceBiome();
}

Biome* GetSnowBiome() {
    return new SnowBiome();
}

Biome* GetBorealBiome() {
    return new BorealBiome();
}

Biome* GetGrasslandBiome() {
    return new GrasslandBiome();
}

Biome* GetTemperateBiome() {
    return new TemperateBiome();
}

Biome* GetSwampBiome() {
    return new SwampBiome();
}

Biome* GetSavannaBiome() {
    return new SavannaBiome();
}

Biome* GetRainForestBiome() {
    return new RainForestBiome();
}

Biome* GetScorchedBiome() {
    return new ScorchedBiome();
}

BiomeFactory::BiomeFactory() {
    m_biome_map[BiomePair(0, 0)] = GetDesertBiome;
    m_biome_map[BiomePair(0, 1)] = GetIceBiome;
    m_biome_map[BiomePair(0, 2)] = GetIceBiome;
    m_biome_map[BiomePair(0, 3)] = GetIceBiome;
    m_biome_map[BiomePair(0, 4)] = GetIceBiome;
    m_biome_map[BiomePair(0, 5)] = GetSnowBiome;
    m_biome_map[BiomePair(0, 6)] = GetSnowBiome;
    m_biome_map[BiomePair(0, 7)] = GetSnowBiome;
    m_biome_map[BiomePair(0, 8)] = GetSnowBiome;
    m_biome_map[BiomePair(0, 9)] = GetSnowBiome;

    m_biome_map[BiomePair(1, 0)] = GetDesertBiome;
    m_biome_map[BiomePair(1, 1)] = GetIceBiome;
    m_biome_map[BiomePair(1, 2)] = GetIceBiome;
    m_biome_map[BiomePair(1, 3)] = GetIceBiome;
    m_biome_map[BiomePair(1, 4)] = GetIceBiome;
    m_biome_map[BiomePair(1, 5)] = GetSnowBiome;
    m_biome_map[BiomePair(1, 6)] = GetSnowBiome;
    m_biome_map[BiomePair(1, 7)] = GetSnowBiome;
    m_biome_map[BiomePair(1, 8)] = GetSnowBiome;
    m_biome_map[BiomePair(1, 9)] = GetSnowBiome;

    m_biome_map[BiomePair(2, 0)] = GetDesertBiome;
    m_biome_map[BiomePair(2, 1)] = GetSnowBiome;
    m_biome_map[BiomePair(2, 2)] = GetBorealBiome;
    m_biome_map[BiomePair(2, 3)] = GetBorealBiome;
    m_biome_map[BiomePair(2, 4)] = GetBorealBiome;
    m_biome_map[BiomePair(2, 5)] = GetBorealBiome;
    m_biome_map[BiomePair(2, 6)] = GetSnowBiome;
    m_biome_map[BiomePair(2, 7)] = GetSnowBiome;
    m_biome_map[BiomePair(2, 8)] = GetSnowBiome;
    m_biome_map[BiomePair(2, 9)] = GetSnowBiome;

    m_biome_map[BiomePair(3, 0)] = GetDesertBiome;
    m_biome_map[BiomePair(3, 1)] = GetDesertBiome;
    m_biome_map[BiomePair(3, 2)] = GetBorealBiome;
    m_biome_map[BiomePair(3, 3)] = GetBorealBiome;
    m_biome_map[BiomePair(3, 4)] = GetBorealBiome;
    m_biome_map[BiomePair(3, 5)] = GetBorealBiome;
    m_biome_map[BiomePair(3, 6)] = GetBorealBiome;
    m_biome_map[BiomePair(3, 7)] = GetBorealBiome;
    m_biome_map[BiomePair(3, 8)] = GetSnowBiome;
    m_biome_map[BiomePair(3, 9)] = GetSnowBiome;

    m_biome_map[BiomePair(4, 0)] = GetDesertBiome;
    m_biome_map[BiomePair(4, 1)] = GetDesertBiome;
    m_biome_map[BiomePair(4, 2)] = GetGrasslandBiome;
    m_biome_map[BiomePair(4, 3)] = GetTemperateBiome;
    m_biome_map[BiomePair(4, 4)] = GetTemperateBiome;
    m_biome_map[BiomePair(4, 5)] = GetTemperateBiome;
    m_biome_map[BiomePair(4, 6)] = GetTemperateBiome;
    m_biome_map[BiomePair(4, 7)] = GetTemperateBiome;
    m_biome_map[BiomePair(4, 8)] = GetTemperateBiome;
    m_biome_map[BiomePair(4, 9)] = GetTemperateBiome;

    m_biome_map[BiomePair(5, 0)] = GetDesertBiome;
    m_biome_map[BiomePair(5, 1)] = GetDesertBiome;
    m_biome_map[BiomePair(5, 2)] = GetGrasslandBiome;
    m_biome_map[BiomePair(5, 3)] = GetGrasslandBiome;
    m_biome_map[BiomePair(5, 4)] = GetGrasslandBiome;
    m_biome_map[BiomePair(5, 5)] = GetTemperateBiome;
    m_biome_map[BiomePair(5, 6)] = GetTemperateBiome;
    m_biome_map[BiomePair(5, 7)] = GetTemperateBiome;
    m_biome_map[BiomePair(5, 8)] = GetTemperateBiome;
    m_biome_map[BiomePair(5, 9)] = GetSwampBiome;

    m_biome_map[BiomePair(6, 0)] = GetDesertBiome;
    m_biome_map[BiomePair(6, 1)] = GetDesertBiome;
    m_biome_map[BiomePair(6, 2)] = GetGrasslandBiome;
    m_biome_map[BiomePair(6, 3)] = GetGrasslandBiome;
    m_biome_map[BiomePair(6, 4)] = GetGrasslandBiome;
    m_biome_map[BiomePair(6, 5)] = GetTemperateBiome;
    m_biome_map[BiomePair(6, 6)] = GetTemperateBiome;
    m_biome_map[BiomePair(6, 7)] = GetRainForestBiome;
    m_biome_map[BiomePair(6, 8)] = GetRainForestBiome;
    m_biome_map[BiomePair(6, 9)] = GetSwampBiome;

    m_biome_map[BiomePair(7, 0)] = GetDesertBiome;
    m_biome_map[BiomePair(7, 1)] = GetDesertBiome;
    m_biome_map[BiomePair(7, 2)] = GetSavannaBiome;
    m_biome_map[BiomePair(7, 3)] = GetSavannaBiome;
    m_biome_map[BiomePair(7, 4)] = GetSavannaBiome;
    m_biome_map[BiomePair(7, 5)] = GetSavannaBiome;
    m_biome_map[BiomePair(7, 6)] = GetRainForestBiome;
    m_biome_map[BiomePair(7, 7)] = GetRainForestBiome;
    m_biome_map[BiomePair(7, 8)] = GetRainForestBiome;
    m_biome_map[BiomePair(7, 9)] = GetSwampBiome;

    m_biome_map[BiomePair(8, 0)] = GetScorchedBiome;
    m_biome_map[BiomePair(8, 1)] = GetDesertBiome;
    m_biome_map[BiomePair(8, 2)] = GetDesertBiome;
    m_biome_map[BiomePair(8, 3)] = GetSavannaBiome;
    m_biome_map[BiomePair(8, 4)] = GetSavannaBiome;
    m_biome_map[BiomePair(8, 5)] = GetRainForestBiome;
    m_biome_map[BiomePair(8, 6)] = GetRainForestBiome;
    m_biome_map[BiomePair(8, 7)] = GetRainForestBiome;
    m_biome_map[BiomePair(8, 8)] = GetRainForestBiome;
    m_biome_map[BiomePair(8, 9)] = GetRainForestBiome;

    m_biome_map[BiomePair(9, 0)] = GetScorchedBiome;
    m_biome_map[BiomePair(9, 1)] = GetScorchedBiome;
    m_biome_map[BiomePair(9, 2)] = GetScorchedBiome;
    m_biome_map[BiomePair(9, 3)] = GetDesertBiome;
    m_biome_map[BiomePair(9, 4)] = GetSavannaBiome;
    m_biome_map[BiomePair(9, 5)] = GetRainForestBiome;
    m_biome_map[BiomePair(9, 6)] = GetRainForestBiome;
    m_biome_map[BiomePair(9, 7)] = GetRainForestBiome;
    m_biome_map[BiomePair(9, 8)] = GetRainForestBiome;
    m_biome_map[BiomePair(9, 9)] = GetRainForestBiome;
}

Biome* BiomeFactory::GetBiome(float32 moisture, float32 temperature) const {
    int32 moisture_ind = static_cast<int32>(10.0f * (moisture - m_moisture_bounds.s) / m_moisture_scale);
    int32 temp_ind = static_cast<int32>(10.0f * (temperature - m_temperature_bounds.s) / m_temperature_scale);

    moisture_ind = max(min(moisture_ind, 9), 0);
    temp_ind = max(min(temp_ind, 9), 0);

    auto biome_gen_func = m_biome_map.find(BiomePair(temp_ind, moisture_ind));
    assert(biome_gen_func != m_biome_map.end());
    return biome_gen_func->second();
}

Biome* BiomeFactory::GetSeaBiome(bool deep_ocean) const {
    if (deep_ocean) {
        return GetDeepOceanBiome();
    } else {
        return GetOceanBiome();
    }
}

}  // BlockyVoyages
}  // World
}  // Biome