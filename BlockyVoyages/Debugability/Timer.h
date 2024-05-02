#pragma once

#include "../types.h"

namespace BlockyVoyages {
namespace Debugability {

class Timer {
public:
    Timer();
    ~Timer();

    float32 GetTime();
private:
    uint64 m_startTick;
    float32 m_perfFreq;

    uint64 GetCurrentTick();
};
}
}