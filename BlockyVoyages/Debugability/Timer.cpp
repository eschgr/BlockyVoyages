#include "stdafx.h"

#include "Timer.h"

#include <Windows.h>

namespace BlockyVoyages {
namespace Debugability {

Timer::Timer() {
    //query performance frequency
    LARGE_INTEGER perfInq;
    // inquire about the performance frequency
    QueryPerformanceFrequency(&perfInq);
    m_perfFreq = 1.0f / static_cast<float32>(perfInq.QuadPart);

    //init class vars
    m_startTick = GetCurrentTick();
}

Timer::~Timer() {}

float32 Timer::GetTime() {
    uint64 curTick = GetCurrentTick();
    return (curTick - m_startTick) * m_perfFreq;
}

uint64 Timer::GetCurrentTick() {
    LARGE_INTEGER perfCnt;
    QueryPerformanceCounter(&perfCnt);
    return perfCnt.QuadPart;
}

}
}