#pragma once

#include "../types.h"
#include "../Graphics/Font.h"
#include "Timer.h"

namespace BlockyVoyages {
namespace Debugability {

class PerformanceMonitor {
public:
    PerformanceMonitor(Timer* timer);
    ~PerformanceMonitor();

    void startFrame();
    void endFrame();

    void ResetStats();

    float32 getPrevFrameTime() { return m_lastFrameTime; }
    float32 getAverageFrameTime() { return m_avgFrameTime; }
    float32 getMovingAverageFrameTime() { return m_movingAverageTime; }
    float32 getMinFrameTime() { return m_minFrameTime; }
    float32 getMaxFrameTime() { return m_maxFrameTime; }
    int64 getFrameCount() { return m_numFrames; }

    void drawStats(int32 x, int32 y, Graphics::Font* font);
private:
    Timer* m_timer;
    float32 m_frameStartTime;

    // Frame statistics
    float32 m_lastFrameTime;
    float32 m_avgFrameTime;
    float32 m_minFrameTime;
    float32 m_maxFrameTime;
    int64 m_numFrames;

    static const int32 kNumWindowFrames = 60;
    float32 m_lastWindowFrameTimes[kNumWindowFrames];
    float32 m_movingAverageTime;
    int32 m_currentWindowFrame;
};
}
}