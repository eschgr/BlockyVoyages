#include "stdafx.h"

#include "PerfMon.h"
#include "Timer.h"

#include <sstream>

#include <assert.h>
#include <Windows.h>

namespace BlockyVoyages {
namespace Debugability {

PerformanceMonitor::PerformanceMonitor(Timer* timer)
    : m_timer(timer),
      m_frameStartTime(0.0f),
      m_lastFrameTime(0.0f),
      m_avgFrameTime(0.0f),
      m_minFrameTime(10000.0f), // 1000s should be sufficiently large. If not, there's a fundamental problem.
      m_maxFrameTime(0.0f),
      m_numFrames(0),
      m_movingAverageTime(0.0f),
      m_currentWindowFrame(0)
{
    memset(m_lastWindowFrameTimes, 0, sizeof(m_lastWindowFrameTimes));
}

PerformanceMonitor::~PerformanceMonitor() {}

void PerformanceMonitor::startFrame() {
    m_frameStartTime = m_timer->GetTime();
}

void PerformanceMonitor::endFrame() {
    // if the frame start time is 0, end was called without a start, so this
    // frame should be ignored.
    if (m_frameStartTime == 0.0f) {
        return;
    }

    // calculate the statistics
    float32 curTime = m_timer->GetTime();
    m_lastFrameTime = curTime - m_frameStartTime;
    ++m_numFrames;

    m_avgFrameTime *= static_cast<float32>(m_numFrames - 1) / m_numFrames;
    m_avgFrameTime += m_lastFrameTime / m_numFrames;

    if (m_lastFrameTime < m_minFrameTime) {
        m_minFrameTime = m_lastFrameTime;
    }
    if (m_lastFrameTime > m_maxFrameTime) {
        m_maxFrameTime = m_lastFrameTime;
    }

    // Calculate the moving average time
    m_movingAverageTime -= m_lastWindowFrameTimes[m_currentWindowFrame];
    m_lastWindowFrameTimes[m_currentWindowFrame] = m_lastFrameTime / kNumWindowFrames;
    m_movingAverageTime += m_lastWindowFrameTimes[m_currentWindowFrame];
    m_currentWindowFrame = (m_currentWindowFrame + 1) % kNumWindowFrames;

    // reset any per frame state
    m_frameStartTime = 0.0f;
}

void PerformanceMonitor::ResetStats() {
    m_avgFrameTime = 0.0f;
    m_minFrameTime = 10000.0f;
    m_maxFrameTime = 0.0f;
    m_numFrames = 0;

    memset(m_lastWindowFrameTimes, 0, sizeof(m_lastWindowFrameTimes));
    m_movingAverageTime = 0.0f;
    m_currentWindowFrame = 0;
}

void PerformanceMonitor::drawStats(int32 x, int32 y, Graphics::Font* font) {
    // print out some stats
    const Vector3f red(1.0f, 0.0f, 0.0f);
    std::stringstream strstream;
    strstream << "FPS: " << 1.0f / m_movingAverageTime;
    font->drawText(x, y, red, strstream.str());
    
    strstream.clear();
    strstream.str("");
    strstream << "Avg FPS: " << 1.0f / m_avgFrameTime;
    font->drawText(x, y + 20, red, strstream.str());

    strstream.clear();
    strstream.str("");
    strstream << "Max FPS: " << 1.0f / m_minFrameTime;
    font->drawText(x, y + 40, red, strstream.str());

    strstream.clear();
    strstream.str("");
    strstream << "Min FPS: " << 1.0f / m_maxFrameTime;
    font->drawText(x, y + 60, red, strstream.str());
}

}
}