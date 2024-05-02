#pragma once

#include "resource.h"

namespace BlockyVoyages {

class MainWindow {
public:
    MainWindow(HINSTANCE hInstance, int cmdShow, int32 width, int32 height, int32 color_depth, bool fullscreen, const TCHAR* winTitle);
    ~MainWindow();

    int32 GetWidth() const { return m_width; }
    int32 GetHeight() const { return m_height; }
    bool IsFullscreen() const { return m_fullscreen; }

    HINSTANCE GetHInstance() const { return m_hInst; }
    HWND GetHWindow() const { return m_hWnd; }

private:
    int32 m_width;
    int32 m_height;
    int32 m_colorDepth;
    bool m_fullscreen;
    TCHAR* m_winTitle;

    HWND m_hWnd;
    HINSTANCE m_hInst;
    TCHAR m_winClass[100];

    ATOM MyRegisterClass();
    bool InitInstance(int cmdShow);
};

extern MainWindow* g_mainWindow;

}