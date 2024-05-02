#include "stdafx.h"

#include "MainWindow.h"

#include "Input\Input.h"

namespace BlockyVoyages {
MainWindow* g_mainWindow = nullptr;

//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE:  Processes messages for the main window.
//
//  WM_PAINT	- Paint the main window
//  WM_DESTROY	- post a quit message and return
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam,
                         LPARAM lParam) {
  switch (message) {
    case WM_DESTROY:
      PostQuitMessage(0);
      break;
    case WM_KEYDOWN:
      switch (wParam) {
        case VK_ESCAPE:
          PostQuitMessage(0);
          break;
        default:
          return DefWindowProc(hWnd, message, wParam, lParam);
      }
      break;
    case WM_ACTIVATE:
      if (Input::g_input) {
        if (wParam) {
          Input::g_input->reacquire();
        } else {
          Input::g_input->unacquire();
        }
      }
      break;
    default:
      return DefWindowProc(hWnd, message, wParam, lParam);
  }
  return 0;
}

MainWindow::MainWindow(HINSTANCE hInstance, int cmdShow, int32 width,
                       int32 height, int32 color_depth, bool fullscreen,
                       const TCHAR* winTitle)
    : m_hInst(hInstance),
      m_width(width),
      m_height(height),
      m_fullscreen(fullscreen),
      m_colorDepth(color_depth) {
  m_winTitle = new TCHAR[wcslen(winTitle) + 1];
  wcscpy(m_winTitle, winTitle);
  LoadString(m_hInst, IDC_BLOCKYVOYAGES, m_winClass,
             sizeof(m_winClass) / sizeof(m_winClass[0]));
  MyRegisterClass();
  InitInstance(cmdShow);
}

MainWindow::~MainWindow() { delete[] m_winTitle; }

//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MainWindow::MyRegisterClass() {
  WNDCLASSEX wcex;
  ZeroMemory(&wcex, sizeof(wcex));
  wcex.cbSize = sizeof(WNDCLASSEX);
  wcex.lpszClassName = m_winClass;
  wcex.cbSize = sizeof(WNDCLASSEX);
  wcex.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
  wcex.lpfnWndProc = WndProc;
  wcex.hInstance = m_hInst;
  wcex.hIcon = LoadIcon(m_hInst, MAKEINTRESOURCE(IDI_BLOCKYVOYAGES));
  wcex.hIconSm = LoadIcon(m_hInst, MAKEINTRESOURCE(IDI_SMALL));
  wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
  wcex.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
  wcex.lpszMenuName = nullptr;
  wcex.cbClsExtra = 0;
  wcex.cbWndExtra = 0;
  return RegisterClassEx(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
bool MainWindow::InitInstance(int cmdShow) {
  RECT clientArea;

  clientArea.left = clientArea.top = 0;
  clientArea.right = m_width;
  clientArea.bottom = m_height;

  if (m_fullscreen) {
    DEVMODE screen_settings;
    ZeroMemory(&screen_settings, sizeof(screen_settings));

    screen_settings.dmSize = sizeof(screen_settings);
    screen_settings.dmPelsWidth = m_width;
    screen_settings.dmPelsHeight = m_height;
    screen_settings.dmBitsPerPel = m_colorDepth;
    screen_settings.dmFields = DM_BITSPERPEL | DM_PELSWIDTH | DM_PELSHEIGHT;

    // Try To Set Selected Mode And Get Results.  NOTE: CDS_FULLSCREEN Gets Rid
    // Of Start Bar.
    LONG result = ChangeDisplaySettings(&screen_settings, CDS_FULLSCREEN);
    if (result != DISP_CHANGE_SUCCESSFUL) {
      m_fullscreen = false;
    }
  }

  DWORD dwExStyle;
  DWORD dwStyle;

  if (m_fullscreen) {
    dwExStyle = WS_EX_APPWINDOW;
    dwStyle = WS_POPUP;
    ShowCursor(FALSE);
  } else {
    dwExStyle = WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;
    dwStyle = WS_OVERLAPPEDWINDOW;
  }

  AdjustWindowRectEx(&clientArea, dwStyle, false, dwExStyle);

  clientArea.right -= clientArea.left;
  clientArea.bottom -= clientArea.top;
  if (m_fullscreen) {
    clientArea.left = clientArea.top = 0;
  } else {
    clientArea.left = clientArea.top = 100;
  }
  m_hWnd =
      CreateWindowEx(dwExStyle, m_winClass, m_winTitle,
                     dwStyle | WS_CLIPSIBLINGS | WS_CLIPCHILDREN,
                     clientArea.left, clientArea.top, clientArea.right,
                     clientArea.bottom, nullptr, nullptr, m_hInst, nullptr);

  if (0 == m_hWnd) {
    LOG(Debugability::DBG_LVL_ALWAYS_TERMINAL, "Failed to create main window.");
    return false;
  }
  ShowWindow(m_hWnd, cmdShow);
  UpdateWindow(m_hWnd);
  LOG(Debugability::DBG_LVL_NORMAL, "Main window created successfully.");
  return true;
}

}  // namespace BlockyVoyages