#pragma once

#include "../types.h"

#include <queue>

#include <Windows.h>
#include <DInput.h>

#include "Controller.h"

namespace BlockyVoyages {
namespace Input {

// This class manages all the input and forwards it to the appropriate objects.
class InputProcessor {
public:
    InputProcessor(HINSTANCE hInst, HWND hwnd);
    ~InputProcessor();

    void process();
    void reacquire();
    void unacquire();

    bool setKeyboardBufferSize(int32 bufferSize);

    void pushKeystroke(ubyte key);
    ubyte popKeystroke();

    // sets the object controller movement keystrokes will be sent to
    void setController(Controller* controller) { m_controller = controller; }

    bool keyDown(ubyte key) { return m_butsActive[key]; }
    bool keyChanged(ubyte key) { return m_butsChanged[key]; }

private:
    // 256 keyboard buttons and 8 mouse buttons
    static const int MAX_BUTTONS = 256 + 8;

    HINSTANCE m_hInst;
    HWND m_hWnd;
    LPDIRECTINPUT8 m_directInput8;

    // some input stuff
    int32 m_keyBufferSize;

    // devices
    LPDIRECTINPUTDEVICE8 m_keyboardDevice;
    LPDIRECTINPUTDEVICE8 m_mouseDevice;

    std::vector<bool> m_butsActive;
    std::vector<bool> m_butsChanged;

    int64 m_mouseX;
    int64 m_mouseY;

    int64 m_mouseDX;
    int64 m_mouseDY;

    std::queue<ubyte> m_keystrokes;

    // The object that input commands will be sent to
    Controller* m_controller;

    void processKeyboard();
    void processMouse();

    inline bool mouseButtonDown(int32 button) const { return m_butsActive[256 + button]; }
};

extern InputProcessor* g_input;

}
}