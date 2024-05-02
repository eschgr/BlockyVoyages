#include "stdafx.h"

#include "Input.h"

#include "../World/Camera.h"

#include "../Debugability/Logger.h"

#pragma comment(lib, "dinput8.lib")
#pragma comment(lib, "dxguid.lib")

namespace BlockyVoyages {
namespace Input {

using namespace Debugability;

InputProcessor* g_input = nullptr;

InputProcessor::InputProcessor(HINSTANCE hInst, HWND hwnd)
    : m_hInst(hInst),
        m_hWnd(hwnd),
        m_directInput8(nullptr),
        m_keyboardDevice(nullptr),
        m_mouseDevice(nullptr),
        m_keyBufferSize(10),
        m_butsActive(MAX_BUTTONS, 0),
        m_butsChanged(MAX_BUTTONS, 0),
        m_controller(nullptr) {
    if (FAILED(DirectInput8Create(m_hInst, 
                                  DIRECTINPUT_VERSION, 
                                  IID_IDirectInput8, 
                                  (void**)&m_directInput8, 
                                  nullptr))) 
    {
        LOG(Debugability::DBG_LVL_CRITICAL, "Unable to create Direct Input Context.");
        return;
    }

    if (FAILED(m_directInput8->CreateDevice(GUID_SysKeyboard, &m_keyboardDevice, nullptr))) {
        LOG(Debugability::DBG_LVL_CRITICAL,
            "Unable to create keyboard device.");
        return;
    }

    if (FAILED(m_keyboardDevice->SetDataFormat(&c_dfDIKeyboard))) {
        LOG(DBG_LVL_CRITICAL, "Unable to set keyboard data format.");
        return;
    }

    if (FAILED(m_keyboardDevice->SetCooperativeLevel(m_hWnd, DISCL_FOREGROUND | DISCL_NONEXCLUSIVE))) {
        LOG(DBG_LVL_CRITICAL, "Unable to set keyboard cooperative level.");
        return;
    }

    if (!setKeyboardBufferSize(m_keyBufferSize)) {
        LOG(DBG_LVL_CRITICAL, "Failed to set keyboard buffer size.");
        return;
    }
    
    if (FAILED(m_directInput8->CreateDevice(GUID_SysMouse, &m_mouseDevice, nullptr))) {
        LOG(DBG_LVL_CRITICAL, "Creation of main mouse failed.");
        return;
    }

    if (FAILED(m_mouseDevice->SetCooperativeLevel(m_hWnd, DISCL_FOREGROUND | DISCL_EXCLUSIVE))) {
        LOG(DBG_LVL_CRITICAL, "Unable to set main mouse cooperative level.");
        return;
    }

    if (FAILED(m_mouseDevice->SetDataFormat(&c_dfDIMouse))) {
        LOG(DBG_LVL_CRITICAL, "Unable to set data format for the main mouse.");
        return;
    }


    // 256 for the keys on the keyboard
    // the mouse can add up to 8 more, so take that into account
    m_butsActive.assign(MAX_BUTTONS, 0);

    m_keyboardDevice->Acquire();
    m_mouseDevice->Acquire();
}

InputProcessor::~InputProcessor(void) {
    // release the keyboard
    m_keyboardDevice->Unacquire();
    m_keyboardDevice->Release();

    m_mouseDevice->Unacquire();
    m_mouseDevice->Release();

    m_directInput8->Release();
    m_directInput8 = nullptr;
}

void InputProcessor::process() {
    processKeyboard();
    processMouse();

    if (!m_controller) {
        return;
    }

    // figure out what to do with the data gotten
    static const float32 rotSpeed = DegToRad(90.0f);

    float32 moveSpeed = 1.4f;

    // player controller processing
    if (m_butsChanged[DIK_LSHIFT] || m_butsChanged[DIK_RSHIFT]) {
        m_controller->toggleRun();
    }
    if (m_butsChanged[DIK_LCONTROL]) {
        m_controller->toggleCrouch();
    }
    if (m_butsActive[DIK_SPACE]) {
        m_controller->jump();
    }
    if (mouseButtonDown(0)) {
        m_controller->useRightHand();
    }
    if (mouseButtonDown(1)) {
        m_controller->useLeftHand();
    }
 
    // Turn the body to the sides
    m_controller->turnRightLeft(-m_mouseDX / 100.0f);

    // Move the head up and down
    m_controller->lookUpDown(-m_mouseDY / 100.0f);

    Vector3f moveDir(0.0f, 0.0f, 0.0f);
    if (m_butsActive[DIK_W]) {
        moveDir.z = -1.0f;
    } else if (m_butsActive[DIK_S]) {
        moveDir.z = 1.0f;
    }
    if (m_butsActive[DIK_A]) {
        moveDir.x = -1.0f;
    } else if (m_butsActive[DIK_D]) {
        moveDir.x = 1.0f;
    }
    m_controller->setTravelDirection(moveDir);

    // game state processing
    if (m_butsActive[DIK_ESCAPE]) {
        PostQuitMessage(0);
    }
}

void InputProcessor::reacquire(void) {
    if (m_keyboardDevice) {
        m_keyboardDevice->Acquire();
    }

    if (m_mouseDevice) {
        m_mouseDevice->Acquire();
    }
}

void InputProcessor::unacquire(void) {
    if (m_keyboardDevice) {
        m_keyboardDevice->Unacquire();
    }

    if (m_mouseDevice) {
        m_mouseDevice->Unacquire();
    }
}

bool InputProcessor::setKeyboardBufferSize(int32 bufferSize) {
    DIPROPDWORD            dipdw;

    dipdw.diph.dwSize = sizeof(DIPROPDWORD);
    dipdw.diph.dwHeaderSize = sizeof(DIPROPHEADER);
    dipdw.diph.dwObj = 0;
    dipdw.diph.dwHow = DIPH_DEVICE;
    dipdw.dwData = bufferSize;
        
    // make an attempt to change the buffer size
    m_keyboardDevice->Unacquire();
    HRESULT hr = m_keyboardDevice->SetProperty(DIPROP_BUFFERSIZE, &dipdw.diph);
    m_keyboardDevice->Acquire();
        
    if (FAILED(hr)) {
        // LOG(DBG_LVL_WARNING, "Failed to change the keyboard buffer size.");
        return false;
    }

    m_keyBufferSize = bufferSize;
    return true;
}

void InputProcessor::pushKeystroke(ubyte key) {
    m_keystrokes.push(key);
}

ubyte InputProcessor::popKeystroke(void) {
    if (!m_keystrokes.empty()) {
        ubyte key = m_keystrokes.back();
        m_keystrokes.pop();
        return key;
    }
    return 0;
}

void InputProcessor::processKeyboard(void) {
    char      buffer[256]; 
    HRESULT  hr; 
    hr = m_keyboardDevice->GetDeviceState(sizeof(buffer), (LPVOID)&buffer); 
    if (FAILED(hr)) { 
        // If it failed, the device has probably been lost. 
        m_keyboardDevice->Acquire();
        if (FAILED(m_keyboardDevice->GetDeviceState(sizeof(buffer), (LPVOID)&buffer))) {
		// LOG(DBG_LVL_WARNING, "Unable to reacquire the keyboard state on processing.");
            return;
        }
    }

    // find which keys are still down
    for(int32 i = 0; i < 256; ++i) {
        bool butActive = (buffer[i] & 0x80) != 0;
        m_butsChanged[i] = butActive != m_butsActive[i];
        m_butsActive[i] = butActive;
    }
}

void InputProcessor::processMouse(void) {
    DIMOUSESTATE mouseState;
    HRESULT  hr; 
    hr = m_mouseDevice->GetDeviceState(sizeof(DIMOUSESTATE), (LPVOID)&mouseState); 
    if (FAILED(hr)) {
        // If it failed, the device has probably been lost. 
        m_mouseDevice->Acquire();
        // if we still failed, just return
        if (FAILED(m_mouseDevice->GetDeviceState(sizeof(DIMOUSESTATE), (LPVOID)&mouseState))) {
            // LOG(DBG_LVL_WARNING, "Unable to reacquire the mouse state on processing.");
            return;
        }
    }

    // find which keys are still down
    for(int32 i = 0; i < 4; ++i) {
        bool butActive = (mouseState.rgbButtons[i] & 0x80) != 0;
        int32 button = 256 + i;
        m_butsChanged[button] = butActive != m_butsActive[button];
        m_butsActive[button] = butActive;
    }

    m_mouseX += mouseState.lX;
    m_mouseY += mouseState.lY;

    m_mouseDX = mouseState.lX;
    m_mouseDY = mouseState.lY;
}
}
}