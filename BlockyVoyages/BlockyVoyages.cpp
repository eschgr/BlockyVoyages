// BlockyVoyages.cpp : Defines the entry point for the application.
//

#include "stdafx.h"
#include "MainWindow.h"

#include "Debugability/Logger.h"
#include "Debugability/PerfMon.h"
#include "Debugability/Timer.h"

#include "Graphics/FontManager.h"
#include "Graphics/OpenGL.h"
#include "Graphics/Shader.h"
#include "Graphics/Texture.h"

#include "Input/Input.h"
#include "Input/CameraController.h"

#include "Math/MathMain.h"

#include "World/Character.h"
#include "World/Map.h"
#include "World/PerspectiveCamera.h"

#include <sstream>

// Constants
// TODO (gesch): Almost all of these should be in a configuration file instead.
const static int kMaxLoadString = 100;

// Graphics Features

const static bool kVsync = false;
const static int kColorDepth = 32;

#if 1
const static int kWidth = 1366;
const static int kHeight = 768;
const static bool kFullscreen = true;
#else
const static int kWidth = 640;
const static int kHeight = 480;
const static bool kFullscreen = false;
#endif

bool fullscreen = kFullscreen;

//*********************************************************
// World features
//*********************************************************
const static BlockyVoyages::Vector2f kRegionLow(-30.0f, -30.0f);
const static BlockyVoyages::Vector2f kRegionHigh(30.0f, 30.0f);

const static float32 kFeatureSize = 0.25f;
const static float32 kMaxHeight = 25.0f;
const static float32 kWorldArea = 250000.0f;  // in km
const static int32 kStartSeed = 69420;

bool enableOctree = false;
bool updateViewable = true;

BlockyVoyages::Vector3f playerStartPos(0.0, kMaxHeight + 10.0f, 0.0f);

using namespace BlockyVoyages;
using BlockyVoyages::Debugability::logger;

// Global Variables:
HINSTANCE hInst;								// current instance
TCHAR szTitle[kMaxLoadString];					// The title bar text
TCHAR szWindowClass[kMaxLoadString];			// the main window class name

World::Camera* mainCam = nullptr;

std::unique_ptr<World::Map> g_world_map;

void CreateWorld(int32 seed) {
    g_world_map.reset(new World::Map(kFeatureSize, kMaxHeight, kWorldArea, seed));

#if defined(_DEBUG)
    Vector2f region_low = kRegionLow;
    Vector2f region_high = kRegionHigh;
#else
    Vector2f half_world_dims = g_world_map->GetWorldDimensions() * 0.5f;
    Vector2f region_low(-half_world_dims.x, -half_world_dims.y);
    Vector2f region_high(half_world_dims.x, half_world_dims.y);
#endif
    g_world_map->addRegion(region_low, region_high);
}

// Forward declarations of functions included in this code module:
int APIENTRY _tWinMain(_In_ HINSTANCE hInstance,
                       _In_opt_ HINSTANCE hPrevInstance,
                       _In_ LPTSTR lpCmdLine,
                       _In_ int nCmdShow)
{
	UNREFERENCED_PARAMETER(hPrevInstance);
	UNREFERENCED_PARAMETER(lpCmdLine);

    // Initialize the logger
    logger = new Debugability::Logger();

 	MSG msg;
    memset(&msg, 0, sizeof(msg));

    int glError = 0;

    // Initialize the main window
    LoadString(hInstance, IDS_APP_TITLE, szTitle, kMaxLoadString);
    g_mainWindow = new MainWindow(hInstance, nCmdShow, kWidth, kHeight, kColorDepth, kFullscreen, szTitle);

    if (0 == g_mainWindow->GetHWindow()) {
		return 0;
	}

    // Create the OpenGL object
    Graphics::g_opengl = new Graphics::OpenGL(g_mainWindow->GetHWindow());
    if (!Graphics::g_opengl->isInitialized()) {
        return 1;
    }
    if (!kVsync) {
        wglSwapIntervalEXT(0);
    }
    BlockyVoyages::Graphics::Texture* texture = new BlockyVoyages::Graphics::Texture("../Data/grass.png", false);
    BlockyVoyages::Graphics::Shader* shader = new BlockyVoyages::Graphics::Shader();
    if (!shader->loadFromFile("../Data/shaders/basic.vs", "../Data/shaders/basic.fs")) {
        LOG(Debugability::DBG_LVL_CRITICAL, "Failed to load basic shader.");
    }
    Graphics::g_fontManager = new Graphics::FontManager();
    Graphics::Font* font = Graphics::g_fontManager->getDefaultFont();

    // Input
    Input::g_input = new Input::InputProcessor(hInstance, g_mainWindow->GetHWindow());

    int32 cur_seed = kStartSeed;
    CreateWorld(cur_seed);

    Debugability::Timer timer;
    Debugability::PerformanceMonitor perfMon(&timer);

    mainCam = new BlockyVoyages::World::PerspectiveCamera(45.0f, static_cast<float32>(kWidth) / kHeight, 0.1f, 1000.0f);
    std::vector<Vector4f> viewPlanes = mainCam->getViewFrustumPlanes();
    Input::CameraController controller(mainCam);
    std::unique_ptr<World::Character> player(new World::Character(g_world_map.get()));
    player->attachCamera(&controller);
    player->setPosition(playerStartPos);
    Input::g_input->setController(player->getController());

    // The GUI View doesn't change, so just use a matrix for it and not a camera
    BlockyVoyages::Matrix44f guiCamMat;
    guiCamMat.OrthographicGL(0.0f, static_cast<float32>(kWidth),
                             0.0f, static_cast<float32>(kHeight),
                             -1.0f, 1.0f);

    // Do some initial setup of OpenGL
    glClearColor(0.25f, 0.25f, 0.5f, 1.0f);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    Vector3f white(1.0f, 1.0f, 1.0f);

    shader->bind();
    // Main message loop:
    while (msg.message != WM_QUIT && !logger->criticalMessageSeen()) {
        if (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
        else {
            perfMon.endFrame();
            perfMon.startFrame();
            // use the time difference from the performance monitor since avoid duplicating
            // the calculation
            float32 dt = perfMon.getPrevFrameTime();

            // process input from the player
            controller.setTimeDelta(dt);
            Input::g_input->process();

            bool world_changed = false;
            if (Input::g_input->keyChanged(DIK_RIGHT) && Input::g_input->keyDown(DIK_RIGHT)) {
                cur_seed++;
                world_changed = true;
            } else if (Input::g_input->keyChanged(DIK_LEFT) && Input::g_input->keyDown(DIK_LEFT)) {
                cur_seed--;
                world_changed = true;
            }
            if (world_changed) {
                CreateWorld(cur_seed);
                
                player.reset(new World::Character(g_world_map.get()));
                player->attachCamera(&controller);
                player->setPosition(playerStartPos);
                Input::g_input->setController(player->getController());
            }

            if (Input::g_input->keyChanged(DIK_K) && Input::g_input->keyDown(DIK_K)) {
                if (player->getAttachedCamera()) {
                    player->attachCamera(NULL);
                    Input::g_input->setController(&controller);
                } else {
                    player->attachCamera(&controller);
                    Input::g_input->setController(player->getController());
                }
            }
            if (Input::g_input->keyChanged(DIK_P) && Input::g_input->keyDown(DIK_P)) {
                perfMon.ResetStats();
            }
            if (Input::g_input->keyChanged(DIK_O) && Input::g_input->keyDown(DIK_O)) {
                enableOctree = !enableOctree;
            }

            // run a simulation step
            player->update(dt);

            if(player->getPosition().y < -100.0f) {
                player->setPosition(playerStartPos);
            }

            // draw the world
            glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
            glEnable(GL_DEPTH_TEST);

            if (updateViewable || viewPlanes.size() != 6) {
                viewPlanes = mainCam->getViewFrustumPlanes();
            }
            g_world_map->Draw(mainCam);
            if (enableOctree) {
                g_world_map->DrawCells(mainCam);
            }
            
            player->draw();
            
            // draw the GUI
            
            shader->setUniform("proj_mat", guiCamMat);
            Matrix44f id;
            id.Identity();
            shader->setUniform("model_view_mat", id);

            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glDisable(GL_DEPTH_TEST);

            perfMon.drawStats(1, 200, font);

            std::stringstream sstream;
            sstream << "Camera Pos: " << mainCam->getPosition().x << ", " << mainCam->getPosition().y << ", " << mainCam->getPosition().z;
            font->drawText(1, 300, white, sstream.str());
            
            sstream.clear();
            sstream.str("");
            sstream << "Player Pos: " << player->getPosition().x << ", " << player->getPosition().y << ", " << player->getPosition().z;
            font->drawText(1, 320, white, sstream.str());

            sstream.clear();
            sstream.str("");
            sstream << "Heading: " << player->getHeading() << "; Pitch: " << player->getPitch();
            font->drawText(0, 340, white, sstream.str());

            glBindTexture(GL_TEXTURE_2D, 0);

            Graphics::g_opengl->present();

#if defined (_DEBUG)
            glError = glGetError();
            if (glError != 0) {
                LOG(Debugability::DBG_LVL_HIGH, "Got OpenGL error %d.", glError);
            }
#endif
        }
    }

    delete Input::g_input;

    // shutdown the engine components
    delete mainCam;
    delete Graphics::g_fontManager;
    delete Graphics::g_opengl;
    // No OpenGL commands should be made after this point

    if (fullscreen) {
        ChangeDisplaySettings(nullptr, 0);            // If So Switch Back To The Desktop
        ShowCursor(true);                             // Show Mouse Pointer
    }

    delete g_mainWindow;

    // This should be last thing before returning.
    delete logger;
	return (int)msg.wParam;
}
