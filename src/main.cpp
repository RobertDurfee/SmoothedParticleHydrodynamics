#include "gl.h"
#include <GLFW/glfw3.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <assert.h>

#include "vertexrecorder.h"
#include "util.h"
#include "camera.h"
#include "glprogram.h"
#include "sph.h"

using namespace std;

namespace {

    void initSystem();
    void stepSystem();
    void drawSystem();
    void freeSystem();
    void resetTime();
    void initRendering();
    void drawAxis();

    const Vector3f LIGHT_POS(3.0f, 3.0f, 5.0f);
    const Vector3f LIGHT_COLOR(120.0f, 120.0f, 120.0f);
    const Vector3f FLOOR_COLOR(1.0f, 0.0f, 0.0f);

    uint64_t start_tick;
    double elapsed_s;
    double simulated_s;
    Camera camera;
    bool gMousePressed = false;
    GLuint program_color;
    GLuint program_light;
    SPH *sph;

    static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (action == GLFW_RELEASE) { // only handle PRESS and REPEAT
            return;
        }
        // Special keys (arrows, CTRL, ...) are documented
        // here: http://www.glfw.org/docs/latest/group__keys.html
        switch (key) {
            case GLFW_KEY_ESCAPE: { // Escape key
                exit(0);
                break;
            } case ' ': {
                Matrix4f eye = Matrix4f::identity();
                camera.SetRotation(eye);
                camera.SetCenter(Vector3f(0, 0, 0));
                break;
            } case 'R': {
                cout << "Resetting simulation\n";
                freeSystem();
                initSystem();
                resetTime();
                break;
            } default: {
                cout << "Unhandled key press " << key << "." << endl;
            }
        }
    }

    static void mouseCallback(GLFWwindow* window, int button, int action, int mods) {
        double xd, yd;
        glfwGetCursorPos(window, &xd, &yd);
        int x = (int)xd;
        int y = (int)yd;
        int lstate = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
        int rstate = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
        int mstate = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE);
        if (lstate == GLFW_PRESS) {
            gMousePressed = true;
            camera.MouseClick(Camera::LEFT, x, y);
        }
        else if (rstate == GLFW_PRESS) {
            gMousePressed = true;
            camera.MouseClick(Camera::RIGHT, x, y);
        }
        else if (mstate == GLFW_PRESS) {
            gMousePressed = true;
            camera.MouseClick(Camera::MIDDLE, x, y);
        }
        else {
            gMousePressed = true;
            camera.MouseRelease(x, y);
            gMousePressed = false;
        }
    }

    static void motionCallback(GLFWwindow* window, double x, double y) {
        if (!gMousePressed) {
            return;
        }
        camera.MouseDrag((int)x, (int)y);
    }

    void setViewport(GLFWwindow* window) {
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        camera.SetDimensions(w, h);
        camera.SetViewport(0, 0, w, h);
        camera.ApplyViewport();
    }

    void drawAxis() {
        glUseProgram(program_color);
        Matrix4f M = Matrix4f::translation(camera.GetCenter()).inverse();
        camera.SetUniforms(program_color, M);
        const Vector3f DKRED(1.0f, 0.5f, 0.5f);
        const Vector3f DKGREEN(0.5f, 1.0f, 0.5f);
        const Vector3f DKBLUE(0.5f, 0.5f, 1.0f);
        const Vector3f GREY(0.5f, 0.5f, 0.5f);
        const Vector3f ORGN(0, 0, 0);
        const Vector3f AXISX(5, 0, 0);
        const Vector3f AXISY(0, 5, 0);
        const Vector3f AXISZ(0, 0, 5);
        VertexRecorder recorder;
        recorder.record_poscolor(ORGN, DKRED);
        recorder.record_poscolor(AXISX, DKRED);
        recorder.record_poscolor(ORGN, DKGREEN);
        recorder.record_poscolor(AXISY, DKGREEN);
        recorder.record_poscolor(ORGN, DKBLUE);
        recorder.record_poscolor(AXISZ, DKBLUE);
        recorder.record_poscolor(ORGN, GREY);
        recorder.record_poscolor(-AXISX, GREY);
        recorder.record_poscolor(ORGN, GREY);
        recorder.record_poscolor(-AXISY, GREY);
        recorder.record_poscolor(ORGN, GREY);
        recorder.record_poscolor(-AXISZ, GREY);
        glLineWidth(3);
        recorder.draw(GL_LINES);
    }

    void initSystem() {
        sph = new SPH();
    }

    void freeSystem() {
        delete sph;
        sph = nullptr;
    }

    void resetTime() {
        elapsed_s = 0;
        simulated_s = 0;
        start_tick = glfwGetTimerValue();
    }

    void stepSystem() {
        assert(elapsed_s - simulated_s >= 0);
        sph->update(elapsed_s - simulated_s);
        simulated_s = elapsed_s;
    }

    void drawSystem() {
        GLProgram glProgram(program_light, program_color, &camera);
        glProgram.updateLight(LIGHT_POS, LIGHT_COLOR.xyz());
        sph->draw(glProgram);
    }

    void initRendering() {
        glClearColor(0, 0, 0, 1);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
}

int main(int argc, char** argv) {
    GLFWwindow* window = createOpenGLWindow(1024, 1024, "Final Project");
    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseCallback);
    glfwSetCursorPosCallback(window, motionCallback);
    initRendering();
    program_color = compileProgram(c_vertexshader, c_fragmentshader_color);
    if (!program_color) {
        printf("Cannot compile program\n");
        return -1;
    }
    program_light = compileProgram(c_vertexshader, c_fragmentshader_light);
    if (!program_light) {
        printf("Cannot compile program\n");
        return -1;
    }
    camera.SetDimensions(600, 600);
    camera.SetPerspective(50);
    camera.SetDistance(10);
    initSystem();
    uint64_t freq = glfwGetTimerFrequency();
    resetTime();
    double now = (double)glfwGetTimerValue();
    double step = (double)(now - start_tick) / freq;
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        setViewport(window);
        if (gMousePressed) {
            drawAxis();
        }
        now = (double)(glfwGetTimerValue() - start_tick) / freq;
        std::cout << "Time between frames: " << now - step << "s" << std::endl;
        step = now;
        elapsed_s += 1.0 / 30.0; // Render 30 fps no matter how long it takes
        stepSystem();
        drawSystem();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glDeleteProgram(program_color);
    glDeleteProgram(program_light);
    return 0;
}
