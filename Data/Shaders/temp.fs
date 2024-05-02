#version 330

out vec4 frag_color;
in float height_;

void main() {
    vec3 color;
    vec3 color0 = vec3(0.0f, 1.0f, 1.0f);
    vec3 color1 = vec3(0.5f, 1.0f, 1.0f);
    vec3 color2 = vec3(1.0f, 1.0f, 1.0f);
    vec3 color3 = vec3(0.0f, 1.0f, 0.5f);
    vec3 color4 = vec3(0.0f, 1.0f, 0.0f);
    vec3 color5 = vec3(0.5f, 1.0f, 0.0f);
    vec3 color6 = vec3(1.0f, 1.0f, 0.0f);
    vec3 color7 = vec3(0.5f, 0.5f, 0.5f);
    vec3 color8 = vec3(1.0f, 0.5f, 0.0f);
    vec3 color9 = vec3(1.0f, 0.0f, 0.0f);
    if (height_ < -20.0f) {
        color = color0;
    } else if (height_ < -10.0f) {
        color = color1;
    } else if (height_ < 0.0f) {
        color = color2;
    } else if (height_ < 10.0f) {
        color = color3;
    } else if (height_ < 20.0f) {
        color = color4;
    } else if (height_ < 30.0f) {
        color = color5;
    } else if (height_ < 40.0f) {
        color = color6;
    } else if (height_ < 50.0f) {
        color = color7;
    } else if (height_ < 60.0f) {
        color = color8;
    } else {
        color = color9;
    }
    frag_color = vec4(color, 1.0f);
}