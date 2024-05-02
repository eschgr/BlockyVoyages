#version 330

out vec4 frag_color;
in float moisture_;

void main() {
    vec3 yellow = vec3(1.0f, 1.0f, 0.0f);
    vec3 green = vec3(0.0f, 1.0f, 0.0f);
    vec3 cyan = vec3(0.0f, 1.0f, 1.0f);
    vec3 red = vec3(1.0, 0.0f, 0.0f);
    vec3 orange = vec3(1.0, 0.5, 0.0);
    
    vec3 color0;
    vec3 color1;
    float moisture = moisture_;
    if (moisture_ >= 7.5f) {
        moisture = moisture - 7.5f;
        color0 = cyan;
        color1 = green;
    } else if (moisture >= 5.0) {
        moisture -= 5.0;
        color0 = green;
        color1 = yellow;
    } else if (moisture >= 2.5) {
        moisture -= 2.5;
        color0 = yellow;
        color1 = orange;
    } else {
        color0 = orange;
        color1 = red;
    }
    moisture = clamp(moisture / 2.5f, 0.0f, 1.0f);
    vec3 color = moisture * color0 + (1.0 - moisture) * color1;
    frag_color = vec4(color, 1.0f);
}