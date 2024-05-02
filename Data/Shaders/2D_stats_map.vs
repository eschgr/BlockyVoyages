#version 330

uniform mat4 model_mat;
uniform mat4 view_proj_mat;
uniform float max_height;

in vec4 position;
in float height;
in float temperature;
in float moisture;

out float moisture_;

void main() {
	gl_Position = view_proj_mat * model_mat * vec4(position.x, height, -position.y, 1.0f);
    moisture_ = moisture;
}