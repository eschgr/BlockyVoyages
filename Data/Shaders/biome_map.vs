#version 330

uniform mat4 model_mat;
uniform mat4 view_proj_mat;

in vec4 position;
in float height;

void main() {
	gl_Position = view_proj_mat * model_mat * vec4(position.x, height, -position.y, 1.0f);
}