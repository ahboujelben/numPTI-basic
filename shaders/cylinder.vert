/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
///              Parts of this code are based on PyMOL Open source project
/// Licence:     Attribution-NonCommercial 4.0 International
////////////////////////////////////////////////////////////////////////////

#version 330 core
#extension GL_EXT_gpu_shader4 : enable

layout(location = 0) in vec3  CylinderPosition;
layout(location = 1) in float CylinderExt;
layout(location = 2) in vec3  CylinderDirection;
layout(location = 3) in vec4  CylinderColor;
layout(location = 4) in float CylinderRadius;

out vec3 cylinder_color_in;
out vec3 cylinder_direction_in;
out float cylinder_radius_in;
out float cylinder_ext_in;


void main()
{  
  cylinder_color_in = CylinderColor.xyz;
  cylinder_direction_in = normalize(CylinderDirection);
  cylinder_radius_in = CylinderRadius;
  cylinder_ext_in = CylinderExt;
  gl_Position = vec4(CylinderPosition,1.0);
}
