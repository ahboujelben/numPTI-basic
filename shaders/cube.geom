/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
////////////////////////////////////////////////////////////////////////////

#version 330
#extension GL_EXT_gpu_shader4 : enable

layout(points) in;
layout(triangle_strip, max_vertices=36) out;

uniform mat4 projection;
uniform mat4 view;

in VertexData
{
    vec3 cubePos;
    vec3 cubeSize;
    vec3 cubeColor;
} vert[];

out FragData
{
    flat vec4 cubeCornerPos;
    flat vec3 Normal;
    flat vec3 cubeColor;
};

void main()
{

    cubeColor=vert[0].cubeColor;
    vec3 size=vert[0].cubeSize;

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(-0.5, -0.5, -0.5),1.0);
    Normal= vec3(0.0,  0.0, -1.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(-0.5, 0.5, -0.5),1.0);
    Normal= vec3(0.0,  0.0, -1.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(0.5, -0.5, -0.5),1.0);
    Normal= vec3(0.0,  0.0, -1.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(0.5, 0.5, -0.5),1.0);
    Normal= vec3(0.0,  0.0, -1.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

/////////////

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(0.5, 0.5, 0.5),1.0);
    Normal= vec3(1.0,  0.0, 0.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(-0.5, 0.5, -0.5),1.0);
    Normal= vec3(0.0,  1.0, 0.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(-0.5, 0.5, 0.5),1.0);
    Normal= vec3(0.0,  1.0, 0.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(-0.5, -0.5, 0.5),1.0);
    Normal= vec3(-1.0,  0.0, 0.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    //////////

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(0.5, 0.5, 0.5),1.0);
    Normal= vec3(0.0,  0.0, 1.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(0.5, -0.5, 0.5),1.0);
    Normal= vec3(0.0,  0.0, 1.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(0.5, -0.5, -0.5),1.0);
    Normal= vec3(1.0,  0.0, 0.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(-0.5, -0.5, 0.5),1.0);
    Normal= vec3(0.0,  -1.0, 0.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

////////

cubeCornerPos = vec4(vert[0].cubePos+size*vec3(-0.5, -0.5, -0.5),1.0);
    Normal= vec3(0.0,  -1.0, 0.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();

    cubeCornerPos = vec4(vert[0].cubePos+size*vec3(-0.5, 0.5, -0.5),1.0);
    Normal= vec3(-1.0,  0.0, 0.0);
    gl_Position = projection * view * cubeCornerPos;
    gl_PrimitiveID = gl_PrimitiveIDIn;
    EmitVertex();
}
