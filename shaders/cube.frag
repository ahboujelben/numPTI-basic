/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
////////////////////////////////////////////////////////////////////////////

#version 330 core
out vec4 FragColor;

in FragData
{
    flat vec4 cubeCornerPos;
    flat vec3 Normal;
    flat vec3 cubeColor;
};

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;
uniform mat4 view;
uniform mat3 normalMatrix;
uniform vec3 phase1Color;
uniform vec3 phase2Color;
uniform vec3 phase3Color;
uniform vec3 phase4Color;

void main()
{
    vec3 objectColor;
    if(cubeColor.x==1)
        objectColor=phase1Color+cubeColor.y*(phase4Color-phase1Color);
    if(cubeColor.x==2)
        objectColor=phase2Color+cubeColor.y*(phase4Color-phase2Color);
    if(cubeColor.x==3)
        objectColor=phase3Color+cubeColor.y*(phase4Color-phase3Color);


    vec3 FragPos=vec3(view*cubeCornerPos);
    // ambient
    float ambientStrength = 0.3;
    vec3 ambient = ambientStrength * lightColor;

    //vec3 lightPosView= vec3(view*vec4(lightPos,1.0));
    // diffuse
    vec3 norm = normalize(normalMatrix*Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    // specular
    float specularStrength = 0.0;
    vec3 viewDir = normalize(- FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = specularStrength * spec * lightColor;

    vec3 result = (ambient + diffuse + specular) * objectColor;
    FragColor = vec4(result, 1.0);
}
