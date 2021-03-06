#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::mat3;
using glm::ivec2;

#define RotationSpeed 0.05f	//Camera rotation speed
#define MoveSpeed 0.05f
#define LightMoveSpeed 0.01f

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
//float lightDepthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
int t;
mat3 cameraR;
vec3 cameraPos( 0, 0, -3.001 );
float f = 1.0f;
float yaw = 0.0f;
vector<Triangle> triangles;

/*LIGHT VALUES - POINT LIGHT*/

vec3 lightPos(0.f,-0.5f,-0.7f);
vec3 lightPower = 5.1f*vec3( 1.0f, 1.0f, 1.0f );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

/*DIRECTIONAL LIGHT*/
vec3 lightDirection(1.f, -0.5f, -0.7f);

/* STRUCTS */

struct Pixel
{
    int x;
    int y;
    float zinv;
    vec3 pos3d;
    float z_value;
    vec3 toLight; //Keeps the distance to light
};

struct Vertex
{
	vec3 position;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void VertexShader( const vec3& v, Pixel& p );
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color, vec3 currentNormal, vec3 currentReflactance );
void ComputePolygonRows( const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels );
void DrawRows( const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels, vec3 color, vec3 currentNormal, vec3 currentReflactance );
void DrawPolygon( const vector<vec3>& vertices, vec3 color, vec3 currentNormal, vec3 currentReflactance );
void PixelShader( const Pixel& p, vec3 currentNormal, vec3 currentReflactance );
vec3 ComputePixelReflectedLight( const Pixel& p, vec3 currentNormal, vec3 currentReflactance );
vec3 ComputePixelDirectionalLight( const Pixel& p, vec3 currentNormal, vec3 currentReflactance );

int main( int argc, char* argv[] )
{
	LoadTestModel(triangles);

	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms.\n";

    Uint8* keystate = SDL_GetKeyState( 0 );

    //Get camera direction
	vec3 right(cameraR[0][0], cameraR[0][1], cameraR[0][2]);
	vec3 down(cameraR[1][0], cameraR[1][1], cameraR[1][2]);
	vec3 forward(cameraR[2][0], cameraR[2][1], cameraR[2][2]);

	//Control camera
    if( keystate[SDLK_UP] )
    {
		cameraPos += MoveSpeed*forward; //Forwards
    }
    if( keystate[SDLK_DOWN] )
    {
		cameraPos-= MoveSpeed*forward;  //Backwards
    }
    if( keystate[SDLK_j] )
    {
		cameraPos -= MoveSpeed*right;  //Left
    }
    if( keystate[SDLK_l] ) {
		cameraPos += MoveSpeed*right;  //Right
	}
    if( keystate[SDLK_i] )
    {
		cameraPos -= MoveSpeed*down;  //Down
    }
    if( keystate[SDLK_k] ) {
		cameraPos += MoveSpeed*down;  //Up
	}

	//Rotate on Y axis
    if( keystate[SDLK_LEFT] )
    {
		yaw += RotationSpeed;
		cameraR = mat3(cos(yaw), 0, sin(yaw), 0, 1, 0, -sin(yaw), 0, cos(yaw));
    }

	if( keystate[SDLK_RIGHT] ) {
		yaw -= RotationSpeed;
		cameraR = mat3(cos(yaw), 0, sin(yaw), 0, 1, 0, -sin(yaw), 0, cos(yaw));
	}

		//Control lights WASD
    if( keystate[SDLK_a] )
    {
		lightPos.x -= LightMoveSpeed;	//Left
		lightDirection.x -= LightMoveSpeed;
    }
    if( keystate[SDLK_d] ) {
		lightPos.x += LightMoveSpeed;	//Right
		lightDirection.x += LightMoveSpeed;
	}
    if( keystate[SDLK_w] )
    {
		lightPos.y -= LightMoveSpeed;	//Up
		lightDirection.y -= LightMoveSpeed;
    }
    if( keystate[SDLK_s] ) {
		lightPos.y += LightMoveSpeed;	//Down
		lightDirection.y += LightMoveSpeed;
	}
	if( keystate[SDLK_x] ) {
		lightPos.z += LightMoveSpeed;	//Down
		lightDirection.z += LightMoveSpeed;
	}
	if( keystate[SDLK_z] ) {
		lightPos.z -= LightMoveSpeed;	//Down
		lightDirection.z -= LightMoveSpeed;
	}
}

void VertexShader( const Vertex& v, Pixel& p ) 
{
	vec3 p_p = vec3(v.position.x, v.position.y, v.position.z);	

	p_p = (p_p - cameraPos) * cameraR;

	if (p_p.z == 0)
		return;

	p.zinv = 1.0f/p_p.z;
	p.x = (int)((f*p_p.x/p_p.z)*(SCREEN_WIDTH/2.0f) + SCREEN_WIDTH/2.0f);
	p.y = (int)((f*p_p.y/p_p.z)*(SCREEN_HEIGHT/2.0f) + SCREEN_HEIGHT/2.0f);	

	//Storing the 3D position of the Vertex to the corresponding 
	//variable in Pixel
	p.pos3d.x = v.position.x;
	p.pos3d.y = v.position.y;
	p.pos3d.z = v.position.z;
	
	p.pos3d = p.pos3d/p_p.z;
	p.z_value = p_p.z;

	p.toLight = normalize(p_p - lightPos);
}

vec3 ComputePixelDirectionalLight( const Pixel& p, vec3 currentNormal, vec3 currentReflactance)
{
	vec3 n = currentNormal;
	float diffuseIntensity = max(0.0f, dot(normalize(n), -lightDirection));
	vec3 outputLight = currentReflactance * (indirectLightPowerPerArea + diffuseIntensity);

	return outputLight;
}

vec3 ComputePixelReflectedLight( const Pixel& p, vec3 currentNormal, vec3 currentReflactance )
{
	float radius = glm::distance(lightPos, p.pos3d);
	vec3 r = normalize(lightPos - p.pos3d);
	vec3 n = currentNormal;
	vec3 term1 = lightPower * max(dot(r,n), 0.f);
	float term2 = 4.0f * M_PI * radius * radius;
	vec3 D = term1/term2;
	
	return currentReflactance * (D + indirectLightPowerPerArea);
}

void PixelShader( Pixel& p, vec3 currentColor, vec3 currentNormal, vec3 currentReflactance )
{
	int x = p.x;
	int y = p.y;
	if( p.zinv > depthBuffer[y][x] )
	{
		depthBuffer[y][x] = p.zinv;
		vec3 R = ComputePixelReflectedLight(p, currentNormal, currentReflactance);
		R = ComputePixelDirectionalLight(p, currentNormal, currentReflactance);
		PutPixelSDL( screen, x, y, currentColor * R);
	}
}

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result )
{
	vec2 _a = vec2(a.x, a.y);
	vec2 _b = vec2(b.x, b.y);

	int N = result.size();
	vec2 step = vec2(_b - _a) / float(max(N-1,1));
	vec2 current( _a );

	float depthStep = (b.zinv - a.zinv) / float(max(N-1,1));
	float currentDepth = a.zinv;

	vec3 posStep = (b.pos3d - a.pos3d) / float(max(N-1,1));
	vec3 currentPos = a.pos3d;
	
	for( int i=0; i<N; ++i )
	{
		result[i].x = current.x;
		result[i].y = current.y;
		result[i].zinv = currentDepth;
		result[i].pos3d = currentPos;

		current += step;
		currentDepth += depthStep;
		currentPos += posStep;
	}
}

void DrawLineSDL( SDL_Surface* surface, Pixel a, Pixel b, vec3 color, vec3 currentNormal, vec3 currentReflactance )
{
	vec2 _a = vec2(a.x, a.y);
	vec2 _b = vec2(b.x, b.y);

	ivec2 delta = abs(_a - _b);
	int pixels = max(delta.x, delta.y) + 1;

	vector<Pixel> result(pixels);
	Interpolate(a, b, result);

	for( uint j = 0; j < result.size(); ++j )
	{
	 	PixelShader(result[j], color, currentNormal, currentReflactance);
	}
}

void ComputePolygonRows( const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels )
{
	// 1. Find max and min y-value of the polygon
	//    and compute the number of rows it occupies.
	int maxY = max( vertexPixels[0].y, max(vertexPixels[1].y, vertexPixels[2].y)); 
	int minY = min( vertexPixels[0].y, min(vertexPixels[1].y, vertexPixels[2].y)); 
	//Number of rows occupied by the triangle
	int rows = maxY - minY + 1;
	// 2. Resize leftPixels and rightPixels
	//    so that they have an element for each row.
	leftPixels.resize(rows);
	rightPixels.resize(rows);
	// 3. Initialize the x-coordinates in leftPixels
	//    to some really large value and the x-coordinates
	//    in rightPixels to some really small value.
	for( int i=0; i<rows; ++i) 
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
		// fill in the y values
		leftPixels[i].y = rightPixels[i].y = i + minY;
	}
	// 4. Loop through all edges of the polygon and use
	//    linear interpolation to find the x-coordinate for
	//    each row it occupies. Update the corresponding
	//    values in rightPixels and leftPixels.
	for(int e=0; e<3; ++e) 
	{
		int ne = (e+1)%3; //next edge
		int edgePixels = abs(vertexPixels[e].y - vertexPixels[ne].y) + 1;
		vector<Pixel> result(edgePixels);
		Interpolate(vertexPixels[e], vertexPixels[ne], result);
		for (uint i=0; i<result.size(); ++i) 
		{
			//Obtain the relative index (1 to rows) of the interpolated point of the edge
			//Otherwise it will have the y coordinate relative to the image space
			int index = result[i].y - minY; 
			if(result[i].x < leftPixels[index].x)
			{
				leftPixels[index].x = result[i].x;
				leftPixels[index].zinv = result[i].zinv;
				leftPixels[index].pos3d = result[i].pos3d;
			}
			if(result[i].x > rightPixels[index].x)
			{
				rightPixels[index].x = result[i].x;
				rightPixels[index].zinv = result[i].zinv;
				rightPixels[index].pos3d = result[i].pos3d;
			}  
		}
	}
}

void DrawRows( const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 color, vec3 currentNormal, vec3 currentReflactance ) 
{
	for(uint i=0; i<leftPixels.size(); i++) 
	{
		if( (leftPixels[i].y < SCREEN_HEIGHT && leftPixels[i].y > 0) || (rightPixels[i].y < SCREEN_HEIGHT && rightPixels[i].y > 0)) 
		{
			DrawLineSDL(screen,leftPixels[i],rightPixels[i],color, currentNormal, currentReflactance);
		}
	}
}

void DrawPolygon( const vector<Vertex>& vertices, vec3 color, vec3 currentNormal, vec3 currentReflactance )
{
    int V = vertices.size();
    vector<Pixel> vertexPixels( V );
    for( int i=0; i<V; ++i )
		VertexShader( vertices[i], vertexPixels[i] );
    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
	DrawRows( leftPixels, rightPixels, color, currentNormal, currentReflactance);
}


void Draw() 
{
	SDL_FillRect( screen, 0, 0 );
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	for( uint i=0; i<SCREEN_HEIGHT; ++i )
	{
		for( uint j=0; j<SCREEN_WIDTH; ++j )
		{
			depthBuffer[i][j] = 0.0f;
		}
	}

	// #pragma omp parallel for
	for( uint i=0; i<triangles.size(); ++i )
	{
		vector<Vertex> vertices(3);

		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;

		vec3 currentNormal = triangles[i].normal;
		vec3 currentReflactance = triangles[i].color;

		DrawPolygon(vertices, triangles[i].color, currentNormal, currentReflactance);
	}

    if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);
    SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
