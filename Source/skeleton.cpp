#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::ivec2;

#define RotationSpeed 0.05f	//Camera rotation speed
#define MoveSpeed 0.2f

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */


mat3 cameraR;
//vec3 camera(0.0f,0.0f, -4.0f);
vec3 cameraPos( 0, 0, -3.001 );
float f = 500.0f;
float yaw = 0.0f;
vector<Triangle> triangles;

void Update();
void Draw();
void VertexShader( const vec3& v, ivec2& p );


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
}

void VertexShader( const vec3& v, ivec2& p ) {	
	float X = v[0] - cameraPos.x;
	float Y = v[1] - cameraPos.y;
	float Z = v[2] - cameraPos.z;
	vec3 p_p = vec3(X,Y,Z);

	p_p = p_p * cameraR;

	if (p_p.z == 0)
		return;

	p.x = int(f*p_p.x/p_p.z) + SCREEN_WIDTH/2.0f;
	p.y = int(f*p_p.y/p_p.z) + SCREEN_HEIGHT/2.0f;	
}

void Draw() {
	SDL_FillRect( screen, 0, 0 );
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);
	for( int i=0; i<triangles.size(); ++i )
	{
		vector<vec3> vertices(3);
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;
		for(int v=0; v<3; ++v)
		{
		    ivec2 projPos;
		    VertexShader( vertices[v], projPos );
		    vec3 color(1,1,1);
		    PutPixelSDL( screen, projPos.x, projPos.y, color );
		}
	}
    if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);
    SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
