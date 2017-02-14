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
float f = 1.0f;
float yaw = 0.0f;
vector<Triangle> triangles;

void Update();
void Draw();
void VertexShader( const vec3& v, ivec2& p );
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color );
void ComputePolygonRows( const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels );
void DrawRows( const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels, vec3 color );
void DrawPolygon( const vector<vec3>& vertices, vec3 color );

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
}

void VertexShader( const vec3& v, ivec2& p ) {
	vec3 p_p = vec3(v[0], v[1], v[2]);	

	p_p = (p_p - cameraPos) * cameraR;

	if (p_p.z == 0)
		return;

	p.x = (int)((f*p_p.x/p_p.z)*(SCREEN_WIDTH/2.0f) + SCREEN_WIDTH/2.0f);
	p.y = (int)((f*p_p.y/p_p.z)*(SCREEN_HEIGHT/2.0f) + SCREEN_HEIGHT/2.0f);	

}

void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result ){
	int N = result.size();
	vec2 step = vec2(b-a) / float(max(N-1,1));
	vec2 current( a );
	for( int i=0; i<N; ++i )
	{
		result[i] = current;
		current += step;
	}
}

void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color ){
	ivec2 delta = abs(a - b);
	int pixels = max(delta.x, delta.y) + 1;
	vector<ivec2> result(pixels);
	Interpolate(a, b, result);
	for( uint j = 0; j < result.size(); ++j){
		PutPixelSDL( screen, result[j].x, result[j].y, color );
	}
}

void ComputePolygonRows( const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels ){
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
		vector<ivec2> result(edgePixels);
		Interpolate(vertexPixels[e], vertexPixels[ne], result);
		for (uint i=0; i<result.size(); ++i) 
		{
			//Obtain the relative index (1 to rows) of the interpolated point of the edge
			//Otherwise it will have the y coordinate relative to the image space
			int index = result[i].y - minY; 
			if(result[i].x < leftPixels[index].x)
			{
				leftPixels[index].x = result[i].x;
			}
			if(result[i].x > rightPixels[index].x)
			{
				rightPixels[index].x = result[i].x;
			}  
		}
	}
}

void DrawRows( const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels, vec3 color ) {
	for(uint i=0; i<leftPixels.size(); i++) {
		if( (leftPixels[i].y < SCREEN_HEIGHT && leftPixels[i].y > 0) || (rightPixels[i].y < SCREEN_HEIGHT && rightPixels[i].y > 0)) {
			DrawLineSDL(screen, leftPixels[i],rightPixels[i],color);
		}
	}
}

void DrawPolygon( const vector<vec3>& vertices, vec3 color )
{
    int V = vertices.size();
    vector<ivec2> vertexPixels( V );
    for( int i=0; i<V; ++i )
		VertexShader( vertices[i], vertexPixels[i] );
    vector<ivec2> leftPixels;
    vector<ivec2> rightPixels;
    ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
	DrawRows( leftPixels, rightPixels, color);
}


void Draw() {
	SDL_FillRect( screen, 0, 0 );
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	//DEBUGGING
	// vector<ivec2> vertexPixels(3);
	// vertexPixels[0] = ivec2(10, 5);
	// vertexPixels[1] = ivec2( 5,10);
	// vertexPixels[2] = ivec2(15,15);
	// vector<ivec2> leftPixels;
	// vector<ivec2> rightPixels;
	// ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
	// for( uint row=0; row<leftPixels.size(); ++row )
	// {
	//     cout << "Start: ("
	//  << leftPixels[row].x << ","
	//  << leftPixels[row].y << "). "
	//  << "End: ("
	//  << rightPixels[row].x << ","
	//  << rightPixels[row].y << "). " << endl;
	// }
	//!DEBUGGING

	// #pragma omp parallel for
	for( uint i=0; i<triangles.size(); ++i )
	{
		vector<vec3> vertices(3);
		vector<ivec2> vertices2D(3);

		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		//WIREFRAME (OLD)
		// vec3 color(1,1,1);
		// for(int v=0; v<3; ++v)
		// {
		//     ivec2 projPos;
		//     VertexShader( vertices[v], projPos );
		//     vertices2D[v] = projPos;
		//     PutPixelSDL( screen, projPos.x, projPos.y, color );
		// }

		// DrawLineSDL( screen, vertices2D[0], vertices2D[1], color );
		// DrawLineSDL( screen, vertices2D[0], vertices2D[2], color );
		// DrawLineSDL( screen, vertices2D[1], vertices2D[2], color );
		//!WIREFRAME (OLD)

		DrawPolygon(vertices, triangles[i].color);

	}
    if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);
    SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
