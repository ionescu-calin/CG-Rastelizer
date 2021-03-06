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
using glm::vec4;
using glm::mat4;

#define RotationSpeed 0.05f	//Camera rotation speed
#define MoveSpeed 0.01f
#define LightMoveSpeed 0.01f

// #define POSTPROCESSING
// #define CHROMATICABERRATION
// #define GRAINEFFECT
// #define CAMERAWARPEFFECT
int cameraSpeed = 0;
/* ----------------------------------------------------------------------------*/
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

struct Edge 
{
	vec3 v1;
	vec3 v2;
	vec3 normal;
};

struct Adjacencies
{
	vector<Triangle> triangles;
	vector<vec3> v1;
	vector<vec3> v2;
};

struct Quad
{
	vec3 v1;
	vec3 v2;
	vec3 v3;
	vec3 v4;
	vec3 normal;
};

/*CLASSES*/
//Used to define an object of triangular surfaces
class Object
{
public:
	std::vector<Triangle> triangles;
	std::vector<Edge> silhouette;
	std::vector<Quad> shadowQuads;

	Object( std::vector<Triangle> triangles )
		:triangles(triangles)
	{

	}
};

/* GLOBAL VARIABLES */
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;

/*BUFFERS*/
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
int stencilBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 frameBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
//int frameBuffering = 1;
int depthBuffering = 1;
int stencilBuffering = 0;
int renderFrontFaces = 0;
int renderBackFaces = 0;

int t;
mat3 cameraR;
vec3 cameraPos( 0, 0, -3.001 );
float f = 1.0f;
float yaw = 0.0f;
vector<Triangle> triangles;
vector<Object> sceneObjects;
vector<Triangle> redTriangles;
vector<Triangle> blueTriangles;
vec3 red( 0.75f, 0.15f, 0.15f );
vec3 blue( 0.15f, 0.15f, 0.75f );
vec3 pink( 1.000, 0.078, 0.576 );

/*LIGHT VALUES - POINT LIGHT*/
vec3 lightPos(0.21, 0.14, 0.17);//(0.01, -0.19, -0.29);//(0.f,-0.5f,-0.7f);//(0.3f, 0.1f, -0.2f);
vec3 lightPower = 1.1f*vec3( 1.0f, 1.0f, 1.0f );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

/*DIRECTIONAL LIGHT*/
vec3 lightDirection(1.f, -0.5f, -0.7f);

/*SHADOW VOLUME*/
vec3 ExtrdudeMagnitude(1.f, 1.f, 1.f);


//Clipping Bounds:
float minX = -1.0f;
float minY = -1.0f;
float minZ = -1.0f;
float maxX = 1.0f;
float maxY = 1.0f;
float maxZ = 1.0f;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void VertexShader( const Vertex& v, Pixel& p );
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color, vec3 currentNormal, vec3 currentReflactance );
void ComputePolygonRows( const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels );
void DrawRows( const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels, vec3 color, vec3 currentNormal, vec3 currentReflactance );
void DrawPolygon( const vector<vec3>& vertices, vec3 color, vec3 currentNormal, vec3 currentReflactance );
void PixelShader( const Pixel& p, vec3 currentNormal, vec3 currentReflactance );
vec3 ComputePixelReflectedLight( const Pixel& p, vec3 currentNormal, vec3 currentReflactance );
vec3 ComputePixelDirectionalLight( const Pixel& p, vec3 currentNormal, vec3 currentReflactance );
vec3 ComputeAmbientLight( vec3 currentReflactance );
void ComputeSilhouettes( vector<Object>& objects );
void ComputeShadowQuads( vector<Object> objects );
void RenderTriangles(vector<Triangle> triangles);
void RederQuad(Quad quad);
void RenderQuads(vector<Quad> quads);
void RenderCaps(vector<Triangle> caps);
void RenderCap(Triangle cap);

int main( int argc, char* argv[] )
{
	LoadTestModel(triangles);

	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer. 

	for( uint i=0; i < triangles.size(); ++i)
	{
		if(triangles[i].color == red)
		{
			redTriangles.push_back(Triangle( triangles[i].v0,  triangles[i].v1,  triangles[i].v2,  triangles[i].color));
		}
		if(triangles[i].color == blue)
		{
			blueTriangles.push_back(Triangle( triangles[i].v0,  triangles[i].v1,  triangles[i].v2,  triangles[i].color));
		} 
	}
	
	Object redCube(redTriangles);
	Object blueCube(blueTriangles);
	sceneObjects.push_back(redCube);
	sceneObjects.push_back(blueCube);

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}


#ifdef POSTPROCESSING

vec3 image[SCREEN_HEIGHT][SCREEN_WIDTH];
void ApplyPostprocessing() {	
	for( int y=0; y<SCREEN_HEIGHT; ++y )
	{
		for( int x=0; x<SCREEN_WIDTH; ++x )
		{
			vec3 c(0.0f, 0.0f, 0.0f);
			//Depth of field using a low pass kernel and the depth buffer
			if(depthBuffer[x][y] > 0.3f && depthBuffer[x][y] < 0.6f) {
					if(x > 1 && y > 1 && y < SCREEN_HEIGHT-2 && x < SCREEN_WIDTH-2) {
					//Hard blur
					c = image[x-2][y-2]/25.0f + image[x-2][y-1]/25.0f + image[x-2][y]/25.0f + image[x-2][y+1]/25.0f + image[x-2][y+2]/25.0f
						+ image[x-1][y-2]/25.0f + image[x-1][y-1]/25.0f + image[x-1][y]/25.0f + image[x-1][y+1]/25.0f + image[x-1][y+2]/25.0f
						+ image[x][y-2]/25.0f + image[x][y-1]/25.0f + image[x][y]/25.0f + image[x][y+1]/25.0f + image[x][y+2]/25.0f
						+ image[x+1][y-2]/25.0f + image[x+1][y-1]/25.0f + image[x+1][y]/25.0f + image[x+1][y+1]/25.0f + image[x+1][y+2]/25.0f
						+ image[x+2][y-2]/25.0f + image[x+2][y-1]/25.0f + image[x+2][y]/25.0f + image[x+2][y+1]/25.0f + image[x+2][y+2]/25.0f;
				} else if (x > 0 && y > 0 && y < SCREEN_HEIGHT-1 && x < SCREEN_WIDTH-1) {
					//Soft blur
					c = image[x][y-1]/9.0f + image[x][y+1]/9.0f + image[x-1][y]/9.0f + image[x+1][y]/9.0f + image[x][y]/9.0f
						+ image[x-1][y-1]/9.0f + image[x-1][y+1]/9.0f + image[x+1][y+1]/9.0f + image[x+1][y-1]/9.0f;
				} 
			} else c = image[x][y];
			
			#ifdef CHROMATICABERRATION
				int r = 3;//rand()%5;
				if(x-r>0 && y-r>0)
					c = (c + image[x-r][y-r].x)/vec3(2.0f, 1.0f, 1.0f);
				if(x+r<SCREEN_WIDTH && y+r < SCREEN_HEIGHT)
					c = (c + image[x+r][y+r].y)/vec3(1.0f, 2.0f, 1.0f);
				if(y-r>0 && x+r<SCREEN_WIDTH)
					c = (c + image[x+r][y-r].z)/vec3(1.0f, 1.0f, 2.0f);
			#endif

			#ifdef GRAINEFFECT
				int g = rand()%2;
				if(g) {
					vec3 grain((float(rand()%100)/100.f), float((rand()%100)/100.f), float((rand()%100)/100.f)); 
					c = (c+grain)/2.0f;
				}
			#endif 

			#ifdef CAMERAWARPEFFECT 
				int cx = SCREEN_WIDTH/2;
				int cy = SCREEN_HEIGHT/2;
				int rr = 200 - cameraSpeed, sr = rr*rr;
				int sx = (x-cx)*(x-cx), sy = (y-cy)*(y-cy);
				int power = sx + sy - sr;
				float maxpower = cx*cx+cy*cy - sr;
				if(power >= 0) {
					float offset = 1.0f - cameraSpeed*(float)power/maxpower;
					c = c*offset;//vec3(0.6f, 0.6f, 0.6f);
				}

			#endif
			PutPixelSDL( screen, x, y, c);
		}
	}
	for( int y=0; y<SCREEN_HEIGHT; ++y )
		for( int x=0; x<SCREEN_WIDTH; ++x ) 
			image[x][y] = vec3(0.0f, 0.0f, 0.0f);
}


#endif


/* 
 * 1) Start with the line clipping algortihm
 * 2) Apply line clipping for all of the edges of the triangle to obtain an array of points
 * 3) reconstruct triangles from that array 
 * 4) Go back and do the normall stuff (go from the homogenous coordinates to the projected)
 TODO: Define xmin, ymin and xmax, ymax
*/

//For the 2D case
const int xmin = -1.f, ymin = -1.f;
const int xmax = 1.f, ymax = 1.f;
typedef int OutCode;

const int INSIDE = 0; // 0000
const int LEFT = 1;   // 0001
const int RIGHT = 2;  // 0010
const int BOTTOM = 4; // 0100
const int TOP = 8;    // 1000

OutCode ComputeOutCode(float x, float y)
{
	OutCode code;

	code = INSIDE;          // initialised as being inside of [[clip window]]

	if (x < xmin)           // to the left of clip window
		code |= LEFT;
	else if (x > xmax)      // to the right of clip window
		code |= RIGHT;
	if (y < ymin)           // below the clip window
		code |= BOTTOM;
	else if (y > ymax)      // above the clip window
		code |= TOP;

	return code;
}

// // TODO: Make it return a list of new vertices
//TODO: Write function that takes the new points and creates triangles:
void CohenSutherlandLineClipAndDraw(float x0, float y0, float x1, float y1, vec4& p0, vec4& p1)
{
	// compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
	OutCode outcode0 = ComputeOutCode(x0, y0);
	OutCode outcode1 = ComputeOutCode(x1, y1);
	bool accept = false;

	while (true) {
		if (!(outcode0 | outcode1)) { // Bitwise OR is 0. Trivially accept and get out of loop
			accept = true;
			break;
		} else if (outcode0 & outcode1) { // Bitwise AND is not 0. (implies both end points are in the same region outside the window). Reject and get out of loop
			break;
		} else {
			// failed both tests, so calculate the line segment to clip
			// from an outside point to an intersection with clip edge
			float x, y;

			// At least one endpoint is outside the clip rectangle; pick it.
			OutCode outcodeOut = outcode0 ? outcode0 : outcode1;

			// Now find the intersection point;
			// use formulas y = y0 + slope * (x - x0), x = x0 + (1 / slope) * (y - y0)
			if (outcodeOut & TOP) {           // point is above the clip rectangle
				x = x0 + (x1 - x0) * (ymax - y0) / (y1 - y0);
				y = ymax;
			} else if (outcodeOut & BOTTOM) { // point is below the clip rectangle
				x = x0 + (x1 - x0) * (ymin - y0) / (y1 - y0);
				y = ymin;
			} else if (outcodeOut & RIGHT) {  // point is to the right of clip rectangle
				y = y0 + (y1 - y0) * (xmax - x0) / (x1 - x0);
				x = xmax;
			} else if (outcodeOut & LEFT) {   // point is to the left of clip rectangle
				y = y0 + (y1 - y0) * (xmin - x0) / (x1 - x0);
				x = xmin;
			}

			// Now we move outside point to intersection point to clip
			// and get ready for next pass.
			if (outcodeOut == outcode0) {
				x0 = x;
				y0 = y;
				outcode0 = ComputeOutCode(x0, y0);
			} else {
				x1 = x;
				y1 = y;
				outcode1 = ComputeOutCode(x1, y1);
			}
		}
	}
	if (accept) {
		// Following functions are left for implementation by user based on
		// their platform (OpenGL/graphics.h etc.)
		// DrawRectangle(xmin, ymin, xmax, ymax);
		// LineSegment(x0, y0, x1, y1);
		// cout<<x0<<" "<<y0<<" "<<x1<<" "<<y1<<endl;
		p0.x = x0; p1.x = x1;
		p0.y = y0; p1.y = y1;

	}


}

bool insideCube(glm::vec4 v)
{
	if(-1.f <= v.x/v.w && v.x/v.w <= 1.f && -1.f <= v.y/v.w && v.y/v.w <= 1.f && -1.f <= v.z/v.w && v.z/v.w <= 1.f) 
		return true;
	return false;
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

	//SetCullingAndClipping();
	//Control camera
    if( keystate[SDLK_UP] )
    {
    	cameraSpeed +=1;
		cameraPos += MoveSpeed*forward; //Forwards
    } 
    if( keystate[SDLK_DOWN] )
    {
    	cameraSpeed +=1;
		cameraPos-= MoveSpeed*forward;  //Backwards
    }
    if(!keystate[SDLK_UP] && !keystate[SDLK_DOWN]) {
    	cameraSpeed = 0;
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

void CheckCountourList( Edge edge, vector<Edge>& contourList )
{
	if( contourList.size() == 0 )
	{
		contourList.push_back(edge);
		return;
	}

	int found = 0;
	for( uint i=0; i<contourList.size(); ++i)
	{
		if( (contourList[i].v1 == edge.v1 && contourList[i].v2 == edge.v2) ||
	   	    (contourList[i].v1 == edge.v2 && contourList[i].v2 == edge.v1) )
		{
			found = 1;
			contourList.erase(contourList.begin() + i);
		}
	}
	
	if(!found) {
		contourList.push_back(edge);
	}
}

void ComputeSilhouettes( vector<Object>& objects, vector<Quad>& quads, vector<Triangle>& caps)
{
	vector<Edge> contourList;
	vector<Edge> triangleEdges;

	for( uint i=0; i<objects.size(); ++i )
	{
		for( uint j=0; j<objects[i].triangles.size(); ++j )
		{
			float x = (objects[i].triangles[j].v0.x + objects[i].triangles[j].v1.x + objects[i].triangles[j].v2.x)/3.f;
			float y = (objects[i].triangles[j].v0.y + objects[i].triangles[j].v1.y + objects[i].triangles[j].v2.y)/3.f;
			float z = (objects[i].triangles[j].v0.z + objects[i].triangles[j].v1.z + objects[i].triangles[j].v2.z)/3.f;

			vec3 averageTrianglePosition = vec3(x,y,z);//(objects[i].triangles[j].v0 + objects[i].triangles[j].v1 + objects[i].triangles[j].v2)/vec3(3.f,3.f,3.f); //don't know if this is correct
			vec3 incidentLightDir = averageTrianglePosition - lightPos;
			if( dot(incidentLightDir, objects[i].triangles[j].normal) >= 0.0f )
			{
				Edge edge1, edge2, edge3;
				//edge 1
				edge1.v1 = objects[i].triangles[j].v0;
				edge1.v2 = objects[i].triangles[j].v1;
				edge1.normal = objects[i].triangles[j].normal;
				triangleEdges.push_back(edge1);
				//edge 2
				edge2.v1 = objects[i].triangles[j].v0;
				edge2.v2 = objects[i].triangles[j].v2;
				edge2.normal = objects[i].triangles[j].normal;
				triangleEdges.push_back(edge2);
				//edge 3
				edge3.v1 = objects[i].triangles[j].v1;
				edge3.v2 = objects[i].triangles[j].v2;
				edge3.normal = objects[i].triangles[j].normal;
				triangleEdges.push_back(edge3);

				for( uint k=0; k<triangleEdges.size(); ++k )
				{
					CheckCountourList(triangleEdges[k], contourList);
				}			
				caps.push_back(objects[i].triangles[j]);

				vec3 extrudedv0, extrudedv1, extrudedv2;
				
				extrudedv0 = objects[i].triangles[j].v0 + ExtrdudeMagnitude * (objects[i].triangles[j].v0 - lightPos);
				extrudedv1 = objects[i].triangles[j].v1 + ExtrdudeMagnitude * (objects[i].triangles[j].v1 - lightPos);	
				extrudedv2 = objects[i].triangles[j].v2 + ExtrdudeMagnitude * (objects[i].triangles[j].v2 - lightPos);	

				caps.push_back(Triangle(extrudedv0, extrudedv1, extrudedv2, objects[i].triangles[j].color));	
			}
			triangleEdges.clear();	
		}
	}

	for( uint k=0; k<contourList.size(); ++k )
	{
		Quad quad;
		quad.v1 = contourList[k].v1;
		quad.v2 = contourList[k].v2;
		quad.v3 = contourList[k].v1 + ExtrdudeMagnitude * (contourList[k].v1 - lightPos);
		quad.v4 = contourList[k].v2 + ExtrdudeMagnitude * (contourList[k].v2 - lightPos);
		quad.normal = contourList[k].normal;
		quads.push_back(quad);
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

vec3 ComputeAmbientLight( vec3 currentReflactance )
{
	return currentReflactance * indirectLightPowerPerArea;
}

void PixelShader( Pixel& p, vec3 currentColor, vec3 currentNormal, vec3 currentReflactance )
{
	int x = p.x;
	int y = p.y;

	if( x < SCREEN_HEIGHT && y < SCREEN_WIDTH )
	{
		if( depthBuffering == 1 )
		{
			if( stencilBuffering == 1 )
			{
				if( p.zinv > depthBuffer[y][x] )
				{
					depthBuffer[y][x] = p.zinv;
					vec3 R = ComputePixelReflectedLight(p, currentNormal, currentReflactance);
					currentColor = currentColor * (R + frameBuffer[y][x]);
					// PutPixelSDL( screen, x, y, currentColor );
					#ifndef POSTPROCESSING
						PutPixelSDL( screen, x, y, currentColor);
					#else 
						image[x][y] = currentColor;
					#endif
				}
			}
			else 
			{
				if( p.zinv > depthBuffer[y][x] )
				{
					depthBuffer[y][x] = p.zinv;
					vec3 R = ComputeAmbientLight(currentReflactance);
					frameBuffer[y][x] = R * currentColor;
				}
			}
		}
		else
		{
			if( p.zinv > depthBuffer[y][x] )
			{
				if(renderFrontFaces)
				{
					stencilBuffer[y][x]++;
					frameBuffer[y][x] = vec3(0.f, 0.f, 0.f);
				}
				else 
				{
					stencilBuffer[y][x]--;
					frameBuffer[y][x] = vec3(0.f, 0.f, 0.f);
				}
			}
		}
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

void RenderTriangles(vector<Triangle> triangles)
{
	//#pragma omp parallel for
	for( uint i=0; i<triangles.size(); ++i )
	{
		if(triangles[i].culled)
			continue;
		vector<Vertex> vertices(3);

		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;

		vec3 currentNormal = triangles[i].normal;
		vec3 currentReflactance = triangles[i].color;
		
		DrawPolygon(vertices, triangles[i].color, currentNormal, currentReflactance);
	}
}

void RenderQuad(Quad quad)
{
	vector<Vertex> vertices(3);

	vertices[0].position = quad.v1;
	vertices[1].position = quad.v2;
	vertices[2].position = quad.v3;

	DrawPolygon(vertices, pink, vec3(0.0f,0.0f,0.0f), pink);

	vertices[0].position = quad.v3;
	vertices[1].position = quad.v2;
	vertices[2].position = quad.v4;

	DrawPolygon(vertices, pink, vec3(0.0f,0.0f,0.0f), pink);
}

void RenderQuads(vector<Quad> quads)
{
	#pragma omp parallel for
	for( uint i=0; i<quads.size(); ++i)
	{
		if(renderFrontFaces)
		{
			if(dot(quads[i].normal, lightPos) >= 0.0f)
			{
				RenderQuad(quads[i]);
			}
		}
		else
		{
			if(dot(quads[i].normal, lightPos) < 0.0f)
			{
				RenderQuad(quads[i]);
			}

		}
	}
}

void RenderCap(Triangle cap)
{
	vector<Vertex> vertices(3);

	vertices[0].position = cap.v0;
	vertices[1].position = cap.v1;
	vertices[2].position = cap.v2;

	DrawPolygon(vertices, pink, vec3(0.0f,0.0f,0.0f), pink);
}

void RenderCaps(vector<Triangle> caps)
{
	#pragma omp parallel for
	for( uint i=0; i<caps.size(); ++i)
	{
		if(renderFrontFaces)
		{
			if(dot(caps[i].normal, lightPos) >= 0.0f)
			{
				RenderCap(caps[i]);
			}
		}
		else
		{
			if(dot(caps[i].normal, lightPos) < 0.0f)
			{
				RenderCap(caps[i]);
			}

		}
	}
}

void Draw() 
{
	SDL_FillRect( screen, 0, 0 );
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	// Initialise depth and stencil buffers
	#pragma omp parallel for
	for( uint i=0; i<SCREEN_HEIGHT; ++i )
	{
		for( uint j=0; j<SCREEN_WIDTH; ++j )
		{
			depthBuffer[i][j] = 0.0f;
			stencilBuffer[i][j] = 0.0f;
		}
	}

	vector<Quad> quads;
	vector<Triangle> caps;
	ComputeSilhouettes(sceneObjects, quads, caps);

	//Clipping setup
	float f = 251.0f;
	vec3 fVec = glm::normalize(vec3(0,0,1.0f)*cameraR);
	float near = cameraPos.z+fVec.z*0.01f, far = cameraPos.z+fVec.z*15.0f;
	float w = (float)SCREEN_WIDTH, h = (float)SCREEN_HEIGHT;
	// Perspective matrix transformation
	mat4 transform = glm::mat4(0.0f);
	// fovy version
	vec3 t(0.0f, -h/2.0f, f);
	vec3 b(0.0f, h/2.0f, f);
	float cy = dot(t,b)/(glm::length(t)*glm::length(b));
	float rfovy = acos(cy);
	float fovy = (180.0f/M_PI)*rfovy;
	float aspect = w/h;
	transform[0][0] = (1.0f/tan(rfovy/2.0f))/aspect;
	transform[1][1] = (1.0f/tan(rfovy/2.0f));
	transform[2][2] = far/(far-near);
	transform[3][2] = near*far/(far-near);
	transform[3][3] = -1.0f;
	
	for( size_t i = 0; i < triangles.size(); ++i )
	{
		triangles[i].culled = false;
	}

	// #pragma omp parallel for
	for( uint i=0; i<triangles.size(); ++i )
	{
		if (glm::dot((triangles[i].v0-cameraPos),triangles[i].normal)>0.0f)
		{
			triangles[i].culled = true;
		}

		vec3 new_v0;
		vec3 new_v1;
		vec3 new_v2;
		if(triangles[i].culled == false) {
			vec3 v0 = triangles[i].v0;
			vec3 v1 = triangles[i].v1;
			vec3 v2 = triangles[i].v2;

			// Go to view space, using this in a matrix form seems to produce different results
			v0 = (v0-cameraPos)*cameraR;
			v1 = (v1-cameraPos)*cameraR;
			v2 = (v2-cameraPos)*cameraR;
			
			// Map to clipping space (Hopefully this is right)
			vec4 tv0 = glm::vec4(v0.x, v0.y, v0.z, 1.0f);
			vec4 tv1 = glm::vec4(v1.x, v1.y, v1.z, 1.0f);
			vec4 tv2 = glm::vec4(v2.x, v2.y, v2.z, 1.0f);
			tv0 = tv0*transform;
			tv1 = tv1*transform;
			tv2 = tv2*transform;
			
			bool culled_v0 = insideCube(tv0);
			bool culled_v1 = insideCube(tv1);
			bool culled_v2 = insideCube(tv2);

			if (!culled_v0 || !culled_v1 || !culled_v2) {
				triangles[i].culled = true;
			}

			tv0 = tv0/tv0[3];
			tv1 = tv1/tv1[3];
			tv2 = tv2/tv2[3];

			//transform the coordinates to view space
			mat4 viewport_transform;

			viewport_transform[0][0] = 0.5;
			viewport_transform[0][1] = 0;
			viewport_transform[0][2] = 0;
			viewport_transform[0][3] = 0.5;

			viewport_transform[1][0] = 0;
			viewport_transform[1][1] = 0.5;
			viewport_transform[1][2] = 0;
			viewport_transform[1][3] = 0.5;

			viewport_transform[2][0] = 0;
			viewport_transform[2][1] = 0;
			viewport_transform[2][2] = 0.5;
			viewport_transform[2][3] = 0.5;

			viewport_transform[3][0] = 0;
			viewport_transform[3][1] = 0;
			viewport_transform[3][2] = 0;
			viewport_transform[3][3] = 1;

			//apply the viewport transform matrix
			vec4 view_tv0 = tv0 * viewport_transform;
			vec4 view_tv1 = tv1 * viewport_transform;
			vec4 view_tv2 = tv2 * viewport_transform;

			new_v0 = vec3(view_tv0[0], view_tv0[1], view_tv0[2]);
			new_v1 = vec3(view_tv1[0], view_tv1[1], view_tv1[2]);
			new_v2 = vec3(view_tv2[0], view_tv2[1], view_tv2[2]);
		}
	}

	// Disable writing to stencil buffer
	stencilBuffering = 0;
	RenderTriangles(triangles);

	// Disable writing to depth buffer
	depthBuffering = 0;

	// Draw front faces
	renderFrontFaces = 1;
	RenderQuads(quads);
	RenderCaps(caps);

	// Draw back faces
	renderFrontFaces = 0;
	RenderQuads(quads);
	RenderCaps(caps);

	// Enable stencil and buffer testing
	stencilBuffering = 1;
	depthBuffering = 1;
	
	//Clear depth buffer
	#pragma omp parallel for
	for( uint i=0; i<SCREEN_HEIGHT; ++i )
	{
		for( uint j=0; j<SCREEN_WIDTH; ++j )
		{
			depthBuffer[i][j] = 0.0f;
		}
	}

	// Render scene preforming depth and stencil test
	RenderTriangles(triangles);

	#ifdef POSTPROCESSING
		ApplyPostprocessing();
	#endif	

    if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);
    SDL_UpdateRect( screen, 0, 0, 0, 0 );

}