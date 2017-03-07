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
};

struct Adjacencies
{
	vector<Triangle> triangles;
	vector<vec3> v1;
	vector<vec3> v2;
};

/*CLASSES*/
//Used to define an object of triangular surfaces
class Object
{
public:
	std::vector<Triangle> triangles;
	std::vector<Edge> silouhette;

	Object( std::vector<Triangle> triangles )
		:triangles(triangles)
	{

	}
};

/* GLOBAL VARIABLES */
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
int t;
mat3 cameraR;
vec3 cameraPos( 0, 0, -3.001 );
float f = 2.0f;
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
vec3 lightPower = 5.1f*vec3( 1.0f, 1.0f, 1.0f );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

/*DIRECTIONAL LIGHT*/
vec3 lightDirection(1.f, -0.5f, -0.7f);

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
void ComputeSilhouettes( vector<Object>& objects );
void ComputeAdjacencies( Triangle triangle, vector<Triangle> triangles, Adjacencies& adjacencies );
//Debugging functions
void DrawSilouhetteEdges(vector<Object> sceneObjects);


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

	ComputeSilhouettes(sceneObjects);
	//debugging
	DrawSilouhetteEdges(sceneObjects);

	while( NoQuitMessageSDL() )
	{
		// Update();
		// Draw();
		ComputeSilhouettes(sceneObjects);
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
		cout << lightPos.x << "\n";
		cout << lightPos.y << "\n";
		cout << lightPos.z << "\n";

	}
}


//Debugging functions
void drawEdge(vector<Pixel> line)
{
	cout << "Starting drawing" << endl;
	for( uint i=0; i<line.size(); ++i )
	{
		PutPixelSDL(screen, line[i].x, line[i].y, pink);
		cout << "Pixel put!" << endl;
	}
	cout << "Finished edge!" << endl;
}

void DrawSilouhetteEdges(vector<Object> sceneObjects)
{
	Vertex vect1, vect2;
	Pixel p1;
	Pixel p2;
	for( uint i=0; i<sceneObjects.size(); ++i )
	{
		for( uint j=0; j<sceneObjects[i].silouhette.size(); ++j )
		{
			vect1.position = sceneObjects[i].silouhette[j].v1;
			vect2.position = sceneObjects[i].silouhette[j].v2;

			VertexShader(vect1, p1);
			VertexShader(vect2, p2);
			
			vec2 _p1 = vec2(p1.x, p1.y);
			vec2 _p2 = vec2(p2.x, p2.y);

			ivec2 delta = abs(_p1 - _p2);
			int pixels = max(delta.x, delta.y) + 1;
			vector<Pixel> result(pixels);

			Interpolate(p1, p2, result);
			drawEdge(result);
		}
	}
	cout << "Finished drawing silouhettes!" << endl;
}
//--


bool isAdjacent(Triangle triangle1, Triangle triangle2, vector<vec3>& v)
{
	int commonVertices = 0;

	if(triangle1.v0 == triangle2.v0 || triangle1.v0 == triangle2.v1 || triangle1.v0 == triangle2.v2)
	{
		commonVertices++;
		v.push_back(triangle1.v0);
	}
	if(triangle1.v1 == triangle2.v0 || triangle1.v1 == triangle2.v1 || triangle1.v1 == triangle2.v2)
	{
		commonVertices++;
		v.push_back(triangle1.v1);
	}
	if( (triangle1.v2 == triangle2.v0 || triangle1.v2 == triangle2.v1 || triangle1.v2 == triangle2.v2) )
	{
		commonVertices++;
		v.push_back(triangle1.v2);
	}

	return (commonVertices == 2);
}

bool isSameTriangle(Triangle triangle1, Triangle triangle2)
{
	return (triangle1.v0 == triangle2.v0 && triangle1.v1 == triangle2.v1 && triangle1.v2 == triangle2.v2);
}

void ComputeAdjacencies( Triangle triangle, vector<Triangle> triangles, Adjacencies& adjacencies )
{	
	vector<vec3> v;
	for( uint i=0; i<triangles.size(); ++i )
	{
		if(!isSameTriangle(triangle, triangles[i]))
		{
			if(isAdjacent(triangle, triangles[i], v))
			{
				adjacencies.triangles.push_back(triangles[i]);
				adjacencies.v1.push_back(v[0]);
				adjacencies.v2.push_back(v[1]);
			}
		}
	}
}

void ComputeSilhouettes( vector<Object>& objects )
{
	Adjacencies adjacencies;

	for( uint i=0; i<objects.size(); ++i )
	{
		for( uint j=0; j<objects[i].triangles.size(); ++j )
		{
			ComputeAdjacencies(objects[i].triangles[j], objects[i].triangles, adjacencies);
			vec3 currentNormal = objects[i].triangles[j].normal;
			float dir1 = dot(normalize(currentNormal), lightPos);
			if(dir1 >= 0)
			{
				for( uint k=0; k<adjacencies.triangles.size(); ++k )
				{	
					vec3 normal = adjacencies.triangles[k].normal;
					float dir2 = dot(normalize(normal), lightPos);
					if(dir1 * dir2 < 0)
					{
						Edge edge;
						edge.v1 = adjacencies.v1[k];
						edge.v2 = adjacencies.v2[k];
						objects[i].silouhette.push_back(edge);
					}
				} 
				adjacencies.triangles.clear();
				adjacencies.v1.clear();
				adjacencies.v2.clear();
			}
		}
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
		// vec3 R = ComputePixelDirectionalLight(p, currentNormal, currentReflactance); //direct light
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
