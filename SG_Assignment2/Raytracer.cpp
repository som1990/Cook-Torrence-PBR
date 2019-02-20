/* 
Cook Torrence Model Physically Based Renderer
Author: Soumitra Goswami
Date: 10/22/2017

Description: This takes the sample raytracer provided by Dr.Keyser and adds physically based materials and 
			physically based Point Light with a spectral distribution of a sunlight.
			The goal of this program was to Use the Spectral Distribute of light and have physical materials 
			by correctly scattering light using the Cook-Torrence BRDF Model.

Input:		This program requires Tab Delimited .csv files of the following:
			1) Sunlight Spectral Distribution Data http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html
			2) CIE 1931 Table of wavelength->XYZ conversion. http://www.cvrl.org/cmfs.htm
			3) Info on Material Refractive Indicies https://refractiveindex.info/

Output:		Processes the information to provide an accurate physical model of the material and renders the image.

*/

#include<GL/glut.h>
#include<fstream>
#include<math.h>
#include<string>
#include<vector>
#include<iostream>
#include<map>
#include<glm/glm.hpp>

using namespace std;

#define ImageW 600
#define ImageH 400

#define PI 3.14159265359

float framebuffer[ImageH][ImageW][3];

struct Coord3D { double x, y, z; };

glm::dvec3 Light;
glm::dvec3 SphCent[6];
double SphRad[6];

//Material Properties
string SphMat[6];
double SphRoughness[6];
double SphReflectivity[6];

//Generating necessary Dictionaries
map<string, map<int, glm::dvec2>> material;
map<float, double> lightInt;
map<float, double> lightInt2;
map<int, Coord3D> waveXYZ;

//MATH FUNCTIONS
//----------------------------------------------------
// Normalizes the vector passed in
void normalize(double& x, double& y, double& z) {
	float temp = sqrt(x*x + y*y + z*z);
	if (temp > 0.0) {
		x /= temp;
		y /= temp;
		z /= temp;
	}
	else {
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
}

// Returns dot product of two vectors
double dot(double x1, double y1, double z1, double x2, double y2, double z2) {
	return (x1*x2 + y1*y2 + z1*z2);
}

// Returns angle between two vectors
float angle(double x1, double y1, double z1, double x2, double y2, double z2) {
	normalize(x1, y1, z1);
	normalize(x2, y2, z2);
	return  acos(dot(x1, y1, z1, x2, y2, z2));
}

//TABLE FUNCTIONS 
//-------------------------------------------------------------------

//Reads the csv table and stores it in a Dictionary with the Wavelengths as key
map<int,glm::dvec2> fillConductor(string file)
{
	ifstream inFile(file);
	int lambda;
	map<int, glm::dvec2> temp;
	while (!inFile.eof())
	{
		inFile >> lambda;
		double n, k;
		inFile >> n;
		inFile >> k;
		temp[lambda] = glm::vec2(n, k);
	}
	inFile.close();
	return temp;
}

/*
void generateMat()
Description:	Generates the supported Materials and stores them in a dictionary.
Input : none.
Output: Creates a library of Materials that can be used by the program.
				
*/

void generateMat()
{

	map<int, glm::dvec2> temp = fillConductor("Aluminum.csv");
	material["Alluminium"] = temp;
	
	temp = fillConductor("Gold.csv");
	material["Gold"] = temp;

	temp = fillConductor("Brass.csv");
	material["Brass"] = temp;

	temp = fillConductor("Silk.csv");
	material["Silk"] = temp;

/*	for (auto& x : temp)
	{
		cout << x.first << ": " << x.second.x << " " << x.second.y << endl;
	}
*/
}

/*
void createTable()
Description:	Generates the supported Light and stores its spectral Distribution in a dictionary. Also Generates
				the CIE 1931 Wavelength->XYZ table to be used by the renderer.
Input : none.
Output: Creates a library of Spectral Distribution of the light and Wavelength->XYZ table.

*/

void createTable()
{
	ifstream inFile("SunlightIntensity.csv");
	float lambda;
	while (!inFile.eof()) {
		inFile >> lambda;
		inFile >> lightInt[lambda];
		inFile >> lightInt2[lambda];
	}
	inFile.close();

	for (auto& x : lightInt)
	{
		//cout << x.first << ":" << x.second << endl;
 	}

	inFile.open("WavelengthXYZ.csv", fstream::in);
	int lam;
	Coord3D temp;
	while (!inFile.eof()) {
		inFile >> lam;
		inFile >> temp.x;
		inFile >> temp.y;
		inFile >> temp.z;
		waveXYZ[lam] = temp;
	}
	inFile.close();
	for (auto& x : waveXYZ)
	{
		//cout << x.first << " : " << x.second.x << " " << x.second.y << " " << x.second.z << endl;
	}
}


//COLOR SPACE FUNCTION
//--------------------------------------------------------------------------------

/*

void convertsRGB(double,double,double,float&,float&,float&)
Description:	Converts the XYZ value of the color to the "standard sRGB" Color Space. You can find more information about
				other color spaces here. http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
Input : double X,Y,Z - The Color in the XYZ Color Space
Output: float& R,G,B - The resultant Color in standard sRGB Color Space.

*/
void convertsRGB(double X, double Y, double Z, float& R, float& G, float& B)
{
	//GLM MATRIX FUNCTION IS COLUMN-MAJOR ORDER WHILE C++ IS ROW-MAJOR ORDER. So the Matrix is TRANSPOSED as compared to the one
	//online to get the right calculations.
	glm::dmat3 conversion =	{ 3.2404542,	-0.9692660,	 0.0556434,
							 -1.5371385,	 1.8760108, -0.2040259,
							 -0.4985314,	 0.0415560,  1.0572252
							};
	glm::dvec3 xyz = { X,Y,Z };
	glm::dvec3 rgb = conversion*xyz;
	R = rgb.r;
	G = rgb.g;
	B = rgb.b;
}

//COOK - TORENCE BRDF FUNCTIONS
//--------------------------------------------------------------------------
/*

void calcF(double,double,double)
Description:	Calculates the Fresnel component based on s-polarized and p-polarized lights. This is material and
				wavelength dependant.
Input :		double n - the Refractive Index of the material for a particular wavelength.
			double k - the coeficcient of extinction of the material for a particular wavelength
			double NdotL - The cosine of the incident angle. Incident Angle is the Angle between the Light Ray and the Normal of the surface.
Output:		double F- The Fresnel component of the Cook-Torrence BRDF Model.

*/

double calcF(double n, double k, double NdotL)
{
	double a = n*n + k*k;
	double cosTheta2 = NdotL*NdotL;
	double b = 2 * n*NdotL;

	double Rp = (a - b + cosTheta2) / (a + b + cosTheta2);
	double Rs = ((a*cosTheta2) - b + 1) / ((a*cosTheta2) + b + 1);

	double F = (Rp + Rs) * 0.5;

	return F;
}

/*

void calcF0(double,double,double)
Description:	Calculates the Fresnel component of a ray normal to surface at the point of intersection. This generates the diffuse
				component of the ray
Input :		double n - the Refractive Index of the material for a particular wavelength.
			double k - the coeficcient of extinction of the material for a particular wavelength
			double NdotL(NOT USED) - The cosine of the incident angle. Incident Angle is the Angle between the Light Ray and the Normal of the surface.
Output:		double F0- The Fresnel component of the Cook-Torrence BRDF Model normal to the surface.

*/

double calcF0(double n, double k, double NdotL)
{
	double F0 = (k*k + (n - 1)*(n - 1)) / (k*k + (n + 1)*(n + 1));
	return F0;
}


/*

void calcFSchlick(double,double,double)
Description:	Calculates the Fresnel component based on the Schlick approximation. Is dependant on material properties, incident ray and wavelength.
Input :		double n - the Refractive Index of the material for a particular wavelength.
			double k - the coeficcient of extinction of the material for a particular wavelength
			double NdotL - The cosine of the incident angle. Incident Angle is the Angle between the Light Ray and the Normal of the surface.
Output:		double F- The Fresnel component of the Cook-Torrence BRDF Model.

*/

double calcFSchlick(double n, double k, double NdotL)
{
	double F0 = (k*k + (n - 1)*(n - 1)) / (k*k + (n + 1)*(n + 1));
	double F = F0 + (1 - F0)*glm::pow((1.0 - NdotL), 5);
	return F;
}

/*

void calcD(double,double)
Description:	Calculates the Distribution Component of the specular BRDF. This is the heart of the specular Component of the shader. The Cook-Torrence
				utilizes the BeckMann Distribution Function to generate it's specular highlight.
Input :			double m - the Roughness of the surface. This specifies the distribution of the reflected rays vectors. The higher the Roughness the more
							spread out the rays will be giving it a more diffused look.
				double NdotH - the consine of the angle between Normal and Half-Angle. Half-Angle is the ray bisecting the LightVector and ViewVector. H=V+L/2
Output:			double D- Amount of Distribution of specular rays . Generates the specular Highlight you see.

*/

double calcD(double m, double NdotH)
{
	//BeckMann Distribution Function
	double tanA = glm::tan(glm::acos(NdotH));
	double c = 1.0 / (m*m*glm::pow(NdotH, 4));
	double exponent = -glm::pow((tanA / m), 2);
	double D = c*glm::exp(exponent);
	return D;
}

/*

void calcD(double,double,double,double)
Description:	Calculates the Geometric Component of the Specular BRDF. This specifies whether the point is in shadow of itself.
Input :			double NdotH - the cosine of the angle between Surface Normal and Half-Angle. Half-Angle is the ray bisecting the LightVector and ViewVector. H=V+L/2
				double NdotV - the cosine of the angle between Surface Normal and View Direction.
				double NdotL - the cosine of the angle between Surface Normal and Light Direction.
				double VdotH - the cosine of the angle between View Direction and Half Angle.
Output:			double G- Returns Geometric Component which tells the renderer how much of the object is in shadow/masked with itself.

*/


double calcG( double NdotH, double NdotV, double NdotL, double VdotH)
{
	double a = 2 * NdotH*NdotV / VdotH;
	double b = 2 * NdotH*NdotL / VdotH;

	double G = glm::min(1.0, glm::min(a, b));
	return G;
}

/*

void calcBRDF(double,double,int, double,double, double, double)
Description:	Collects the Fresnel, Distribution and Geometric component of the Specular BRDF as well as the Diffuse Component of the BRDF to generate the 
				resultant reflectance of the particular wavelength after it scatters on the surface. This will always be from 0-1. Tells the renderer how much
				energy  of Light the material reflected and/or Transmitted. In our case we only calculate Reflectance.  
Input :			double n -
				double k - 
				int SphNum - Tells the renderer which Material to get it's properties from.
				double NdotH - the cosine of the angle between Surface Normal and Half-Angle. Half-Angle is the ray bisecting the LightVector and ViewVector. H=V+L/2
				double NdotV - the cosine of the angle between Surface Normal and View Direction.
				double NdotL - the cosine of the angle between Surface Normal and Light Direction.
				double VdotH - the cosine of the angle between View Direction and Half Angle.
				double s - Reflectivity of the surface is taken from the global Material Reflectivity list. This specifies the distribution of Specular Component and 
							Diffuse Component of the Material.
				double m - the Roughness of the surface is taken from the global Material Roughness list. This specifies the distribution of the reflected rays vectors. 
							The higher the Roughness the more spread out the rays will be giving it a more diffused look.

Output:			Returns the Reflectance of the Material . Percentage of Energy Reflected out of the surface based on physical Properties of the material.

*/

double calcBRDF(double n, double k, int SphNum, double NdotH, double NdotV, double NdotL, double VdotH)
{
	double s = SphReflectivity[SphNum];
	double d = 1.0 - s;				//energy conservation
	double Rs = 0;
	double Rd = 0;
	
	double m = SphRoughness[SphNum];
	
	double F = calcFSchlick(n,k,NdotL);
	double D = calcD(m, NdotH);
	double G = calcG(NdotH, NdotV, NdotL, VdotH);
	
	double a = F / PI;
	double b = D / NdotL;
	double c = G / NdotV;

	Rs = a * b * c;
	Rd = calcF0(n,k,NdotL);
	

//	cout  << " F/PI: " << a << " D/NDotL: " << b << " G/NdotV: " << c << endl;
	double R = s*Rs + d*Rd;

	return R;
}

//PIXEL COLOR FUNCTIONS
//---------------------------------------------------------------------------------
// Get Color for a point on the surface
void GetColor(glm::dvec3 view,   // Normalized Vector pointing FROM eye TO surface
	glm::dvec3 normal, // Normalized Vector giving surface normal
	glm::dvec3 light,  // Normalized Vector pointing FROM surface TO light
	int SphNum,     // Sphere Number (0-5)
	float& R,       // Return these values for surface color.
	float& G,
	float& B) {

	double X=0, Y=0, Z=0;
	
	//Generating Half-Angle . 
	//We invert the view because V needs to be FROM Surface TO Eye.
	glm::dvec3 h;
	h.x = (-view.x + light.x) * 0.5;
	h.y = (-view.y + light.y) * 0.5;
	h.z = (-view.z + light.z) * 0.5;
	normalize(h.x, h.y, h.z);

	//Generating the required angles.
	double NdotH = glm::dot(normal,h); // Normal to Light
	double VdotH = glm::dot(-view, h); // View to Half-Angle
	double NdotV = glm::dot(normal, -view); // View to Normal
	double NdotL = glm::dot(normal, light); // Normal to Light
	map<int, glm::dvec2> mat = material.find(SphMat[SphNum])->second; //Getting the Material
	int lambda = 390;
	//shoot each wavelength ray from the light based on its spectral distribution. In our case we are using Sunlight.
	for (int i = lambda; i <= 750; i=i+5)
	{
		//Gets the intensity of this particular wavelength based on the spectral Distribution Function of the light.
		map<float, double>::iterator it = lightInt2.find(i);
		double lInt = it->second;
		//gets the resultant xbar, ybar and zbar of the current wavelength.
		map<int, Coord3D>::iterator xyz = waveXYZ.find(i);
		double xBar = xyz->second.x;
		double yBar = xyz->second.y;
		double zBar = xyz->second.z;
		map<int, glm::dvec2>::iterator mIt = mat.find(i);
		double n = mIt->second.x;
		double k = mIt->second.y;
//		cout << "i:" << lambda;
		double R = calcBRDF(n, k, SphNum, NdotH, NdotV, NdotL, VdotH);
//		cout << "i:" << lambda << " R: " << R << endl;
		
		//COOK-TORRENCE Basic Reflectance Model 
		double Ir = lInt*NdotL*0.03*R;
		
		X += xBar*Ir;
		Y += yBar*Ir;
		Z += zBar*Ir;
	}
	//Converts the resultant XYZ to Pixel Data our Renderer can display. 
	convertsRGB(X, Y, Z, R, G, B);
	
}

//OPENGL FUNCTIONS
//------------------------------------------------------------------------------------
// Draws the scene
void drawit(void)
{
	glDrawPixels(ImageW, ImageH, GL_RGB, GL_FLOAT, framebuffer);
	glFlush();
}

// Sets pixel x,y to the color RGB
void setFramebuffer(int x, int y, float R, float G, float B)
{
	if (R <= 1.0)
		if (R >= 0.0)
			framebuffer[y][x][0] = R;
		else
			framebuffer[y][x][0] = 0.0;
	else
		framebuffer[y][x][0] = 1.0;
	if (G <= 1.0)
		if (G >= 0.0)
			framebuffer[y][x][1] = G;
		else
			framebuffer[y][x][1] = 0.0;
	else
		framebuffer[y][x][1] = 1.0;
	if (B <= 1.0)
		if (B >= 0.0)
			framebuffer[y][x][2] = B;
		else
			framebuffer[y][x][2] = 0.0;
	else
		framebuffer[y][x][2] = 1.0;
}

//RAYTRACE FUNCTION
//----------------------------------------------------------------------------
void display(void)
{
	int i, j, k;
	float R, G, B;
	glm::dvec3 refpt;
	glm::dvec3 view;
	glm::dvec3 normal;
	glm::dvec3 light;
	glm::dvec3 intpt;
	double xstep = 12.0 / ImageW;
	double ystep = 8.0 / ImageH;
	double t;
	double a, b, c;
	int intsphere;

	refpt.x = -6.0 + xstep / 2.0;
	refpt.y = -4.0 + ystep / 2.0;
	refpt.z = -10.0;

	for (i = 0; i<ImageW; i++, refpt.x += xstep) {
		for (j = 0; j<ImageH; j++, refpt.y += ystep) {
			// Compute the view vector
			view.x = refpt.x; view.y = refpt.y; view.z = refpt.z;
			normalize(view.x, view.y, view.z);

			// Find intersection with sphere (if any) - only 1 sphere can intesect.
			intsphere = -1;
			for (k = 0; (k<6) && (intsphere == -1); k++) {
				a = 1.0;  // Since normalized;
				b = 2.0*view.x*(-SphCent[k].x) + 2.0*view.y*(-SphCent[k].y) + 2.0*view.z*(-SphCent[k].z);
				c = SphCent[k].x*SphCent[k].x + SphCent[k].y*SphCent[k].y + SphCent[k].z*SphCent[k].z -
					SphRad[k] * SphRad[k];
				if ((b*b - 4 * a*c) >= 0.0) {  // We have an intersection with that sphere
											   // Want nearest of two intersections
					t = (-b - sqrt(b*b - 4 * a*c)) / 2.0;
					intsphere = k;
				}
			}

			if (intsphere != -1) { // We had an intersection with a sphere
				intpt.x = t*view.x; intpt.y = t*view.y; intpt.z = t*view.z;
				normal.x = (intpt.x - SphCent[intsphere].x) / SphRad[intsphere];
				normal.y = (intpt.y - SphCent[intsphere].y) / SphRad[intsphere];
				normal.z = (intpt.z - SphCent[intsphere].z) / SphRad[intsphere];
				normalize(normal.x, normal.y, normal.z);

				light.x = Light.x - intpt.x;
				light.y = Light.y - intpt.y;
				light.z = Light.z - intpt.z;
				normalize(light.x, light.y, light.z);
				GetColor(view, normal, light, intsphere, R, G, B);

			}
			else {
				R = G = B = 0.0;
			}
			setFramebuffer(i, j, R, G, B);
		}
		refpt.y = -4.0 + ystep / 2.0;
	}

	drawit();
}

//INITIALIZATION FUNCTIONS

void init(void)
{
	int i, j;

	// Initialize framebuffer to clear
	for (i = 0; i<ImageH; i++) {
		for (j = 0; j<ImageW; j++) {
			framebuffer[i][j][0] = 0.0;
			framebuffer[i][j][1] = 0.0;
			framebuffer[i][j][2] = 0.0;
		}
	}
	createTable();
	generateMat();
	// Create Sphere data
	SphCent[0].x = -3.0;
	SphCent[0].y = 1.5;
	SphCent[0].z = -10.0;
	SphCent[1].x = 0.0;
	SphCent[1].y = 1.5;
	SphCent[1].z = -10.0;
	SphCent[2].x = 3.0;
	SphCent[2].y = 1.5;
	SphCent[2].z = -10.0;
	SphCent[3].x = -3.0;
	SphCent[3].y = -1.5;
	SphCent[3].z = -10.0;
	SphCent[4].x = 0.0;
	SphCent[4].y = -1.5;
	SphCent[4].z = -10.0;
	SphCent[5].x = 3.0;
	SphCent[5].y = -1.5;
	SphCent[5].z = -10.0;

	//Initialize the ratio of Diffuse to Specular here
	SphReflectivity[0] = SphReflectivity[1] = SphReflectivity[2] = 0.6;
	SphReflectivity[3] = SphReflectivity[4] = 0.75;
	SphReflectivity[5] = 0.4;

	//Initialize the roughness here.
	SphRoughness[0] = 0.1;
	SphRoughness[1] = 0.5;
	SphRoughness[2] = 0.8;
	SphRoughness[3] = 0.4;
	SphRoughness[4] = 0.4;
	SphRoughness[5] = 0.4;

	//Initialize the Material Type Here.
	SphMat[0] = SphMat[1] = SphMat[2] = "Alluminium";
	SphMat[3] = "Brass";
	SphMat[4] = "Gold";
	SphMat[5] = "Silk";

	for (i = 0; i < 6; i++)
	{
		SphRad[i] = 1.0;
		
		//SphMat[i] = "Gold";
	}

	// Set Light Position
	Light.x = 10.0;
	Light.y = 6.0;
	Light.z = 0.0;

	// Eye is at origin, looking down -z axis, y axis is up,
	// Looks at 8x6 window centered around z = -10.
}


int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(ImageW, ImageH);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Soumitra Goswami - 641 Assignment 2");
	init();
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}