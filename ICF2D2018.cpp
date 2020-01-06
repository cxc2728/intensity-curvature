//  This file contains sourcecode distributed as freeware. 
//  The intellectual property of the sourcecode is shown 
//  here to belong to Carlo Ciulla.

// Disclaimer: 

// The website here named www.sourcecodewebsiteCarloCiulla.com [1] does not intend 
// to convey the meaning of profit making for what pertains to the content
// provided. --->>> Instead, when the content is downloaded, the user(s) are
// kindly invited to donate money to charity organizations involved in 
// helping people in need of food and water. <<<---


// The Novel Re-sampling Locations have been sized to be a fraction of 
// the pixel size. The programs presented here confirm both concepts and 
// implications brought to knowledge through the unifying theory [1].

// Reference:

// [1] Carlo Ciulla "Improved Signal and Image Interpolation in Biomedical Applications: 
// The Case of Magnetic Resonance Imaging (MRI)." Medical Information Science 
// Reference - IGI Global Publisher - March 2009; ISBN: 978 - 160566202 - 2.

//  Project Title: Intensity-Curvature Functional: HPF-Filter2D, G42D, B32D, H32D, H42D, LGR2D, SRE2D
#define _CRT_SECURE_NO_WARNINGS

#include < iostream >
#include < fstream >
#include < string >
#include < io.h >
#include < dos.h >
#include < conio.h >
#include < stdlib.h >
#include < sstream >
#include < stdio.h >
#include < iomanip >
#include < istream >
#include < math.h >


#define ft_SCALE 255 
// this constant is necessary in order to control the 
// contrast brightness of the Fourier Transformations 
// (direct and inverse)


using namespace std;

void OnFourierTransform(char imageFilename[], int rcxres, int rcyres);
void OnInverseFourierTransform(char filename[], int rcyres, int rcxres);
void OnInverseFourierTransform(char filename1[], char filename2[], char filename3[], char filename4[], int rcxres, int rcyres, char filename5[]);


class SRE2D2013 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data {

		double **Signal; // pointer to the matrix entry 

		double **ICF_SRE2D; // pointer to the matrix entry

		double **SRE2D_ictBI; // pointer to the matrix entry

		double **SRE2D_ictAI; // pointer to the matrix entry

		double **ICF_G42D; // pointer to the matrix entry

		double **g42D_ictBI; // pointer to the matrix entry

		double **g42D_ictAI; // pointer to the matrix entry

		double **ICF_B32D; // pointer to the matrix entry

		double **b32D_ictBI; // pointer to the matrix entry

		double **b32D_ictAI; // pointer to the matrix entry

		double **ICF_H32D; // pointer to the matrix entry

		double **H32D_ictBI; // pointer to the matrix entry

		double **H32D_ictAI; // pointer to the matrix entry

		double **ICF_H42D; // pointer to the matrix entry

		double **H42D_ictBI; // pointer to the matrix entry

		double **H42D_ictAI; // pointer to the matrix entry

		double **ICF_LGR2D; // pointer to the matrix entry

		double **LGR2D_ictBI; // pointer to the matrix entry

		double **LGR2D_ictAI; // pointer to the matrix entry

		double **theta_x; // pointer to the matrix entry

		double **theta_y; // pointer to the matrix entry

		double **omega_f; // pointer to the matrix entry

		double **Filter2DH; // pointer to the matrix entry

	}*pointer; // pointer to the matrices

public:

	SRE2D2013(int x, int y) : n1(x), n2(y) { };// constructor 
	
	void allocateData();

	void save();

	~SRE2D2013() { } // destructor

};

void SRE2D2013::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;
			
	 pointer->Signal = new double*[this->n1];

	 pointer->ICF_SRE2D = new double*[this->n1];

	 pointer->SRE2D_ictBI = new double*[this->n1];

	 pointer->SRE2D_ictAI = new double*[this->n1];

	 pointer->ICF_G42D = new double*[this->n1];
 
	 pointer->g42D_ictBI = new double*[this->n1];

	 pointer->g42D_ictAI = new double*[this->n1];

	 pointer->ICF_B32D = new double*[this->n1];
 
	 pointer->b32D_ictBI = new double*[this->n1];

	 pointer->b32D_ictAI = new double*[this->n1];

	 pointer->ICF_H32D = new double*[this->n1];
 
	 pointer->H32D_ictBI = new double*[this->n1];

	 pointer->H32D_ictAI = new double*[this->n1];

	 pointer->ICF_H42D = new double*[this->n1];
 
	 pointer->H42D_ictBI = new double*[this->n1];

	 pointer->H42D_ictAI = new double*[this->n1];

	 pointer->ICF_LGR2D = new double*[this->n1];
 
	 pointer->LGR2D_ictBI = new double*[this->n1];

	 pointer->LGR2D_ictAI = new double*[this->n1];

	 pointer->theta_x = new double*[this->n1];

	 pointer->theta_y = new double*[this->n1];

	 pointer->omega_f = new double*[this->n1];

	 pointer->Filter2DH = new double*[this->n1]; 


	 for( int v=0; v < this->n1; v++ ) { // (1)
		 
		 pointer->Signal[v] = new double[this->n2];

		 pointer->ICF_SRE2D[v] = new double[this->n2];

		 pointer->SRE2D_ictBI[v] = new double[this->n2];

		 pointer->SRE2D_ictAI[v] = new double[this->n2];

		 pointer->ICF_G42D[v] = new double[this->n2];

		 pointer->g42D_ictBI[v] = new double[this->n2];

		 pointer->g42D_ictAI[v] = new double[this->n2];

		 pointer->ICF_B32D[v] = new double[this->n2];

		 pointer->b32D_ictBI[v] = new double[this->n2];

		 pointer->b32D_ictAI[v] = new double[this->n2];

		 pointer->ICF_H32D[v] = new double[this->n2];

		 pointer->H32D_ictBI[v] = new double[this->n2];

		 pointer->H32D_ictAI[v] = new double[this->n2];

		 pointer->ICF_H42D[v] = new double[this->n2];

		 pointer->H42D_ictBI[v] = new double[this->n2];

		 pointer->H42D_ictAI[v] = new double[this->n2];

		 pointer->ICF_LGR2D[v] = new double[this->n2];

		 pointer->LGR2D_ictBI[v] = new double[this->n2];

		 pointer->LGR2D_ictAI[v] = new double[this->n2];

		 pointer->theta_x[v] = new double[this->n2];

		 pointer->theta_y[v] = new double[this->n2];

		 pointer->omega_f[v] = new double[this->n2];

		 pointer->Filter2DH[v] = new double[this->n2];

	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for(int v=0; v < this->n1; v++ ) { // (a)

			for( int f=0; f < this->n2 ; f++ ) { // (b)
		 
			pointer->Signal[v][f] = (double)0.0;

			pointer->ICF_SRE2D[v][f] = (double)0.0;

			pointer->SRE2D_ictBI[v][f] = (double)0.0;

			pointer->SRE2D_ictAI[v][f] = (double)0.0;

			pointer->ICF_G42D[v][f] = (double)0.0;

			pointer->g42D_ictBI[v][f] = (double)0.0;

			pointer->g42D_ictAI[v][f] = (double)0.0;

			pointer->ICF_B32D[v][f] = (double)0.0;

			pointer->b32D_ictBI[v][f] = (double)0.0;

			pointer->b32D_ictAI[v][f] = (double)0.0;

			pointer->ICF_H32D[v][f] = (double)0.0;

			pointer->H32D_ictBI[v][f] = (double)0.0;

			pointer->H32D_ictAI[v][f] = (double)0.0;

			pointer->ICF_H42D[v][f] = (double)0.0;

			pointer->H42D_ictBI[v][f] = (double)0.0;

			pointer->H42D_ictAI[v][f] = (double)0.0;

			pointer->ICF_LGR2D[v][f] = (double)0.0;

			pointer->LGR2D_ictBI[v][f] = (double)0.0;

			pointer->LGR2D_ictAI[v][f] = (double)0.0;

			pointer->theta_x[v][f] = (double)0.0;

			pointer->theta_y[v][f] = (double)0.0;

			pointer->omega_f[v][f] = (double)0.0;

			pointer->Filter2DH[v][f] = (double)0.0;

		    } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


void SRE2D2013::save() { // saveImages

	FILE * savedata;
	char outputFile[128];
	
	sprintf(outputFile, "%s","Signal.img");

	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->Signal[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ICF_SRE2D.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ICF_SRE2D[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","SRE2D_ICTBI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->SRE2D_ictBI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","SRE2D_ICTAI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->SRE2D_ictAI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ICF_G42D.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ICF_G42D[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","g42D_ICTBI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->g42D_ictBI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","g42D_ICTAI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->g42D_ictAI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ICF_B32D.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ICF_B32D[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","b32D_ICTBI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->b32D_ictBI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","b32D_ICTAI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->b32D_ictAI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ICF_H32D.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ICF_H32D[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","H32D_ICTBI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->H32D_ictBI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","H32D_ICTAI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->H32D_ictAI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ICF_H42D.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ICF_H42D[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","H42D_ICTBI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->H42D_ictBI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","H42D_ICTAI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->H42D_ictAI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ICF_LGR2D.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ICF_LGR2D[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","LGR2D_ICTBI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->LGR2D_ictBI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","LGR2D_ICTAI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->LGR2D_ictAI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","theta_x.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->theta_x[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","theta_y.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->theta_y[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","omega_f.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->omega_f[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","Filter2DH.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->Filter2DH[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

} // saveImages



int main ( int argc, char * argv[] ) {

	char outputFile[128]="ICF2D2018.log";

	FILE * savedata;

	double MAX = 5000000000000000000.0;

if (argc < 12) { std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the pixel size along the X direction (double)" << endl;
				 std::cout << "Please enter the pixel size along the Y direction (double)" << endl;
				 std::cout << "Please enter the misplacement along the X direction (double)" << endl;
				 std::cout << "Please enter the misplacement along the Y direction (double)" << endl;
				 std::cout << "Please enter the XY rotation angle (double)" << endl;
				 std::cout << "Please enter the cut-off frequency of the traditional HPF (double) > 0 and < 1" << endl;
				 std::cout << "Please enter the a constant (double)" << endl;
				 std::cout << "Please enter the b constant (double)" << endl;
				 std::cout << endl;
				 exit(0); }

else { // run the program (begin)
	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // processing (begin)

	int n1 = atoi(argv[3]);
	int n2 = atoi(argv[2]);

	double XPixelSize = atof(argv[4]);
	double YPixelSize = atof(argv[5]);

	double x_misplacement_X = atof(argv[6]);
	double y_misplacement_Y = atof(argv[7]);

	double theta = atof(argv[8]);

	char imageFileName[128];
	
	sprintf(imageFileName, "%s", argv[1]);
	
	double cutOFF_frequency = atof(argv[9]);

	double the_A_constant = atof(argv[10]);
	double theBconstant = atof(argv[11]);

	std::cout << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
	std::cout << "The pixel size along the X direction is: " << atof(argv[4]) << endl;
	std::cout << "The pixel size along the Y direction is: " << atof(argv[5]) << endl;
	std::cout << "The XY rotation angle is: " << atof(argv[8]) << endl;
	std::cout << "The cutOFF frequency of the traditional HPF is: " << atof(argv[9]) << endl;
	std::cout << "The value of the a constant is: " << atof(argv[10]) << endl;
	std::cout << "The value of the b constant is: " << atof(argv[11]) << endl;


	
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The XY rotation angle is: ", theta);
	fprintf(savedata,"%s%lf\n", "The cutOFF frequency of the traditional HPF is: ", cutOFF_frequency);
	fprintf(savedata,"%s%lf\n", "The value of the a constant is: ", atof(argv[10]));
	fprintf(savedata,"%s%lf\n", "The value of the b constant is: ", atof(argv[11]));


    double misplacement_X = ((double)1.0 - ( cos( (double)theta ) + sin( (double)theta ) ) + x_misplacement_X);
    double misplacement_Y = ((double)1.0 - ( -sin( (double)theta ) + cos( (double)theta ) ) + y_misplacement_Y);

      misplacement_X = ((double)misplacement_X/XPixelSize);
      misplacement_Y = ((double)misplacement_Y/YPixelSize);

	  //////////////////***********//////////////////////
	  // Above formula scales the misplacement to the  //
	  // pixel size the same way the following formula //
	  // would do: (min - misplacement)/(min - max)    //  
	  //////////////////***********//////////////////////


	SRE2D2013 SRE(n1,n2);

	SRE.allocateData();

	/// read image file (begin)
	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		exit(0);

	} else { // else

	double number;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
			
		fread(&number,sizeof(double),1,pf);
		
		SRE.pointer->Signal[i1][i2] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl;

    // compute omega_f & theta_x & theta_y (begin)
	double omega;
	 // 11-10-2017 Calculate gradients (begins)
	for (int i1=0; i1 < n1-1; i1++) {// x dim

	    for (int i2=0; i2 < n2-1; i2++) { // y dim


		omega = ((double) SRE.pointer->Signal[i1+1][i2+1] + SRE.pointer->Signal[i1][i2] - 
		                  SRE.pointer->Signal[i1+1][i2] - SRE.pointer->Signal[i1][i2+1] );


		SRE.pointer->theta_x[i1][i2] = ( (double) SRE.pointer->Signal[i1+1][i2] - SRE.pointer->Signal[i1][i2] ); 
		
		SRE.pointer->theta_y[i1][i2] = ( (double) SRE.pointer->Signal[i1][i2+1] - SRE.pointer->Signal[i1][i2] ); 

		SRE.pointer->omega_f[i1][i2] = ((double) SRE.pointer->Signal[i1+1][i2+1] + SRE.pointer->Signal[i1][i2] - 
												 SRE.pointer->Signal[i1+1][i2] - SRE.pointer->Signal[i1][i2+1] );
	   

		} // y dim
        
	}  // x dim
	 // 11-10-2017 Calculate gradients (ends)
	// compute omega_f & theta_x & theta_y (end)   	
	
	
	double MAX = 5000000000000000000.0;
	double max=-MAX;
	double min=MAX;

	
	///  scale data (begins)
	/// compute max and min of data (begin)
	// theta_x
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->theta_x[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->theta_x[i1][i2];
              
		if( SRE.pointer->theta_x[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->theta_x[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->theta_x[i1][i2] = (double)0.0;

           else SRE.pointer->theta_x[i1][i2] = (double) ft_SCALE * (double) fabs( (min - SRE.pointer->theta_x[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim 
	// scale (end)
	// scale data (ends)

    max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->theta_y[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->theta_y[i1][i2];
              
		if( SRE.pointer->theta_y[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->theta_y[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->theta_y[i1][i2] = (double)0.0;

           else SRE.pointer->theta_y[i1][i2] = (double) ft_SCALE * (double) fabs ( (min - SRE.pointer->theta_y[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim 

    max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->omega_f[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->omega_f[i1][i2];
              
		if( SRE.pointer->omega_f[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->omega_f[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->omega_f[i1][i2] = (double)0.0;

           else SRE.pointer->omega_f[i1][i2] = (double) ft_SCALE * (double) fabs ( (min - SRE.pointer->omega_f[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim 
	// scale (end)
	// compute omega_f & theta_x & theta_y (end) 

    std::cout << "omega_f & theta_x & theta_y calculated" << endl;

	std::cout << "Now calculating the traditional High Pass Filter" << endl;

	int XNEI = (int)2;
	double FW = ((double)n1*n2);
    int n7 = ( (int)floor( (double)n1/2.0) );  
	int n8 = ( (int)floor( (double)n2/2.0) );

	// build the FILTER function www.wikipedia.org & www.dspguide.com/ch14/5.htm (begin)
	double highPassFilter;
	double pi = 3.141592;
	double x = 0.0, y = 0.0, w = 0.0;
	double deltaT = 1.0;

	for (int pp =-n7+XNEI; pp < n7-XNEI; pp++) {

            for (int qq =-n8; qq < n8; qq++) {


                x = SRE.pointer->Signal[pp + n7][qq + n8];
				
				w = SRE.pointer->Signal[pp + n7 - 1][qq + n8];

				y = SRE.pointer->Filter2DH[pp + n7 - 1][qq + n8];

				double RC = (double) 1.0 / ((double) 2.0 * pi * cutOFF_frequency);

				highPassFilter = ((double) RC) / ((double) RC + deltaT);

				SRE.pointer->Filter2DH[pp + n7][qq + n8] = ((double) highPassFilter * y ) + 
					                                      (((double)highPassFilter) * ((double)x - w ));
				          
			}

        } 
	
	
	// scale data (begin)
	// Filter2DH
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->Filter2DH[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->Filter2DH[i1][i2];
              
		if( SRE.pointer->Filter2DH[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->Filter2DH[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->Filter2DH[i1][i2] = (double)0.0;

           else SRE.pointer->Filter2DH[i1][i2] = (double) ft_SCALE * (double) fabs( (min - SRE.pointer->Filter2DH[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim 
	
	///  scale data (end)
	
	 std::cout << "High Pass Filter calculated" << endl;
	// build the FILTER function www.wikipedia.org & www.dspguide.com/ch14/5.htm (end)

	// calculate ICF_SRE2D (begin)
	std::cout << "Compute ICF of the SRE2D model function" << endl;
	
	double k = 0.0, s = 0.0;
	for (int i1=0; i1 < n1-1; i1++) {// x dim

	    for (int i2=0; i2 < n2-1; i2++) { // y dim

		  
        k = (double) SRE.pointer->Signal[i1][i2] * misplacement_X * misplacement_Y + 
				   ( misplacement_X * misplacement_Y * misplacement_X / 2.0) * 
				   ( SRE.pointer->Signal[i1+1][i2] - SRE.pointer->Signal[i1][i2] ) + 
				   ( misplacement_X * misplacement_Y * misplacement_Y / 2.0) * 
				   ( SRE.pointer->Signal[i1][i2+1] - SRE.pointer->Signal[i1][i2] ) +
				   ( misplacement_X * misplacement_X * misplacement_Y * misplacement_Y / 4.0) * 
				   ( SRE.pointer->Signal[i1+1][i2+1] + SRE.pointer->Signal[i1][i2] - SRE.pointer->Signal[i1+1][i2] - SRE.pointer->Signal[i1][i2+1] );


        s = (double) SRE.pointer->Signal[i1][i2] * misplacement_X * misplacement_Y; 

		SRE.pointer->SRE2D_ictBI[i1][i2] = (double)s;

		SRE.pointer->SRE2D_ictAI[i1][i2] = (double)k;
       
		   if ( (double)k != 0.0 ) SRE.pointer->ICF_SRE2D[i1][i2] = (double)s/k;  
		   else SRE.pointer->ICF_SRE2D[i1][i2] = (double) 0.0;

		} // y dim
        
	}  // x dim


	std::cout << "ICF of the SRE2D model function calculated" << endl;
	// calculate ICF_SRE2D(end)

	// scale ICF_SRE2D (begin) 
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->ICF_SRE2D[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->ICF_SRE2D[i1][i2];
              
		if( SRE.pointer->ICF_SRE2D[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->ICF_SRE2D[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->ICF_SRE2D[i1][i2] = (double)0.0;

           else SRE.pointer->ICF_SRE2D[i1][i2] = (double) ft_SCALE * (double) fabs( (min - SRE.pointer->ICF_SRE2D[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim
	// scale ICF_SRE2D (end)


	/// allocate memory & store ALPHAs (begin)
	double * ALPHA3 = 0;
	double * ALPHA2 = 0;
	double * ALPHA1 = 0;

	/// allocate ALPHAs (begin)
	if ((ALPHA3 = (double *) calloc( n1*n2, sizeof(double)) ) == NULL)
	{

	std::cout << "Not enough memory to allocate Image data, Exit." << endl;
    fprintf(savedata,"%s\n", "Not enough memory to allocate Image data, Exit.");

	fclose(savedata);	
	SRE.~SRE2D2013();
	exit(0);

	} else { // else

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		int index = ((int) i2*n1 + i1); 

		*(ALPHA3 + index) = (double)0.0;
    	

		} // y dim
        
	}  // x dim 

	}//else

	if ((ALPHA2 = (double *) calloc( n1*n2, sizeof(double)) ) == NULL)
	{

	std::cout << "Not enough memory to allocate Image data, Exit." << endl;
    fprintf(savedata,"%s\n", "Not enough memory to allocate Image data, Exit.");

	fclose(savedata);
	
	free(ALPHA3);
	SRE.~SRE2D2013();
	exit(0);

	} else { // else


	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
		
		int index = ((int) i2*n1 + i1); 

		*(ALPHA2 + index) = (double)0.0;


		} // y dim
        
	}  // x dim 

	}//else

	if ((ALPHA1 = (double *) calloc( n1*n2, sizeof(double)) ) == NULL)
	{

	std::cout << "Not enough memory to allocate Image data, Exit." << endl;
    fprintf(savedata,"%s\n", "Not enough memory to allocate Image data, Exit.");

	fclose(savedata);
	
	free(ALPHA3);
	free(ALPHA2);
	SRE.~SRE2D2013();
	exit(0);

	} else { // else


	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
		
		int index = ((int) i2*n1 + i1); 

		*(ALPHA1 + index) = (double)0.0;


		} // y dim
        
	}  // x dim 

	}//else
	/// allocate ALPHAs (end)
	
	 // make the convolutions (begin)
	    for (int i = 1; i < n1-1; i++) {
            for (int j = 1; j < n2-1; j++) {

			int index = ((int) j*n1 + i);

			*(ALPHA1 + index) =((double) -SRE.pointer->Signal[i+1][j] - SRE.pointer->Signal[i-1][j] - SRE.pointer->Signal[i-1][j+1] + SRE.pointer->Signal[i+1][j-1] - 
				                          SRE.pointer->Signal[i][j-1] - SRE.pointer->Signal[i-1][j-1] - SRE.pointer->Signal[i][j+1] + SRE.pointer->Signal[i+1][j+1]);

			*(ALPHA2 + index) =((double) -SRE.pointer->Signal[i+1][j] - SRE.pointer->Signal[i-1][j] - SRE.pointer->Signal[i-1][j+1] + SRE.pointer->Signal[i+1][j-1] - 
				                          SRE.pointer->Signal[i][j-1] - SRE.pointer->Signal[i-1][j-1] - SRE.pointer->Signal[i][j+1] + SRE.pointer->Signal[i+1][j+1]);
            }  
                                 
        }
        

	    for (int i = 1; i < n1-1; i++) {
            for (int j = 1; j < n2-1; j++) {

			int index = ((int) j*n1 + i);

			  *(ALPHA3 + index) = ((double) -SRE.pointer->Signal[i+1][j] + SRE.pointer->Signal[i-1][j] - SRE.pointer->Signal[i][j+1] + SRE.pointer->Signal[i][j-1]);
            }
         } // make the convolutions (ends)   
	/// allocate memory & store ALPHAs (end)

	// save ALPHAs (begins)
	sprintf(outputFile, "%s","ALPHA1.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output ALPHA1.img output file" << endl;

	} else  { // (save)

	int index;

	for( int v=0; v < n1; v++ ) { // (a)

		for( int f=0; f < n2; f++ ) {

			index = ((int) f*n1 + v);
	
			fwrite(& *(ALPHA1 + index),sizeof(double),1,savedata);

		}

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ALPHA2.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output ALPHA2.img output file" << endl;

	} else  { // (save)

	int index;

	for( int v=0; v < n1; v++ ) { // (a)

		for( int f=0; f < n2; f++ ) {

			index = ((int) f*n1 + v);
	
			fwrite(& *(ALPHA2 + index),sizeof(double),1,savedata);

		}

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ALPHA3.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output ALPHA3.img output file" << endl;

	} else  { // (save)

	int index;

	for( int v=0; v < n1; v++ ) { // (a)

		for( int f=0; f < n2; f++ ) {

			index = ((int) f*n1 + v);
	
			fwrite(& *(ALPHA3 + index),sizeof(double),1,savedata);

		}

	} // (a)

	fclose(savedata);

	} // (save)
	// save ALPHAs (ends)

	std::cout << "Compute ICF of the G4 model function" << endl;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		int index = ((int) i2*n1 + i1); 

		double alpha3 =  (double)*(ALPHA3 + index);

		double alpha2 =  (double)*(ALPHA2 + index);

		double x5 = (double)pow((double)misplacement_X, 5.0);

		double x4 = (double)pow((double)misplacement_X, 4.0);

		double x3 = (double)pow((double)misplacement_X, 3.0);

		double x2 = (double)pow((double)misplacement_X, 2.0);

		double x1 = (double)misplacement_X;


		double y5 = (double)pow((double)misplacement_Y, 5.0);

		double y4 = (double)pow((double)misplacement_Y, 4.0);

		double y3 = (double)pow((double)misplacement_Y, 3.0);

		double y2 = (double)pow((double)misplacement_Y, 2.0);

		double y1 = (double)misplacement_Y;

		/// Calculation of the Intensity-curvature Functional (begin) ///
		double E0 = (double) 4.0 * the_A_constant * SRE.pointer->Signal[i1][i2] * x1 * y1 * (alpha2 + 2.0 * alpha3);

		double M4 = ((double) 1.0/5.0 * y1 * x5 + 1.0/2.0 * x4 * y2 + 2.0/3.0 * x3 * y3 + 1.0/2.0 * x2 * y4 + 1.0/5.0 * x1 * y5);

		double M3 = ((double)(1.0/4.0 * y1 * x4 + 1.0/6.0 * x3 * y2 + 1.0/3.0 * x3 * y2 + 
			                  1.0/3.0 * x2 * y3 + 1.0/6.0 * x2 * y3 + 1.0/4.0 * x1 * y4));

		double M2 = ((double)(1.0/3.0 * y1 * x3 + 1.0/2.0 * x2 * y2 + 1.0/3.0 * x1 * y3));

		double M1 = ((double) 1.0/2.0 * y1 * x2 + 1.0/2.0 * x1 * y2);

		double EIN = (double) 4.0 * (SRE.pointer->Signal[i1][i2] * (alpha2 * the_A_constant * (6.0 * (y1 * x2/2.0 + x1 * y2/2.0) + x1 * y1) +
			                  2.0 * alpha3 * the_A_constant * x1 * y1) + 6.0 * M4 * (alpha2 * the_A_constant * alpha2 * the_A_constant) + M3 *
							  (4.0 * (alpha2 * the_A_constant * alpha2 * the_A_constant) + 2.0 * alpha3 * alpha2 * the_A_constant * the_A_constant +
							  6.0 * alpha3 * alpha2 * the_A_constant * the_A_constant) + M2 * (2.0 * (alpha2 * the_A_constant * alpha2 * the_A_constant) + 
							  alpha3 * alpha2 * the_A_constant * the_A_constant + 13.0 * alpha3 * alpha2 * the_A_constant * the_A_constant + 2.0 * 
							  (alpha3 * the_A_constant * alpha3 * the_A_constant)) + M1 * (25.0/4.0 * (alpha2 * the_A_constant * alpha2 * the_A_constant) +
							  1.0/2.0 * alpha3 * alpha2 * the_A_constant * the_A_constant + 8.0 * alpha3 * alpha2 * the_A_constant * the_A_constant + 
							  4.0 * (alpha3 * the_A_constant * alpha3 * the_A_constant)) + ((alpha2 * the_A_constant * alpha2 * the_A_constant) + 
							  2.0 * alpha3 * alpha2 * the_A_constant * the_A_constant + alpha3 * alpha2 * the_A_constant * the_A_constant + 2.0 * 
							  (alpha3 * the_A_constant * alpha3 * the_A_constant) ) * x1 * y1);


		SRE.pointer->g42D_ictBI[i1][i2] = (double)E0;

		SRE.pointer->g42D_ictAI[i1][i2] = (double)EIN;
		
		if ( EIN != 0.0 ) SRE.pointer->ICF_G42D[i1][i2] = ((double)E0/EIN);
		else			  SRE.pointer->ICF_G42D[i1][i2] = (double)0.0;
		/// Calculation of the Intensity-curvature Functional (end) ///

		} // y dim
        
	}  // x dim	
	
	// scale ICF_G42D (begin) 
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->ICF_G42D[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->ICF_G42D[i1][i2];
              
		if( SRE.pointer->ICF_G42D[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->ICF_G42D[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->ICF_G42D[i1][i2] = (double)0.0;

           else SRE.pointer->ICF_G42D[i1][i2] = (double) ft_SCALE * (double) fabs( (min - SRE.pointer->ICF_G42D[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim
	// scale ICF_G42D (end)

	std::cout << "ICF of the G4 model function calculated" << endl;

	/// allocate memory & store OMEGAs (begin)
	double * OMEGAA = 0;
	double * OMEGAB = 0;
	double * OMEGAAB = 0;

	/// allocate OMEGAs (begin)
	if ((OMEGAA = (double *) calloc( n1*n2, sizeof(double)) ) == NULL)
	{

	std::cout << "Not enough memory to allocate Image data, Exit." << endl;
    fprintf(savedata,"%s\n", "Not enough memory to allocate Image data, Exit.");

	fclose(savedata);
	
	SRE.~SRE2D2013();
	exit(0);

	} else { // else

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		int index = ((int) i2*n1 + i1); 

		*(OMEGAA + index) = (double)0.0;
    	

		} // y dim
        
	}  // x dim 

	}//else

	if ((OMEGAB = (double *) calloc( n1*n2, sizeof(double)) ) == NULL)
	{

	std::cout << "Not enough memory to allocate Image data, Exit." << endl;
    fprintf(savedata,"%s\n", "Not enough memory to allocate Image data, Exit.");

	fclose(savedata);
	
	free(OMEGAA);
	SRE.~SRE2D2013();
	exit(0);

	} else { // else


	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
		
		int index = ((int) i2*n1 + i1); 

		*(OMEGAB + index) = (double)0.0;


		} // y dim
        
	}  // x dim 

	}//else

	if ((OMEGAAB = (double *) calloc( n1*n2, sizeof(double)) ) == NULL)
	{

	std::cout << "Not enough memory to allocate Image data, Exit." << endl;
    fprintf(savedata,"%s\n", "Not enough memory to allocate Image data, Exit.");

	fclose(savedata);
	
	free(OMEGAA);
	free(OMEGAB);
	SRE.~SRE2D2013();
	exit(0);

	} else { // else


	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
		
		int index = ((int) i2*n1 + i1); 

		*(OMEGAAB + index) = (double)0.0;


		} // y dim
        
	}  // x dim 

	}//else
	/// allocate OMEGAs (end)

	 // make the convolutions (begin) 
	 for (int i = 1; i < n1-1; i++) {
            for (int j = 1; j < n2-1; j++) {

			int index = ((int) j*n1 + i);
    
			*(OMEGAA + index) = ((double) -SRE.pointer->Signal[i+1][j] - SRE.pointer->Signal[i-1][j] - SRE.pointer->Signal[i-1][j+1] + SRE.pointer->Signal[i+1][j-1] - 
				                           SRE.pointer->Signal[i][j-1] - SRE.pointer->Signal[i-1][j-1] - SRE.pointer->Signal[i][j+1] + SRE.pointer->Signal[i+1][j+1]);
         
							   // ALPHA2
			                   /* ((double) -SRE.pointer->Signal[i+1][j] - SRE.pointer->Signal[i-1][j] - SRE.pointer->Signal[i-1][j+1] + SRE.pointer->Signal[i+1][j-1] - 
				                             SRE.pointer->Signal[i][j-1] - SRE.pointer->Signal[i-1][j-1] - SRE.pointer->Signal[i][j+1] + SRE.pointer->Signal[i+1][j+1]); */
            }  
                                 
        }
        
        for (int i = 1; i < n1-1; i++) {
            for (int j = 1; j < n2-1; j++) {

			int index = ((int) j*n1 + i);
    
			*(OMEGAB + index) = ((double) -SRE.pointer->Signal[i+1][j] - SRE.pointer->Signal[i-1][j] - SRE.pointer->Signal[i-1][j+1] + SRE.pointer->Signal[i+1][j-1] - 
				                           SRE.pointer->Signal[i][j-1] - SRE.pointer->Signal[i-1][j-1] - SRE.pointer->Signal[i][j+1] + SRE.pointer->Signal[i+1][j+1]);

			                 // ALPHA2
                             /* ((double) -SRE.pointer->Signal[i+1][j] - SRE.pointer->Signal[i-1][j] - SRE.pointer->Signal[i-1][j+1] + SRE.pointer->Signal[i+1][j-1] - 
				                           SRE.pointer->Signal[i][j-1] - SRE.pointer->Signal[i-1][j-1] - SRE.pointer->Signal[i][j+1] + SRE.pointer->Signal[i+1][j+1]); */
            }  
                                 
        }
        
        for (int i = 1; i < n1-1; i++) {
            for (int j = 1; j < n2-1; j++) {

              int index = ((int) j*n1 + i);

              *(OMEGAAB + index) = ((double) -SRE.pointer->Signal[i+1][j] + SRE.pointer->Signal[i-1][j] - 
				                              SRE.pointer->Signal[i][j+1] + SRE.pointer->Signal[i][j-1]); 

			                    // ALPHA3
			                    /* ((double) -SRE.pointer->Signal[i+1][j] + SRE.pointer->Signal[i-1][j] - 
								              SRE.pointer->Signal[i][j+1] + SRE.pointer->Signal[i][j-1]); */
         
            }
         }  

	// make the convolutions (ends)
	/// allocate memory & store OMEGAs (end)

	// save OMEGAs (begins)
	sprintf(outputFile, "%s","OMEGAA.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)

	int index;

	for( int v=0; v < n1; v++ ) { // (a)

		for( int f=0; f < n2; f++ ) {

			index = ((int) f*n1 + v);
	
			fwrite(& *(OMEGAB + index),sizeof(double),1,savedata);

		}

	} // (a)

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","OMEGAB.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)

	int index;

	for( int v=0; v < n1; v++ ) { // (a)

		for( int f=0; f < n2; f++ ) {

			index = ((int) f*n1 + v);
	
			fwrite(& *(OMEGAB + index),sizeof(double),1,savedata);

		}

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","OMEGAAB.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)

	int index;

	for( int v=0; v < n1; v++ ) { // (a)

		for( int f=0; f < n2; f++ ) {

			index = ((int) f*n1 + v);
	
			fwrite(& *(OMEGAAB + index),sizeof(double),1,savedata);

		}

	} // (a)

	fclose(savedata);

	} // (save)
	// save OMEGAs (ends)

	std::cout << "Compute ICF of the B3 model function" << endl;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		int index = ((int) i2*n1 + i1); 

		double omegaA =  (double)*(OMEGAA + index);

		double omegaB =  (double)*(OMEGAB + index);
		
		double omegaAB = (double)*(OMEGAAB + index);

		/// Calculation of the Intensity-Curvature Functional (begin) ///
		double E0 = (double)SRE.pointer->Signal[i1][i2] * 
			        (2.0 * the_A_constant * omegaA + 2.0 * theBconstant * omegaB);

		double EIN = (double)SRE.pointer->Signal[i1][i2] * 
			         (omegaA * (2.0 * the_A_constant * misplacement_X * misplacement_Y) + 
					  omegaAB * (the_A_constant * misplacement_X * misplacement_Y * misplacement_Y) +
					  omegaB * (2.0 * theBconstant * misplacement_X * misplacement_Y) +
					  omegaAB * (theBconstant * misplacement_X * misplacement_X * misplacement_Y) +
					  2.0 * omegaAB * (the_A_constant * misplacement_X * misplacement_X * misplacement_Y +
                                       theBconstant * misplacement_X * misplacement_Y * misplacement_Y)

					 ) + (omegaA * (the_A_constant * pow(misplacement_X, 3.0) * misplacement_Y/3.0 + 
					                misplacement_X * misplacement_Y * misplacement_Y/2.0) * omegaA * 2.0 * the_A_constant)
						 + 
						 (omegaA * (the_A_constant * the_A_constant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/3.0 + 
					                the_A_constant * misplacement_X * pow(misplacement_Y, 3.0) * 2.0/3.0) * omegaAB
						 ) +
						 (omegaA * (the_A_constant * pow(misplacement_X, 3.0) * misplacement_Y/3.0 + 
					                misplacement_X * pow(misplacement_Y, 2.0)/2.0) * omegaB * 2.0 * theBconstant
						 ) +
						 (omegaA * (the_A_constant * theBconstant * pow(misplacement_X, 4.0) * misplacement_Y/2.0 + 
					                theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0) * omegaAB
						 ) +
						 (omegaA * (the_A_constant * the_A_constant * pow(misplacement_X, 4.0) * misplacement_Y/2.0 + 
					                the_A_constant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0 +
									the_A_constant * theBconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/3.0 + 
									theBconstant * misplacement_X * pow(misplacement_Y, 3.0) * 2.0/3.0 ) * 2.0 * omegaAB
						 ) + 

						 (omegaB * (theBconstant * pow(misplacement_Y, 3.0) * misplacement_X/3.0 - 
					                pow(misplacement_X, 2.0) * misplacement_Y/2.0) * omegaA * 2.0 * the_A_constant)
						 + 
						 (omegaB * (the_A_constant * theBconstant * pow(misplacement_Y, 4.0) * misplacement_X/2.0 - 
					                the_A_constant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0) * omegaAB)
						 + 
						 (omegaB * (theBconstant * pow(misplacement_Y, 3.0) * misplacement_X/3.0 - 
					                pow(misplacement_X, 2.0) * misplacement_Y/2.0) * omegaB * 2.0 * theBconstant)
						 + 
						 (omegaB * (theBconstant * theBconstant * pow(misplacement_Y, 3.0) * pow(misplacement_X, 2.0)/3.0 - 
					                theBconstant * pow(misplacement_X, 3.0) * misplacement_Y * 2.0/3.0) * omegaAB)
						 +
						 (omegaB * (the_A_constant * theBconstant * pow(misplacement_Y, 3.0) * pow(misplacement_X, 2.0)/3.0 - 
					                the_A_constant * pow(misplacement_X, 3.0) * misplacement_Y * 2.0/3.0 + 
									theBconstant * theBconstant * misplacement_X * pow(misplacement_X, 4.0)/2.0 -
									theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0) * 2.0 * omegaAB)
						 +

						 (omegaAB * (the_A_constant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/6.0 + 
					                 theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 3.0)/6.0 ) * 
									 omegaA * 2.0 * the_A_constant 
						 ) + 
						 (omegaAB * (the_A_constant * the_A_constant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0) * 2.0/9.0 + 
					                 the_A_constant * theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 4.0)/4.0 ) * 
									 omegaAB 
						 ) + 
						 (omegaAB * (the_A_constant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/6.0 + 
					                 theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 3.0)/6.0 ) * 
									 omegaB * 2.0 * theBconstant
						 ) + 
						 (omegaAB * (the_A_constant * theBconstant * pow(misplacement_X, 4.0) * pow(misplacement_Y, 2.0)/4.0 + 
					                 theBconstant * theBconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0) * 2.0/9.0 ) * 
									 omegaAB
						 ) + 
						 (omegaAB * (the_A_constant * the_A_constant * pow(misplacement_X, 4.0) * pow(misplacement_Y, 2.0)/4.0 + 
					                 the_A_constant * theBconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0) * 2.0/9.0 +
									 the_A_constant * theBconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0) * 2.0/9.0 + 
									 theBconstant * theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 4.0)) * 2.0 * omegaAB 
						 ); 

		SRE.pointer->b32D_ictBI[i1][i2] = (double)E0;

		SRE.pointer->b32D_ictAI[i1][i2] = (double)EIN;

		if ( EIN != 0.0 ) SRE.pointer->ICF_B32D[i1][i2] = ((double)E0/EIN);
		else			  SRE.pointer->ICF_B32D[i1][i2] = (double)0.0;
		/// Calculation of the Intensity-Curvature Functional (end) ///

				} // y dim
        
	}  // x dim	
	
	
	// scale ICF_B32D (begin) 
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->ICF_B32D[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->ICF_B32D[i1][i2];
              
		if( SRE.pointer->ICF_B32D[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->ICF_B32D[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->ICF_B32D[i1][i2] = (double)0.0;

           else SRE.pointer->ICF_B32D[i1][i2] = (double) ft_SCALE * (double) fabs( (min - SRE.pointer->ICF_B32D[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim
	// scale ICF_B32D (end)
	
	std::cout << "ICF of the B3 model function calculated" << endl;

	std::cout << "Compute ICF of the H3 model function" << endl;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		 int index = ((int) i2*n1 + i1); 

	     double alpha1 =  (double)*(ALPHA1 + index);

	     double alpha2 =  (double)*(ALPHA2 + index);


		/// Calculation of the Intensity-Curvature Functional (begin) ///
		double E0 = (double)4.0 * misplacement_X * misplacement_Y * SRE.pointer->Signal[i1][i2] *
                    ((double) -4.0 * alpha1 * the_A_constant + 2.0 * the_A_constant * alpha2);

		double EIN = (double) 4.0 * ((double) -4.0 * alpha1 * the_A_constant + 2.0 * the_A_constant * alpha2) *
			       (  misplacement_X * misplacement_Y * SRE.pointer->Signal[i1][i2] +
					  alpha1 * (-2.0 * the_A_constant * (misplacement_X * misplacement_X * misplacement_X * misplacement_Y/3.0 + 
					  misplacement_X * misplacement_X * misplacement_Y * misplacement_Y/2.0 + 
					  misplacement_X * misplacement_Y * misplacement_Y * misplacement_Y/3.0) +
					  1.0/2.0 * misplacement_X * misplacement_Y * ( the_A_constant + 1.0 )) +
					  
					  alpha2 * (the_A_constant * (misplacement_X * misplacement_X * misplacement_X * misplacement_Y/3.0 +
					  misplacement_X * misplacement_X * misplacement_Y * misplacement_Y/2.0 + 
					  misplacement_X * misplacement_Y * misplacement_Y * misplacement_Y/3.0) - 
					  (2.0 * the_A_constant + 1.0/2.0) * (misplacement_X * misplacement_X * misplacement_Y/2.0 + 
					  misplacement_X * misplacement_Y * misplacement_Y/2.0) + 
					  (3.0/4.0) * misplacement_X * misplacement_Y * ( the_A_constant + 1.0 ) )   
					);

		SRE.pointer->H32D_ictBI[i1][i2] = (double)E0;

		SRE.pointer->H32D_ictAI[i1][i2] = (double)EIN;
		
		if ( EIN != 0.0 ) SRE.pointer->ICF_H32D[i1][i2] = ((double)E0/EIN);
		else			  SRE.pointer->ICF_H32D[i1][i2] = (double)0.0;
		/// Calculation of the Intensity-Curvature Functional (end) ///

			} // y dim
        
	}  // x dim	

	// scale ICF_H32D (begin) 
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->ICF_H32D[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->ICF_H32D[i1][i2];
              
		if( SRE.pointer->ICF_H32D[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->ICF_H32D[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->ICF_H32D[i1][i2] = (double)0.0;

           else SRE.pointer->ICF_H32D[i1][i2] = (double) ft_SCALE * (double) fabs( (min - SRE.pointer->ICF_H32D[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim
	// scale ICF_H32D (end)

	std::cout << "ICF of the H3 model function calculated" << endl;

	std::cout << "Compute ICF of the H4 model function" << endl;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			int index = ((int) i2*n1 + i1); 

			double alpha3 =  (double)*(ALPHA3 + index);

			double alpha2 =  (double)*(ALPHA2 + index);	

		/// Calculation of the Intensity-Curvature Functional (begin) ///
		double E0 = (double) 4.0 * misplacement_X * misplacement_Y * SRE.pointer->Signal[i1][i2] * 
			        ((double) - 2.0 * alpha3 + 2.0 * alpha2 );

		double EIN = (double) 4.0 * (SRE.pointer->Signal[i1][i2] * alpha3 * (3.0 * (pow(misplacement_X, 2.0) * 
			                  misplacement_Y/2.0 + misplacement_X * pow(misplacement_Y, 2.0)/2.0) - 2.0 * misplacement_X *
							  misplacement_Y)) + SRE.pointer->Signal[i1][i2] * alpha2 * ( (-pow(misplacement_X, 2.0) * 
							  misplacement_Y/2.0 - misplacement_X * pow(misplacement_Y, 2.0)/2.0 + 2.0 * misplacement_X * misplacement_Y)) +

							  4.0 * ( ( (3.0/2.0 * alpha3 * alpha3 - 1.0/2.0 * alpha3 * alpha2 - 3.0/6.0 * alpha3 * alpha2 + 1.0/6.0 * alpha2 * alpha2) *
							  (pow(misplacement_X, 5.0) * misplacement_Y/5.0 + 3.0/8.0 * pow(misplacement_X, 4.0) *
							  pow(misplacement_Y, 2.0) + pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0)/3.0 + pow(misplacement_X, 2.0) *
							  pow(misplacement_Y, 4.0)/8.0 + pow(misplacement_X, 4.0) * pow(misplacement_Y, 2.0)/8.0 +
							  pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0)/3.0 + 3.0/8.0 * pow(misplacement_X, 2.0) *
							  pow(misplacement_Y, 4.0) + misplacement_X * pow(misplacement_Y, 5.0)/5.0) +

							  (-4.0 * alpha3 * alpha3 + 2.0 * alpha3 * alpha2 + 4.0/3.0 * alpha3 * alpha2 - 4.0/3.0 * alpha2 * alpha2) * 
							  (pow(misplacement_X, 4.0) * misplacement_Y/4.0 + pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/2.0 +
							  pow(misplacement_X, 2.0) * pow(misplacement_Y, 3.0)/2.0 + misplacement_X * pow(misplacement_Y, 4.0)/4.0) +
							  (2.0 * alpha3 * alpha3 - 2.0 * alpha3 * alpha2 - 8.0 * alpha3 * alpha2 + 4.0 * alpha2 * alpha2) *
							  (pow(misplacement_X, 3.0) * misplacement_Y/3.0 + pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0 + 
							  misplacement_X * pow(misplacement_Y, 3.0)/3.0) + (2.0 * alpha3 * alpha3 - 2.0/3.0 * alpha3 * alpha2 + 8.0 * alpha3 *
							  alpha2 - 16.0/3.0 * alpha2 * alpha2) * (pow(misplacement_X, 2.0) * misplacement_Y/2.0 + misplacement_X *
							  pow(misplacement_Y, 2.0)/2.0) + (-4.0/3.0 * alpha3 * alpha3 + 4.0/3.0 * alpha3 * alpha2 - 8.0/3.0 * alpha3 * alpha2 +
							  8.0/3.0 * alpha2 * alpha2) * (misplacement_X * misplacement_Y)) );
         
		SRE.pointer->H42D_ictBI[i1][i2] = (double)E0;

		SRE.pointer->H42D_ictAI[i1][i2] = (double)EIN;

		if ( EIN != 0.0 ) SRE.pointer->ICF_H42D[i1][i2] = ((double)E0/EIN);
		else			  SRE.pointer->ICF_H42D[i1][i2] = (double)0.0;

		/// Calculation of the Intensity-Curvature Functional (end) ///

				} // y dim
        
	}  // x dim	

	// scale ICF_H42D (begin) 
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->ICF_H42D[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->ICF_H42D[i1][i2];
              
		if( SRE.pointer->ICF_H42D[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->ICF_H42D[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->ICF_H42D[i1][i2] = (double)0.0;

           else SRE.pointer->ICF_H42D[i1][i2] = (double) ft_SCALE * (double) fabs( (min - SRE.pointer->ICF_H42D[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim
	// scale ICF_H42D (end)
			
	std::cout << "ICF of the H4 model function calculated" << endl;

	std::cout << "Compute ICF of the Lagrange model function" << endl;
	
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
			
		int index = ((int) i2*n1 + i1); 

		double alpha3 =  (double)*(ALPHA3 + index);

		double alpha2 =  (double)*(ALPHA2 + index);

		/// Calculate the Intensity-Curvature Functional (begin) ///
		double E0 = (double) 4.0 * misplacement_X * misplacement_Y * SRE.pointer->Signal[i1][i2] * 
			        ((double) - 2.0 * alpha3 + 2.0 * alpha2 );

		double EIN = (double) 4.0 * (SRE.pointer->Signal[i1][i2] * alpha3 * (3.0 * (pow(misplacement_X, 2.0) * 
			                  misplacement_Y/2.0 + misplacement_X * pow(misplacement_Y, 2.0)/2.0) - 2.0 * misplacement_X *
							  misplacement_Y)) + SRE.pointer->Signal[i1][i2] * alpha2 * ( (-pow(misplacement_X, 2.0) * 
							  misplacement_Y/2.0 - misplacement_X * pow(misplacement_Y, 2.0)/2.0) + 2.0 * misplacement_X * misplacement_Y) +

							  4.0 * ( (3.0/2.0 * alpha3 * alpha3 - 1.0/2.0 * alpha3 * alpha2 - 3.0/6.0 * alpha3 * alpha2 + 1.0/6.0 * alpha2 * alpha2) *
							  (pow(misplacement_X, 5.0) * misplacement_Y/5.0 + 3.0/8.0 * pow(misplacement_X, 4.0) *
							  pow(misplacement_Y, 2.0) + pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0)/3.0 + pow(misplacement_X, 2.0) *
							  pow(misplacement_Y, 4.0)/8.0 + pow(misplacement_X, 4.0) * pow(misplacement_Y, 2.0)/8.0 +
							  pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0)/3.0 + 3.0/8.0 * pow(misplacement_X, 2.0) *
							  pow(misplacement_Y, 4.0) + misplacement_X * pow(misplacement_Y, 5.0)/5.0) +
							  (-4.0 * alpha3 * alpha3 + 2.0 * alpha3 * alpha2 + 10.0/3.0 * alpha3 * alpha2 - 8.0/6.0 * alpha2 * alpha2) * 
							  (pow(misplacement_X, 4.0) * misplacement_Y/4.0 + pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/2.0 +
							  pow(misplacement_X, 2.0) * pow(misplacement_Y, 3.0)/2.0 + misplacement_X * pow(misplacement_Y, 4.0)/4.0) +
							  (1.0/2.0 * alpha3 * alpha3 - 3.0/2.0 * alpha3 * alpha2 - 45.0/6.0 * alpha3 * alpha2 + 23.0/6.0 * alpha2 * alpha2) *
							  (pow(misplacement_X, 3.0) * misplacement_Y/3.0 + pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0 + 
							  misplacement_X * pow(misplacement_Y, 3.0)/3.0) + (4.0 * alpha3 * alpha3 - 2.0 * alpha3 * alpha2 + 40.0/6.0 * alpha3 *
							  alpha2 - 28.0/6.0 * alpha2 * alpha2) * (pow(misplacement_X, 2.0) * misplacement_Y/2.0 + misplacement_X *
							  pow(misplacement_Y, 2.0)/2.0) + (-2.0 * alpha3 * alpha3 + 2.0 * alpha2 * alpha2) * (misplacement_X * misplacement_Y));
                     
		SRE.pointer->LGR2D_ictBI[i1][i2] = (double)E0;

		SRE.pointer->LGR2D_ictAI[i1][i2] = (double)EIN;

		if ( EIN != 0.0 ) SRE.pointer->ICF_LGR2D[i1][i2] = ((double)E0/EIN);
		else			  SRE.pointer->ICF_LGR2D[i1][i2] = (double)0.0;
		/// Calculate the Intensity-Curvature Functional (end) ///
	
		} // y dim
        
	}  // x dim	

	// scale ICF_LGR2D (begin) 
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->ICF_LGR2D[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->ICF_LGR2D[i1][i2];
              
		if( SRE.pointer->ICF_LGR2D[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->ICF_LGR2D[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->ICF_LGR2D[i1][i2] = (double)0.0;

           else SRE.pointer->ICF_LGR2D[i1][i2] = (double) ft_SCALE * (double) fabs( (min - SRE.pointer->ICF_LGR2D[i1][i2]) / (min - max) );
			
		} // y dim
        
	}  // x dim
	// scale ICF_LGR2D (end)

	std::cout << "ICF of the Lagrange model function calculated" << endl;
	
	SRE.save(); // save all of the images

	
	std::cout << "Please Wait...Now running the Fourier Transformation on the Signal" << endl;

	OnFourierTransform(imageFileName, n2, n1);

	std::cout << "Please Wait...Now running the Fourier Transformation on the HP Filter" << endl;

	OnFourierTransform("Filter2DH.img", n2, n1);

	std::cout << "Please Wait...Now running the Fourier Transformation on the ICF_SRE2D" << endl;

	OnFourierTransform("ICF_SRE2D.img", n2, n1);

	std::cout << "Please Wait...Now running the Fourier Transformation on the ICF_B32D" << endl;

	OnFourierTransform("ICF_B32D.img", n2, n1);
	
	std::cout << "Please Wait...Now running the Fourier Transformation on the ICF_G42D" << endl;

	OnFourierTransform("ICF_G42D.img", n2, n1);
		
	std::cout << "Please Wait...Now running the Fourier Transformation on the ICF_H32D" << endl;

	OnFourierTransform("ICF_H32D.img", n2, n1);

	std::cout << "Please Wait...Now running the Fourier Transformation on the ICF_H42D" << endl;

	OnFourierTransform("ICF_H42D.img", n2, n1);

	std::cout << "Please Wait...Now running the Fourier Transformation on the ICF_LGR2D" << endl;

	OnFourierTransform("ICF_LGR2D.img", n2, n1);
	
	

	std::cout << "Please Wait...Now running the inverse Fourier Transformation procedures" << endl;

	char imageFileR[300];
	char imageFileI[300];

	sprintf(imageFileR, "%s%s", "K-SpaceR-", argv[1]);
	sprintf(imageFileI, "%s%s", "K-SpaceI-", argv[1]);

		OnInverseFourierTransform(imageFileR,
			                      "K-SpaceR-Filter2DH.img",
			                      imageFileI,
			                      "K-SpaceI-Filter2DH.img", n2, n1, "Filter2DH");
				
		OnInverseFourierTransform(imageFileR,
			                      "K-SpaceR-ICF_B32D.img",
			                      imageFileI,
			                      "K-SpaceI-ICF_B32D.img", n2, n1, "ICF_B32D");

		OnInverseFourierTransform(imageFileR,
			                      "K-SpaceR-ICF_G42D.img",
			                      imageFileI,
			                      "K-SpaceI-ICF_G42D.img", n2, n1, "ICF_G42D");

		OnInverseFourierTransform(imageFileR,
			                      "K-SpaceR-ICF_H32D.img",
			                      imageFileI,
			                      "K-SpaceI-ICF_H32D.img", n2, n1, "ICF_H32D");

		OnInverseFourierTransform(imageFileR,
			                      "K-SpaceR-ICF_H42D.img",
			                      imageFileI,
			                      "K-SpaceI-ICF_H42D.img", n2, n1, "ICF_H42D");

		OnInverseFourierTransform(imageFileR,
			                      "K-SpaceR-ICF_LGR2D.img",
			                      imageFileI,
			                      "K-SpaceI-ICF_LGR2D.img", n2, n1, "ICF_LGR2D");

		OnInverseFourierTransform(imageFileR,
			                      "K-SpaceR-ICF_SRE2D.img",
			                      imageFileI,
			                      "K-SpaceI-ICF_SRE2D.img", n2, n1, "ICF_SRE2D");


	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	fclose(savedata);
	
	free(ALPHA3);
	free(ALPHA2);
	free(ALPHA1);

	free(OMEGAA);
	free(OMEGAB);
	free(OMEGAAB);

	delete SRE.pointer;
	SRE.~SRE2D2013();
	} // processing (end)

	} // run the program (end)

	return 0;
} // end of main 



void OnFourierTransform(char imageFilename[], int rcxres, int rcyres)
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;

	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;

	double * kSpaceR = 0;
	double * kSpaceI = 0;
	double * Signal = 0;

	FILE * logfile;
	
	char logfilename[128]="Fourier-T.log";

  	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

		printf("%s\n %s\n" , "Unable to open log File", "Now Exit");

		exit(0);
	
	} else { // allocate memory 


	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}

	if ((Signal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);
		exit(0);

	}

	} // allocate memory 

	//// read image data and initialize pointers
	double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);
				
				*(kSpaceR+index) = (double) 0.0;

				*(kSpaceI+index) = (double) 0.0;

			}

		}

	FILE * pf;
	char SignalFilename[128];
	double readData;
	
	sprintf(SignalFilename, "%s", imageFilename);

	if ((pf = fopen(SignalFilename,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // read data


	for (i=0; i<rcyres; i++)
	{ ///read signal data
		for (j=0; j<rcxres; j++)
		{

			index = ((j*rcyres)+i);
          
            fread(&readData,sizeof(double),1,pf);

			*(Signal+index) = (double)readData;

		}
	} ///read signal data

	fprintf(logfile,"%s\n", "Signal Read in DOUBLE (64bits) format");

	fclose (pf);
	} // save data

	double phase, complexR, complexI;
	
	///// Fourier Transform //////
	for (i=0; i<NofYpixels; i++)
	{ ///calculate k-space data

		for (j=0; j<NofXpixels; j++)
		{

	
			dx = ((int) i - NofYpixels/2);
		    dy = ((int) j - NofXpixels/2);

			k2 = ((int)(dy*NofYpixels)+dx); 

			w = ((j*NofYpixels)+i);

			for (int s=0; s<NofYpixels; s++)
			{ ///calculate k-space data 
				for (int p=0; p<NofXpixels; p++)
				{ 
					

		     		ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);
 
				    k3 = ((int)(ds*NofXpixels)+dp); 

					t = ((p*NofYpixels)+s);

					phase = ((double) (2.0 * pi * k2 * k3) / (NofXpixels*NofYpixels));

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (begin)**/
					complexR = (double) cos( (double)phase ) + (double) sin( (double)phase ); 

					complexI = -(double) sin( (double)phase ) + (double) cos( (double)phase ); 
					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (end)**/
				
					*(kSpaceR+w) += (double) *(Signal+t) * (double) complexR;

					*(kSpaceI+w) -= (double) *(Signal+t) * (double) complexI;

			}

		}///calculate k-space data 

			    
		}
	} ///calculate k-space data

	///// Fourier Transform //////
	double savedata = 0.0;
	char FTfilename[128];

	sprintf(FTfilename, "%s%s", "K-SpaceR-", imageFilename);

    fprintf(logfile, "%s\t%s\n", "Now Saving K-Space Signal (Real) in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");


	 // FIFO memory deallocation method
 	 free(kSpaceR);
 	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(kSpaceR+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "K-Space Signal (Real) Saved");

	fclose (pf);
	} // save data



	sprintf(FTfilename, "%s%s", "K-SpaceI-", imageFilename);

    fprintf(logfile, "%s\t%s\n", "Now Saving K-Space Signal (Imaginary) in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(kSpaceI+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "K-Space Signal (Imaginary) Saved");

	fclose (pf);
	
	} // save data

	sprintf(FTfilename, "%s%s", "K-SpaceM-", imageFilename);

    fprintf_s(logfile, "%s\t%s\n", "Now Saving K-Space Magnitude of the Signal in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf_s(logfile, "%s\n", "Cannot open file to save K-Space Magnitude of the Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data	

		// save a zero image (begin)
		for (int s=0; s<NofYpixels; s++)
		{ 
			for (int p=0; p<NofXpixels; p++)
			{ 

			savedata = (double)0.0;
          
            fwrite(&savedata,sizeof(double),1,pf);

			}
		} // save a zero image (end)

	fclose(pf);
	
	}
		
	if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf_s(logfile, "%s\n", "Cannot open file to save K-Space Magnitude of the Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data
		
		// K-Space Magnitude (begin)
		for (int s=0; s<(int)NofYpixels; s++)
		{ 
			for (int p=0; p<(int)NofXpixels; p++)
			{ 
			
		
			index = ((p*NofYpixels)+s);

			savedata = (double) sqrt( (double)*(kSpaceR+index)*(double)*(kSpaceR+index) + 
		   		                      (double)*(kSpaceI+index)*(double)*(kSpaceI+index) );
          
            fwrite(&savedata,sizeof(double),1,pf);
			
		}
	} // K-Space Magnitude (end)

	fprintf_s(logfile,"%s\n", "K-Space Magnitude of the Signal Saved");

	fclose (pf);
	} // save data

	printf("%s\n", "FT Processing Completed");
    fprintf_s(logfile,"%s\n", "FT Processing Completed");

	fclose(logfile);

	// FIFO memory deallocation method
	free(kSpaceR);
	free(kSpaceI);
	free(Signal);

}

void OnInverseFourierTransform(char filename[], int rcxres, int rcyres)
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;
	
	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;
	
	double phase;

	//2010
	double emittingSource = 0.98; 
	double scale = ((double)rcxres*rcyres*emittingSource); 
	//2010

	FILE * logfile;
	char logfilename[128]="INV-FourierT.log";

	FILE *image;
	char imageFilename[256];

	double * kSpaceR = 0;
	double * kSpaceI = 0;
	double * reconSignal = 0;

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	 exit(0);
	
	} else { // allocate memory

		
	printf("%s\n", "Now INV FT Processing...");
    fprintf(logfile,"%s\n", "Now INV FT Processing...");

	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}


	if ((reconSignal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
	
		fprintf(logfile,"%s\n", "Not enough memory to allocate Imaginary Image data: Exit");
	
		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);

		exit(0);

	}

	} // allocate memory

	
	//// read image data and initialize pointers
    sprintf(imageFilename, "%s%s", "K-SpaceR-", filename);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceR+index) = (double) number;

		
			}

		}

		fclose(image);

	}// read data and initialize pointers


    char imageFilename2[128];

	sprintf(imageFilename2, "%s%s", "K-SpaceI-", filename);


    if ((image = fopen(imageFilename2,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceI+index) = (double) number;

			}

		}

		fclose(image);


		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				*(reconSignal+index) = (double)0.0;
					
			}

		}


	}// read data and initialize pointers

	double real = 0.0, imaginary = 0.0;
	
	///// INV Fourier Transform //////
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{
		
	    	dx = ((int) i - NofYpixels/2);
		    dy = ((int) j - NofXpixels/2);
		
	  	    k2 = ((int)(dx*NofXpixels)+dy);

			w = ((j*NofYpixels)+i);

			real = (double)0.0;
			imaginary = (double)0.0;

			
			for (int s=0; s<NofYpixels; s++)
			{ ///process k-space data

				for (int p=0; p<NofXpixels; p++)
				{ 

					ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);

					k3 = ((int)(dp*NofYpixels)+ds);  
				
					t = ((p*NofYpixels)+s);
					
					phase = ((double) (2.0 * pi * k2 * k3) / (NofXpixels*NofYpixels));

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (begin)**/
					real += ((double) *(kSpaceR+t) * cos( (double) phase)) + ((double) *(kSpaceI+t) * (double) sin((double)phase));

					imaginary += -((double) *(kSpaceR+t) * sin((double)phase)) + ((double) *(kSpaceI+t) * cos((double)phase)); 
					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (end)**/
			}

		}///process k-space data 

			*(reconSignal+w) =  (double) sqrt( ((double) real * real)  + ((double) imaginary * imaginary) );

			*(reconSignal+w) /= (double)scale;
		}
	} ///process k-space data



	// scale *(reconSignal) (begin) 
	double max=-*(reconSignal);
	double min=*(reconSignal);

	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{

			w = ((j*NofYpixels)+i);

				if( (double)*(reconSignal+w) > (double)max ) 
			
					max = (double)*(reconSignal+w);
              
				if( (double)*(reconSignal+w) < (double)min ) 
			
					min = (double)*(reconSignal+w);
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{

				if ( max == min ) (double)*(reconSignal+w);

				else (double)*(reconSignal+w)  =  (double)min - ( ( (double)*(reconSignal+w) * (min - max) ) / (double) ft_SCALE);
			
		} // y dim
        
	}  // x dim
	// scale *(reconSignal+w) (end)

	double savedata = 0.0;
	FILE * pf;
	char reconFilename[128];

	sprintf(reconFilename, "%s%s", "reconSignal-", filename);


    fprintf(logfile, "%s\t%s\n", "Now Saving Reconstructed Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(reconSignal);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(reconSignal+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "Reconstructed Signal Saved");

	fclose (pf);
	} // save data


    printf("%s\n", "Inverse FT Processing Completed");
    fprintf(logfile,"%s\n", "Inverse FT Processing Completed");

	fclose(logfile);
		
	// FIFO memory deallocation method
	free(kSpaceR);
	free(kSpaceI);
	free(reconSignal);

}

void OnInverseFourierTransform(char filename1[], char filename2[], char filename3[], char filename4[], int rcxres, int rcyres, char filename5[])
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;
	
	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;
	
	double phase;

	
	//2010
	double emittingSource = 1.0; 
	double scale = ((double)rcxres*rcyres*emittingSource); 
	//2010
	

	FILE *image;
	char imageFilename[256];

	double * kSpaceR = 0;
	double * kSpaceI = 0;

	double * kSpaceR_Subtract = 0;
	double * kSpaceI_Subtract = 0;
	double * reconSignal = 0;

    // allocate memory
	printf("%s\n", "Now INV FT Processing...");

	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		std::cout << "Not enough memory to allocate Real Image data: Exit" << endl;
   
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		std::cout << "Not enough memory to allocate Imaginary Image data: Exit" << endl;
   
		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}

	if ((kSpaceR_Subtract = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		std::cout << "Not enough memory to allocate Image data: Exit" << endl;
        free(kSpaceR);
		free(kSpaceI);
		exit(0);

	}

	if ((kSpaceI_Subtract = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		std::cout << "Not enough memory to allocate Real Image data: Exit" << endl;
   
		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);
		free(kSpaceR_Subtract);
		exit(0);

	}

	if ((reconSignal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
	
		std::cout << "Not enough memory to allocate Imaginary Image data: Exit" << endl;
	
		// FIFO memory deallocation method
	    free(kSpaceR);
		free(kSpaceI);
		free(kSpaceR_Subtract);
		free(kSpaceI_Subtract);
		exit(0);

	}

	
	//// read image data and initialize pointers (begins)
    sprintf(imageFilename, "%s", filename1);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 std::cout << "Cannot open Image File: " << imageFilename << endl;

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceR+index) = (double) number;

		
			}

		}

		fclose(image);

	}// read data and initialize pointers


    char imageFilename2[128];

	sprintf(imageFilename2, "%s", filename2);


    if ((image = fopen(imageFilename2,"rb+"))==NULL)
	{
	
	 std::cout << "Cannot open Image File: " << imageFilename2 << endl;

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceR_Subtract+index) = (double) number;

			}

		}

		fclose(image);


		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				*(reconSignal+index) = (double)0.0;
					
			}

		}


	}// read data and initialize pointers

	char imageFilename3[128];

	sprintf(imageFilename3, "%s", filename3);

    if ((image = fopen(imageFilename3,"rb+"))==NULL)
	{
	
	 std::cout << "Cannot open Image File: " << imageFilename << endl;
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceI+index) = (double) number;

		
			}

		}

		fclose(image);

	}

	char imageFilename4[128];

	sprintf(imageFilename4, "%s", filename4);

    if ((image = fopen(imageFilename4,"rb+"))==NULL)
	{
	
	 std::cout << "Cannot open Image File: " << imageFilename << endl;
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceI_Subtract+index) = (double) number;

		
			}

		}

		fclose(image);

	} //// read image data and initialize pointers (ends)

	double real = 0.0, imaginary = 0.0;
	
	///// INV Fourier Transform //////
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{
		
	    	dx = ((int) i - NofYpixels/2);
		    dy = ((int) j - NofXpixels/2);
		
	  	    k2 = ((int)(dx*NofXpixels)+dy);

			w = ((j*NofYpixels)+i);

			real = (double)0.0;
			imaginary = (double)0.0;

			
			for (int s=0; s<NofYpixels; s++)
			{ ///process k-space data

				for (int p=0; p<NofXpixels; p++)
				{ 

					ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);

					k3 = ((int)(dp*NofYpixels)+ds);  
				
					t = ((p*NofYpixels)+s);
					
					phase = ((double) (2.0 * pi * k2 * k3) / (NofXpixels*NofYpixels));

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (begin)**/
					real += ((double) ((double)*(kSpaceR+t) - *(kSpaceR_Subtract+t) ) * (double) cos( (double) phase)) + 
						    ((double) ((double)*(kSpaceI+t) - *(kSpaceI_Subtract+t) ) * (double) sin((double) phase));

					imaginary += -((double) ((double)*(kSpaceR+t) - *(kSpaceR_Subtract+t) ) * sin((double)phase)) + 
						          ((double) ((double)*(kSpaceI+t) - *(kSpaceI_Subtract+t) ) * cos((double)phase));  
					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (end)**/

			}

		}///process k-space data 

			*(reconSignal+w) =  (double) sqrt( ((double) real * real)  + ((double) imaginary * imaginary) );

	     	*(reconSignal+w) /= (double)scale;

		}
	} ///process k-space data

	double savedata = 0.0;
	FILE * pf;
	char reconFilename[128];

	sprintf(reconFilename, "%s%s%s", "RecontructedSignal-", filename5 , ".img");

	std::cout << "Now Saving Reconstructed Signal in File: " << reconFilename << endl;

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 std::cout << "Cannot open file to save K-Space Signal" << endl;

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(reconSignal+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	std::cout << "Reconstructed Signal Saved" << endl;

	fclose (pf);
	} // save data
 
    std:: cout << "Inverse FT Processing Completed" << endl;
	
	// FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);

}