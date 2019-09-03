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

//  Project Title: Intensity-Curvature Functional Based High Pass Filters: HPF (traditional), ICF, & Swarm Particle Optimization
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
void OnInverseFourierTransform(char filename1[], char filename2[], char filename3[], char filename4[], int rcxres, int rcyres, char filename5 []);

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

		double **ICF; // pointer to the matrix entry

		double **Y; // pointer to the matrix entry

		double **X; // pointer to the matrix entry

		double **TF; // pointer to the matrix entry

		double **Xt; // pointer to the matrix entry

		double **theta_x; // pointer to the matrix entry

		double **theta_y; // pointer to the matrix entry

		double **omega_f; // pointer to the matrix entry

		double **Filter2DH; // pointer to the matrix entry

		double **SwarmH; // pointer to the matrix entry

		double **SwarmHPF; // pointer to the particle velocity

		double **particleBestPosition; // pointer to the particle best position

	    double **ParticlePosition; // pointer to the particle position

		double **g; // swarm's best known position

		double **store_particleBestPosition; // storage 

		double **store_g; // storage

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

	 pointer->ICF = new double*[this->n1];

	 pointer->Y = new double*[this->n1];

	 pointer->X = new double*[this->n1];

	 pointer->TF = new double*[this->n1];

	 pointer->Xt = new double*[this->n1];

	 pointer->theta_x = new double*[this->n1];

	 pointer->theta_y = new double*[this->n1];

	 pointer->omega_f = new double*[this->n1];

	 pointer->Filter2DH = new double*[this->n1]; 

	 pointer->SwarmH = new double*[this->n1]; 

	 pointer->SwarmHPF = new double*[this->n1];

	 pointer->particleBestPosition = new double*[this->n1];

	 pointer->ParticlePosition = new double*[this->n1];
	
	 pointer->g = new double*[this->n1];

	 pointer->store_particleBestPosition = new double*[this->n1];

	 pointer->store_g = new double*[this->n1];


	 for( int v=0; v < this->n1; v++ ) { // (1)
		 
		 pointer->Signal[v] = new double[this->n2];

		 pointer->ICF[v] = new double[this->n2];

		 pointer->Y[v] = new double[this->n2];

		 pointer->X[v] = new double[this->n2];

		 pointer->TF[v] = new double[this->n2];

		 pointer->Xt[v] = new double[this->n2];

		 pointer->theta_x[v] = new double[this->n2];

		 pointer->theta_y[v] = new double[this->n2];

		 pointer->omega_f[v] = new double[this->n2];

		 pointer->Filter2DH[v] = new double[this->n2];

		 pointer->SwarmH[v] = new double[this->n2];

		 pointer->SwarmHPF[v] = new double[this->n2];

		 pointer->particleBestPosition[v] = new double[this->n2];

		 pointer->ParticlePosition[v] = new double[this->n2];

		 pointer->g[v] = new double[this->n2];

		 pointer->store_particleBestPosition[v] = new double[this->n2];

		 pointer->store_g[v] = new double[this->n2];


	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for(int v=0; v < this->n1; v++ ) { // (a)

			for( int f=0; f < this->n2 ; f++ ) { // (b)
		 
			pointer->Signal[v][f] = (double)0.0;

			pointer->ICF[v][f] = (double)0.0;

			pointer->Y[v][f] = (double)0.0;

			pointer->X[v][f] = (double)0.0;

			pointer->TF[v][f] = (double)0.0;

			pointer->Xt[v][f] = (double)0.0;

			pointer->theta_x[v][f] = (double)0.0;

			pointer->theta_y[v][f] = (double)0.0;

			pointer->omega_f[v][f] = (double)0.0;

			pointer->Filter2DH[v][f] = (double)0.0;

			pointer->SwarmH[v][f] = (double)0.0;
			
			pointer->SwarmHPF[v][f] = (double) 0.0;

			pointer->particleBestPosition[v][f] = (double)0.0;
	
			pointer->ParticlePosition[v][f] = (double) 0.0;

			pointer->g[v][f] = (double) 0.0;

			pointer->store_particleBestPosition[v][f] = (double) 0.0;

			pointer->store_g[v][f] = (double) 0.0;


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
		

	sprintf(outputFile, "%s","ICF.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ICF[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","Y.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->Y[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","X.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->X[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","Xt.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->Xt[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)
	
	sprintf(outputFile, "%s","TF.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->TF[v][f],sizeof(double),1,savedata);

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

	sprintf(outputFile, "%s","SwarmH.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->SwarmH[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	// SwarmHPF.img is the high pass filtered signal:
	sprintf(outputFile, "%s","SwarmHPF.img"); 
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->SwarmHPF[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","particleBestPosition.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->particleBestPosition[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ParticlePosition.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ParticlePosition[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)
	
	sprintf(outputFile, "%s","G.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->g[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

} // saveImages



int main ( int argc, char * argv[] ) {

	char outputFile[128]="SRE2D2018.log";

	FILE * savedata;

	double MAX = 5000000000000000000.0;

if (argc < 14) { std::cout << endl;
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
				 std::cout << "Please enter the tolerance (double) > 0 and < 1" << endl;
				 std::cout << "Please enter the value of omega (double) > 0 and < 1" << endl;
				 std::cout << "Please enter the value of Phi_P (double) > 0 and < 1" << endl;
				 std::cout << "Please enter the value of Phi_G (double) > 0 and < 1" << endl;
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

	double tolerance = atof(argv[10]);
	double omegaSwarm = atof(argv[11]);
	double Phi_P = atof(argv[12]);
	double Phi_G = atof(argv[13]);
	double error;

	std::cout << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
	std::cout << "The pixel size along the X direction is: " << atof(argv[4]) << endl;
	std::cout << "The pixel size along the Y direction is: " << atof(argv[5]) << endl;
	std::cout << "The XY rotation angle is: " << atof(argv[8]) << endl;
	std::cout << "The cutOFF frequency of the traditional HPF is: " << atof(argv[9]) << endl;
	std::cout << "The value of the tolerance is: " << atof(argv[10]) << endl;
	std::cout << "The value of omega is: " << atof(argv[11]) << endl;
	std::cout << "The value of Phi_P is: " << atof(argv[12]) << endl;
	std::cout << "The value of Phi_G is: " << atof(argv[13]) << endl;

	
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The XY rotation angle is: ", theta);
	fprintf(savedata,"%s%lf\n", "The cutOFF frequency of the traditional HPF is: ", cutOFF_frequency);
	fprintf(savedata,"%s%lf\n", "The value of the tolerance is: ", atof(argv[10]));
	fprintf(savedata,"%s%lf\n", "The value of omega is: ", atof(argv[11]));
	fprintf(savedata,"%s%lf\n", "The value of Phi_P is: ", atof(argv[12]));
	fprintf(savedata,"%s%lf\n", "The value of Phi_G is: ", atof(argv[13]));
	

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

           else SRE.pointer->theta_x[i1][i2] = (double) ft_SCALE * (min - SRE.pointer->theta_x[i1][i2]) / (min - max) ;
			
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

           else SRE.pointer->theta_y[i1][i2] = (double) ft_SCALE * (min - SRE.pointer->theta_y[i1][i2]) / (min - max) ;
			
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

           else SRE.pointer->omega_f[i1][i2] = (double) ft_SCALE * (min - SRE.pointer->omega_f[i1][i2]) / (min - max) ;
			
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

           else SRE.pointer->Filter2DH[i1][i2] = (double) ft_SCALE * (min - SRE.pointer->Filter2DH[i1][i2]) / (min - max) ;
			
		} // y dim
        
	}  // x dim 
	
	///  scale data (end)
	
	 std::cout << "High Pass Filter calculated" << endl;
	// build the FILTER function www.wikipedia.org & www.dspguide.com/ch14/5.htm (end)

	// calculate ICF (begin)
	std::cout << "Compute Intensity-Curvature Functional" << endl;
	
	for (int i1=0; i1 < n1-1; i1++) {// x dim

	    for (int i2=0; i2 < n2-1; i2++) { // y dim


			SRE.pointer->X[i1][i2] = (double)SRE.pointer->Signal[i1][i2];
	 
			SRE.pointer->Y[i1][i2] = (double)SRE.pointer->Signal[i1][i2];


		} // y dim
	} // x dim
	
	double k, s;
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
       
		   if ( (double)s == 0.0 && (double)k == 0.0 ) SRE.pointer->ICF[i1][i2] = (double)1.0; // de l'Hopital
           
		   else if ( (double)s == 0.0 && (double)k != 0.0 ) SRE.pointer->ICF[i1][i2] = (double)0.0;
		   
		   else if ( (double)s != 0.0 && (double)k == 0.0 ) SRE.pointer->ICF[i1][i2] = (double)0.0;

		   else  if ( (double)s != 0.0 && (double)k != 0.0 ) SRE.pointer->ICF[i1][i2] = (double)s/k;

		   /// calculation of the the input function: X, 
		   /// the output function: Y, and the transfer function: TF (begins)
		   double s_Y = ( (double) misplacement_X  * misplacement_Y  * SRE.pointer->X[i1][i2] );

		   double k_Y =  ( (double) misplacement_X  * misplacement_Y  * SRE.pointer->X[i1][i2] +
			               (double) misplacement_X  * misplacement_X  * misplacement_Y  * 
			              ((double) (SRE.pointer->Signal[i1+1][i2] - SRE.pointer->X[i1][i2]) / 2.0) +
                           (double) misplacement_X  * misplacement_Y  * misplacement_Y  * 
						  ((double) (SRE.pointer->Signal[i1][i2+1] - SRE.pointer->X[i1][i2]) / 2.0) +
			             ( (double) misplacement_X  * misplacement_Y  * misplacement_X  * misplacement_Y  *
	                     ( (double)SRE.pointer->X[i1][i2] + SRE.pointer->Signal[i1+1][i2+1] - 
						           SRE.pointer->Signal[i1+1][i2] - SRE.pointer->Signal[i1][i2+1]) / 4.0) );

		   if ( (double)s_Y == 0.0 && (double)k_Y == 0.0 ) SRE.pointer->Y[i1][i2] = (double)1.0; // de l'Hopital
           
		   else if ( (double)s_Y == 0.0 && (double)k_Y != 0.0 ) SRE.pointer->Y[i1][i2] = (double)0.0;
		   
		   else if ( (double)s_Y != 0.0 && (double)k_Y == 0.0 ) SRE.pointer->Y[i1][i2] = (double)0.0;

		   else  if ( (double)s_Y != 0.0 && (double)k_Y != 0.0 ) SRE.pointer->Y[i1][i2] = (double)s_Y / (double) k_Y;

		    double k_X = (  (double) 1.0 - ( (double) SRE.pointer->Y[i1][i2] *
			                          ( (double) 1.0 - misplacement_X / 2.0 - (double) misplacement_Y / 2.0 
									  + (double) misplacement_X  * misplacement_Y / 4.0 ) ) );

			double s_X = ( ((double) misplacement_X  * SRE.pointer->Signal[i1+1][i2] / 2.0) +
						   ((double) misplacement_Y  * SRE.pointer->Signal[i1][i2+1] / 2.0) + 
							(double) misplacement_X  * misplacement_Y * 
						   ((double) SRE.pointer->Signal[i1+1][i2+1] - SRE.pointer->Signal[i1+1][i2] - SRE.pointer->Signal[i1][i2+1] / 4.0) ) *
						    (double) SRE.pointer->Y[i1][i2];
			
			if ( (double)s_X == 0.0 && (double)k_X == 0.0 ) SRE.pointer->X[i1][i2] = (double)1.0; // de l'Hopital
           
		    else if ( (double)s_X == 0.0 && (double)k_X != 0.0 ) SRE.pointer->X[i1][i2] = (double)0.0;
		   
		    else if ( (double)s_X != 0.0 && (double)k_X == 0.0 ) SRE.pointer->X[i1][i2] = (double)0.0;

		    else  if ( (double)s_X != 0.0 && (double)k_X != 0.0 ) SRE.pointer->X[i1][i2] = (double)s_X / (double) k_X;	   	 

			  SRE.pointer->Xt[i1][i2] = (double) SRE.pointer->Signal[i1][i2];

			  if ( (double)SRE.pointer->Y[i1][i2] == 0.0 && (double)SRE.pointer->Xt[i1][i2] == 0.0 ) 
				  
				  SRE.pointer->TF[i1][i2] = (double)1.0; // de l'Hopital
           
			  else if ( (double)SRE.pointer->Y[i1][i2] == 0.0 && (double)SRE.pointer->Xt[i1][i2] != 0.0 ) 
				  
				  SRE.pointer->TF[i1][i2] = (double)0.0;
		   
		      else if ( (double)SRE.pointer->Y[i1][i2] != 0.0 && (double)SRE.pointer->Xt[i1][i2] == 0.0 ) 
				  
				  SRE.pointer->TF[i1][i2] = (double)0.0;

		      else if ( (double)SRE.pointer->Y[i1][i2] != 0.0 && (double)SRE.pointer->Xt[i1][i2] != 0.0 ) 
				  
				  SRE.pointer->TF[i1][i2] = (double)SRE.pointer->Y[i1][i2] / (double) SRE.pointer->Xt[i1][i2];
		 /// calculation of the the input function: X, 
		/// the output function: Y, and the transfer function: TF (ends)
		  

		} // y dim
        
	}  // x dim


	std::cout << "Intensity-Curvature Functional Calculated" << endl;
	// calculate ICF(end)

	// scale ICF (begin) 
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->ICF[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->ICF[i1][i2];
              
		if( SRE.pointer->ICF[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->ICF[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->ICF[i1][i2] = (double)0.0;

           else SRE.pointer->ICF[i1][i2] = (double) ft_SCALE * (min - SRE.pointer->ICF[i1][i2]) / (min - max) ;
			
		} // y dim
        
	}  // x dim
	// scale ICF (end)

	double RC;
	double condition_fpi_fg;
	
	// Particle swarm optimization (begin)
	// www.en.wikipedia.org/wiki/Particle_swarm_optimization (begins)

	// initialize (begins)
	// Initialize the particle's position with the signal + random number:
	for (int i1=0; i1 < n1; i1++) { // x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
		
			SRE.pointer->ParticlePosition[i1][i2] = (double)SRE.pointer->Signal[i1][i2] + ((double)rand() / (double) RAND_MAX); // xi
			

			} // for each dimension (y)

		} // for each dimension (x)

	    // Initialize the particle's best known position to its initial position (the signal + random number):
    	for (int i1=0; i1 < n1; i1++) {// x dim
		   	
			for (int i2=0; i2 < n2; i2++) { // y dim
		
				SRE.pointer->particleBestPosition[i1][i2] = (double) SRE.pointer->ParticlePosition[i1][i2]; // pi <- xi
	
			   // initialize 'g': the swarm's best known position (the signal + random number):
			   SRE.pointer->g[i1][i2] = (double) SRE.pointer->Signal[i1][i2] + ((double) rand() / (double) RAND_MAX);

			  
			   // condition calculated for if statement (1)
			   condition_fpi_fg = (double) SRE.pointer->particleBestPosition[i1][i2] - (double) SRE.pointer->g[i1][i2]; 
			   // expected to be negative
			
			   	if( condition_fpi_fg <= (double) 0.0 ) { // if (1)

			    // Update the swarm's best known position:
				SRE.pointer->g[i1][i2] = SRE.pointer->particleBestPosition[i1][i2]; // g <- pi

	         } // if (1)

			} // for each dimension (y)

		} // for each dimension (x)
		// initialize (ends)

			
    // Initialize the particle's velocity (SwarmHPF) with the signal plus a random number
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		 SRE.pointer->SwarmHPF[i1][i2] = (double) SRE.pointer->Signal[i1][i2] + ((double) rand() / (double) RAND_MAX); // vi
	
		} // for each dimension (y)

	} // for each dimension (x)


		int Iterations = 0;
		
		do { // do while loop begins here
		
			Iterations++;
			
		    // Update the particle's velocity:
		    for (int i1=1; i1 < n1; i1++) {// x dim
       	
				for (int i2=1; i2 < n2; i2++) { // y dim

					//Pick random numbers:
					double random_p = ((double) rand() / (double) RAND_MAX); // rp
					double random_g = ((double) rand() / (double) RAND_MAX); // rg
					
					RC = (double) 1.0 / ((double) 2.0 * pi * cutOFF_frequency);

				    highPassFilter = ((double) RC) / ((double) RC + deltaT);

					SRE.pointer->SwarmH[i1][i2] = 
						
				    ((double)highPassFilter * SRE.pointer->SwarmH[i1-1][i2]) + 
				   (((double) Phi_P * random_p) * 
				    ((double)SRE.pointer->particleBestPosition[i1][i2] - (double)SRE.pointer->ParticlePosition[i1-1][i2])) +
				   (((double) Phi_G * random_g) * 
				   ((double)SRE.pointer->g[i1][i2] - (double)SRE.pointer->ParticlePosition[i1-1][i2])); // vi <-
					
					/* equation of the HPF 
					//SRE.pointer->SwarmH[i1][i2] = ((double)highPassFilter * y ) + 
				                                   (((double)highPassFilter) * ((double)x - w )); */
					
					//Update the particle's position:
					SRE.pointer->ParticlePosition[i1][i2] += (double) SRE.pointer->SwarmH[i1][i2]; // xi <- xi + vi
														   	
					// expected to be negative:
					double condition_fxi_fpi = (double)SRE.pointer->ParticlePosition[i1][i2] - 
						                       (double)SRE.pointer->particleBestPosition[i1][i2];

					// expected to be negative:
					double condition_fpi_fg  = (double)SRE.pointer->particleBestPosition[i1][i2] - 
						                       (double)SRE.pointer->g[i1][i2];
 
					if( condition_fxi_fpi <= (double)0.0 ) { // if (2)

						// Update the particle's best known position using the Delta Rule (first order partial derivative of the cost function (**)
						SRE.pointer->particleBestPosition[i1][i2] = SRE.pointer->store_particleBestPosition[i1][i2] -
							                                          (double) 2.0 * ((double) SRE.pointer->SwarmH[i1][i2] - 
						                                              (double) SRE.pointer->SwarmHPF[i1][i2]) *
																     ((double) Phi_P * random_p);
					}

					if( condition_fpi_fg <= (double)0.0) { // if (3)

						// Update the Swarm's best known position using the Delta Rule (first order partial derivative of the cost function (**)   
						SRE.pointer->g[i1][i2] = SRE.pointer->store_g[i1][i2] - 
							                       (double) 2.0 * ((double) SRE.pointer->SwarmH[i1][i2] - 
						                           (double) SRE.pointer->SwarmHPF[i1][i2]) *
									              ((double) Phi_G * random_g);
					}

				} // for each dimension (y)

			} // for each dimension (x)

				
					// calculate the error of fit (cost function (**) (begins)
					error = 0.0;
					for (int i1=0; i1 < n1; i1++) {// x dim
       	
						for (int i2=0; i2 < n2; i2++) { // y dim

					error += pow( ((double) SRE.pointer->SwarmH[i1][i2] - 
						           (double) SRE.pointer->SwarmHPF[i1][i2]), 2.0); 
					
									} // for each dimension (y)

					} // for each dimension (x)
					// calculate the error of fit (cost function (**) (ends)

				    // print the error of fit: (begins)
					std::cout << "Error of fit: " << error << endl;
					fprintf(savedata,"%s%lf\n", "The error of fit is: ", error);
					if( error <= tolerance ) {	
				    std::cout << "The error of fit reached the desired tolerance: " << error << endl;
					std::cout << "Tolerance: " << tolerance << endl;
					break; // (***)
					}
					// print the error of fit: (ends)
					
					// update the value of the particleVelocity (begins). 
				    // The particleVelocity (SwarmHPF) is not the same as SwarmH because the
					// break instruction in the above calculation of the Error 
					// of fit (see (***) forces the program to exit the while loop before
					// the update here reported: 
					for (int i1=0; i1 < n1; i1++) { // x dim
       	
						for (int i2=0; i2 < n2; i2++) { // y dim
						
							SRE.pointer->SwarmHPF[i1][i2] = (double)SRE.pointer->SwarmH[i1][i2];

							} // for each dimension (y)

					} // for each dimension (x)
					// update the value of particleVelocity (ends)

					// store particle Best Position and store g; of the current iteration so that can be used at the next iteration
			    	for (int i1=0; i1 < n1; i1++) {// x dim
       	
							for (int i2=0; i2 < n2; i2++) { // y dim
					
							SRE.pointer->store_particleBestPosition[i1][i2] = (double)SRE.pointer->particleBestPosition[i1][i2];
							
							SRE.pointer->store_g[i1][i2] = (double)SRE.pointer->g[i1][i2];

							} // for each dimension (y)

					} // for each dimension (x)
		
		} while (error > tolerance && Iterations < 100); // do while loop ends here
	   // www.en.wikipedia.org/wiki/Particle_swarm_optimization (ends)
	   // Particle swarm optimization (ends)


	std::cout << "Swarm Filter calculated" << endl;


	// scale Swarm (begin) 
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRE.pointer->SwarmHPF[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->SwarmHPF[i1][i2];
              
		if( SRE.pointer->SwarmHPF[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->SwarmHPF[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) SRE.pointer->SwarmHPF[i1][i2] = (double)0.0;

           else SRE.pointer->SwarmHPF[i1][i2] = (double) ft_SCALE * (min - SRE.pointer->SwarmHPF[i1][i2]) / (min - max) ;
			
		} // y dim
        
	}  // x dim
	// scale ICF (end)

	SRE.save(); // save all of the images

	
	std::cout << "Please Wait...Now running the Fourier Transformations on the Signal" << endl;

	OnFourierTransform(imageFileName, n2, n1);

	OnInverseFourierTransform(imageFileName, n2, n1);

	
	std::cout << "Please Wait...Now running the Fourier Transformations on theta_x" << endl;

	OnFourierTransform("theta_x.img", n2, n1);

	OnInverseFourierTransform("theta_x.img", n2, n1);

	
	std::cout << "Please Wait...Now running the Fourier Transformations on theta_y" << endl;

	OnFourierTransform("theta_y.img", n2, n1);

	OnInverseFourierTransform("theta_y.img", n2, n1);


	std::cout << "Please Wait...Now running the Fourier Transformations on omega_f" << endl;

	OnFourierTransform("omega_f.img", n2, n1);

	OnInverseFourierTransform("omega_f.img", n2, n1);

	
	std::cout << "Please Wait...Now running the Fourier Transformations on the HP Filter" << endl;

	OnFourierTransform("Filter2DH.img", n2, n1);

	OnInverseFourierTransform("Filter2DH.img", n2, n1);

		
	std::cout << "Please Wait...Now running the Fourier Transformations on the ICF" << endl;

	OnFourierTransform("ICF.img", n2, n1);

	OnInverseFourierTransform("ICF.img", n2, n1);

	
	std::cout << "Please Wait...Now running the Fourier Transformations on the Swarm" << endl;

	OnFourierTransform("SwarmHPF.img", n2, n1);

	OnInverseFourierTransform("SwarmHPF.img", n2, n1);


	std::cout << "Please Wait...Now running the inverse Fourier Transformation procedures" << endl;


		OnInverseFourierTransform("K-SpaceR-Filter2DH.img",
			                      "K-SpaceR-ICF.img",
			                      "K-SpaceI-Filter2DH.img",
			                      "K-SpaceI-ICF.img", n2, n1, "Filter2DH");

		OnInverseFourierTransform("K-SpaceR-SwarmHPF.img",
			                      "K-SpaceR-ICF.img",
			                      "K-SpaceI-SwarmHPF.img",
			                      "K-SpaceI-ICF.img", n2, n1, "SwarmHPF");

		OnInverseFourierTransform("K-SpaceR-theta_x.img",
			                      "K-SpaceR-ICF.img",
			                      "K-SpaceI-theta_x.img",
			                      "K-SpaceI-ICF.img", n2, n1, "theta_x");

		OnInverseFourierTransform("K-SpaceR-theta_y.img",
			                      "K-SpaceR-ICF.img",
			                      "K-SpaceI-theta_y.img",
			                      "K-SpaceI-ICF.img", n2, n1, "theta_y");

		OnInverseFourierTransform("K-SpaceR-omega_f.img",
			                      "K-SpaceR-ICF.img",
			                      "K-SpaceI-omega_f.img",
			                      "K-SpaceI-ICF.img", n2, n1, "omega_f");


	std::cout << "Please Wait...Now running the Fourier Transformations on RecontructedSignal-Filter2DH" << endl;

	OnFourierTransform("RecontructedSignal-Filter2DH.img", n2, n1);

	OnInverseFourierTransform("RecontructedSignal-Filter2DH.img", n2, n1);

	
	std::cout << "Please Wait...Now running the Fourier Transformations on RecontructedSignal-SwarmHPF" << endl;

	OnFourierTransform("RecontructedSignal-SwarmHPF.img", n2, n1);

	OnInverseFourierTransform("RecontructedSignal-SwarmHPF.img", n2, n1);

		
	std::cout << "Please Wait...Now running the Fourier Transformations on RecontructedSignal-theta_x" << endl;

	OnFourierTransform("RecontructedSignal-theta_x.img", n2, n1);

	OnInverseFourierTransform("RecontructedSignal-theta_x.img", n2, n1);


	std::cout << "Please Wait...Now running the Fourier Transformations on RecontructedSignal-theta_y" << endl;

	OnFourierTransform("RecontructedSignal-theta_y.img", n2, n1);

	OnInverseFourierTransform("RecontructedSignal-theta_y.img", n2, n1);


	std::cout << "Please Wait...Now running the Fourier Transformations on RecontructedSignal-omega_f" << endl;

	OnFourierTransform("RecontructedSignal-omega_f.img", n2, n1);

	OnInverseFourierTransform("RecontructedSignal-omega_f.img", n2, n1);


	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	fclose(savedata);
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

void OnInverseFourierTransform(char filename1[], char filename2[], char filename3[], char filename4[], int rcxres, int rcyres, char filename5 [])
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