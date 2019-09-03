//  This file contains sourcecode distributed as freeware. 
//  The intellectual property of the sourcecode is shown 
//  here to belong to Carlo Ciulla.

// Reference (book):

// [2] Carlo Ciulla "Signal Resilient to Interpolation: An Exploration on the 
// Approximation Properties of the Mathematical Functions" - CreateSpace Publisher - 
// June 2013; pp. 1 - 336, ISBN: 978-1477567487.

// Disclaimer:

// The website here named www.sourcecodewebsiteCarloCiulla.com 2013 [2] does not intend 
// to convey the meaning of profit making for what pertains to the content
// provided. --->>> Instead, when the content is downloaded, the user(s) are
// kindly invited to donate money to charity organizations involved in 
// helping people in need of food and water. <<<---


//  Project Title: Bivariate SRI Interpolation
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


#define NSCALE 0.1	
#define ESCALE 0.1
#define TH 0

using namespace std;


class SRI2D2013 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data {

		double **OriginalSignal; // pointer to the matrix entry 

		double **Signal; // pointer to the matrix entry
		
		double **ClassicCurvature; // pointer to the matrix entry 
	
		double **ResilientCurvature; // pointer to the matrix entry 

		double ** ThePixICResilient; // pointer to the matrix entry 

		double **shifted_OriginalSignal; // pointer to the matrix entry 
	
		double **SE; // pointer to the matrix entry 

		double **ThePixICClassic; // pointer to the matrix entry 

		double **the_SRI_shifted_OriginalSignal; // pointer to the matrix entry

		double **HybridSignal; // pointer to the matrix entry

		double **SE_SRI; // pointer to the matrix entry 

		double **ResilientSignal;	// pointer to the matrix entry

		double **ICF; // pointer to the matrix entry

		double **ictBI; // pointer to the matrix entry

		double **ictAI; // pointer to the matrix entry

	}*pointer; // pointer to the matrices

public:

	SRI2D2013(int x, int y) : n1(x), n2(y) { };// constructor 
	
	void allocateData();

	void save();

	~SRI2D2013() { } // destructor

};

void SRI2D2013::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;
			
	 pointer->OriginalSignal = new double*[this->n1];

	 pointer->Signal = new double*[this->n1];

	 pointer->ClassicCurvature = new double*[this->n1];

	 pointer->ResilientCurvature = new double*[this->n1];

	 pointer->ThePixICResilient = new double*[this->n1];

	 pointer->shifted_OriginalSignal = new double*[this->n1];

	 pointer->SE = new double*[this->n1];

	 pointer->ThePixICClassic = new double*[this->n1];

	 pointer->the_SRI_shifted_OriginalSignal = new double*[this->n1];

	 pointer->HybridSignal = new double*[this->n1];

	 pointer->SE_SRI = new double*[this->n1];

	 pointer->ResilientSignal = new double*[this->n1];

	 pointer->ICF = new double*[this->n1];

	 pointer->ictBI = new double*[this->n1];

	 pointer->ictAI = new double*[this->n1];


	 for( int v=0; v < this->n1; v++ ) { // (1)
		 
		 pointer->OriginalSignal[v] = new double[this->n2];

		 pointer->Signal[v] = new double[this->n2];

		 pointer->ClassicCurvature[v] = new double[this->n2];

		 pointer->ResilientCurvature[v] = new double[this->n2];

		 pointer->ThePixICResilient[v] = new double[this->n2];

		 pointer->shifted_OriginalSignal[v] = new double[this->n2];

		 pointer->SE[v] = new double[this->n2];

		 pointer->ThePixICClassic[v] = new double[this->n2];

		 pointer->the_SRI_shifted_OriginalSignal[v] = new double[this->n2];

  	     pointer->HybridSignal[v] = new double[this->n2];

		 pointer->SE_SRI[v] = new double[this->n2];

		 pointer->ResilientSignal[v] = new double[this->n2];

		 pointer->ICF[v] = new double[this->n2];

		 pointer->ictBI[v] = new double[this->n2];

		 pointer->ictAI[v] = new double[this->n2];


	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for( int v=0; v < this->n1; v++ ) { // (a)

			for( int f=0; f < this->n2 ; f++ ) { // (b)
		 
			pointer->OriginalSignal[v][f] = (double)0.0;

			pointer->Signal[v][f] = (double)0.0;

			pointer->ClassicCurvature[v][f] = (double)0.0;

			pointer->ResilientCurvature[v][f] = (double)0.0;

			pointer->ThePixICResilient[v][f] = (double)0.0;

			pointer->shifted_OriginalSignal[v][f] = (double)0.0;

			pointer->SE[v][f] = (double)0.0;

			pointer->ThePixICClassic[v][f] = (double)0.0;

			pointer->the_SRI_shifted_OriginalSignal[v][f] = (double)0.0;

			pointer->HybridSignal[v][f] = (double)0.0;

			pointer->SE_SRI[v][f] = (double)0.0;

			pointer->ResilientSignal[v][f] = (double)0.0;

			pointer->ICF[v][f] = (double)0.0;

			pointer->ictBI[v][f] = (double)0.0;

			pointer->ictAI[v][f] = (double)0.0;

			 } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


void SRI2D2013::save() { // saveImages

	FILE * savedata;
	char outputFile[128];
	
	sprintf(outputFile, "%s","OriginalSignal.img");

	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->OriginalSignal[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	
	sprintf(outputFile, "%s","ClassicCurvature.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ClassicCurvature[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)
	


	sprintf(outputFile, "%s","ResilientCurvature.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ResilientCurvature[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","ThePixICResilient.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ThePixICResilient[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	
	sprintf(outputFile, "%s","shifted_OriginalSignal.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->shifted_OriginalSignal[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)



	sprintf(outputFile, "%s","SE.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->SE[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	
	sprintf(outputFile, "%s","ThePixICClassic.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ThePixICClassic[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","the_SRI_shifted_OriginalSignal.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->the_SRI_shifted_OriginalSignal[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","SE_SRI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->SE_SRI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","ResilientSignal.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ResilientSignal[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","HybridSignal.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->HybridSignal[v][f],sizeof(double),1,savedata);

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

	sprintf(outputFile, "%s","ICTBI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ictBI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","ICTAI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ictAI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


} // saveImages



int main ( int argc, char * argv[] ) {

	char outputFile[128]="SRI2D.log";

	FILE * savedata;

	double MAX = 5000000000000000000.0;

if (argc < 13) { std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Before running the program, please make sure that the image is padded of 'n >= 0'";
				 std::cout << "number of pixels along X and Y" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the pixel size along the X direction (double)" << endl;
				 std::cout << "Please enter the pixel size along the Y direction (double)" << endl;
				 std::cout << "Please enter the misplacement along the X direction (double)" << endl;
				 std::cout << "Please enter the misplacement along the Y direction (double)" << endl;
				 std::cout << "Please enter the XY rotation angle (double)" << endl;
				 std::cout << "Please enter the a constant (double)" << endl;
				 std::cout << "Please enter the b constant (double)" << endl;
				 std::cout << "Please enter that 'n >= 0' number of pixels along X, Y and Z which will pad the image" << endl;
				 std::cout << "Processing will not be correct if you enter a value of n which is greater than" << endl;
				 std::cout << "the actual number of pixels along X and Y which was used to pad the image" << endl;
				 std::cout << "before running the program." << endl;
				 std::cout << "Please type 'n' to scale the Image Data or 'e' to exponentialize the Image Data" << endl;
				 std::cout << endl;
				 exit(0); }

else { // run the program (begin)

	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // processing (begin)

	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);

	double XPixelSize = atof(argv[4]);
	double YPixelSize = atof(argv[5]);

	double x_misplacement_X = atof(argv[6]);
	double y_misplacement_Y = atof(argv[7]);

	double theta = atof(argv[8]);

	double theAconstant = atof(argv[9]);

	double theBconstant = atof(argv[10]);

	char imageFileName[128];

	int PAD = atoi(argv[11]);

	unsigned char Fyn = 'n';

	Fyn = *argv[12];

	if ( PAD < 0 ) { // if test PAD
	
		std::cout << "Please make sure to pad the image with 'n >= 0' number of pixels along X, Y and Z" << endl;
	
		fprintf(savedata,"%s\n", "Please make sure to pad the image with 'n >= 0' number of pixels along X, Y and Z");

		fclose(savedata);
	
		exit(0);
	
	} // if test PAD

	sprintf(imageFileName, "%s", argv[1]);

	std::cout << endl;
	std::cout << "The image file name is: " << imageFileName << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
	std::cout << "The pixel size along the X direction is: " << atof(argv[4]) << endl;
	std::cout << "The pixel size along the Y direction is: " << atof(argv[5]) << endl;
	std::cout << "The XY rotation angle is: " << atof(argv[8]) << endl;
	std::cout << "The a constant is: " << atof(argv[9]) << endl;
	std::cout << "The b constant is: " << atof(argv[10]) << endl;
	std::cout << "The 'n' number of pixels along X and Y which you entered to pad the image is: " << atoi(argv[11]) << endl;
	std::cout << "Processing will not be correct if you enter a value of n which is greater than" << endl;
	std::cout << "the actual number of pixels along X and Y which was used to pad the image" << endl;
	std::cout << "before running the program." << endl;
	
	if ( tolower(Fyn) == 'e' )
	std::cout << "The program exponentialize the Image Data before processing: " << Fyn << endl;
	else if ( tolower(Fyn) == 'n' ) 
	std::cout << "The program scales the Image Data before processing: " << Fyn << endl;
	std::cout << endl;

	fprintf(savedata,"%s%s\n", "The image file name is: " , imageFileName);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The XY rotation angle is: ", theta);
	fprintf(savedata,"%s%lf\n", "The a constant is: ", theAconstant);
	fprintf(savedata,"%s%lf\n", "The b constant is: ", theBconstant);
	fprintf(savedata,"%s%d\n", "The 'n' number of pixels along X and Y which you entered to pad the image is: ", PAD);
	fprintf(savedata,"%s\n", "Processing will not be correct if you enter a value of n which is greater than");
	fprintf(savedata,"%s\n", "the actual number of pixels along X and Y which was used to pad the image");
	fprintf(savedata,"%s\n", "before running the program.");
	
	if ( tolower(Fyn) == 'e' )
	fprintf(savedata,"%s%c\n", "The program exponentialize the Image Data before processing: ", Fyn);
	else if ( tolower(Fyn) == 'n' )
	fprintf(savedata,"%s%c\n", "The program scales the Image Data before processing: ", Fyn);
	fprintf(savedata,"\n");

    double misplacement_X = ((double)1.0 - ( cos( (double)theta ) + sin( (double)theta ) ) + x_misplacement_X);
    double misplacement_Y = ((double)1.0 - ( -sin( (double)theta ) + cos( (double)theta ) ) + y_misplacement_Y);

      misplacement_X = ((double)misplacement_X/XPixelSize);
      misplacement_Y = ((double)misplacement_Y/YPixelSize);

	  //////////////////***********//////////////////////
	  // Above formula scales the misplacement to the  //
	  // pixel size the same way the following formula //
	  // would do: (min - misplacement)/(min - max)    //  
	  //////////////////***********//////////////////////

	SRI2D2013 SRI(n1,n2);

	SRI.allocateData();

	/// read image file (begin)
	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		SRI.~SRI2D2013();
		exit(0);

	} else { // else

	double number;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
			
		fread(&number,sizeof(double),1,pf);
		
		SRI.pointer->OriginalSignal[i1][i2] = (double)number;

		SRI.pointer->Signal[i1][i2] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl;

	double max=-MAX;
	double min=MAX;

	// make decision as to if standarzize or scale the Image data (begin)
	if ( tolower(Fyn) == 'e' ) {
		
	double Std = 0, Average = 0;
	long int counting = 0;
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) {

			Average += (double) SRI.pointer->OriginalSignal[i1][i2];

			counting++;

			}
			
		} // y dim
        
	}  // x dim 

	Average /= ((double)counting);

	counting = 0;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) {

			Std = (double)Std + ((double) (Average - SRI.pointer->OriginalSignal[i1][i2]) * 
				                          (Average - SRI.pointer->OriginalSignal[i1][i2]) );

			counting++;

			}

		} // y dim
        
	}  // x dim 

	Std = (double) sqrt( (double) Std / ((double)counting) );


	// standardize (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		   if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 

           SRI.pointer->Signal[i1][i2] = (double) ESCALE * ( (double) exp((double) (SRI.pointer->OriginalSignal[i1][i2] - Average) / Std) );
       
		   } else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) {
  	
			   SRI.pointer->Signal[i1][i2] = (double)0.0;
			
		   }

		} // y dim
        
	}  // x dim 
	// standardize (end)

	} else if ( tolower(Fyn) == 'n' )
	{ // make decision as to if standarzize or scale the Image data 

	/// compute max and min of data (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( SRI.pointer->OriginalSignal[i1][i2] > (double)max ) 
			
			max = (double)SRI.pointer->OriginalSignal[i1][i2];
              
		if( SRI.pointer->OriginalSignal[i1][i2] < (double)min ) 
			
			min = (double)SRI.pointer->OriginalSignal[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		   if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 

           if ( max == min ) SRI.pointer->Signal[i1][i2] = (double)0.0;

           else SRI.pointer->Signal[i1][i2] = (double) NSCALE * (min - SRI.pointer->OriginalSignal[i1][i2]) / (min - max) ;
       
		   } else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) {
  	
			   SRI.pointer->Signal[i1][i2] = (double)0.0;
			
		   }

		} // y dim
        
	}  // x dim 
	
	// scale (end)
	std::cout << "Image data scaled" << endl;
	} // make decision as to if standarzize or scale the Image data (end)

	/// run JAVA (begin)
	char command[200];
	sprintf(command, "%s%s%s%d%s%d%s%c", "java -jar make2DOmega2013.jar ", 
		                                 imageFileName, " ", n1, " ", n2, " ", Fyn); 
	system(command);
	/// run JAVA (end)

	
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
	
	SRI.~SRI2D2013();
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
	SRI.~SRI2D2013();
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
	SRI.~SRI2D2013();
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

	/// read omega (begin)
	FILE * readOmega;
	char inputFile[128];
	double readdata;
	
	sprintf(inputFile, "%s","omegaA.img");

	if ((readOmega = fopen(inputFile,"rb"))==NULL)
	{

		std::cout << "Cannot open input file, Now Exit..." << endl;
		fprintf(savedata,"%s\n", "Cannot open input file, Now Exit...");

		free(OMEGAA);
		free(OMEGAB);
		free(OMEGAAB);
		SRI.~SRI2D2013();
		exit(0);
		
	} else  { // (read)


	for( int i1=0; i1 < n1; i1++ ) { // (a)

		for( int i2=0; i2 < n2; i2++ ) {
	
		int index = ((int) i2*n1 + i1); 

		fread(&readdata, sizeof(double),1, readOmega);

		*(OMEGAA + index) = (double)readdata;

		}

	} // (a)

	fclose(readOmega);

	} // (read)

	sprintf(inputFile, "%s","omegaB.img");

	if ((readOmega = fopen(inputFile,"rb"))==NULL)
	{

		std::cout << "Cannot open input file, Now Exit..." << endl;
		fprintf(savedata,"%s\n", "Cannot open input file, Now Exit...");

	    free(OMEGAA);
		free(OMEGAB);
		free(OMEGAAB);
		SRI.~SRI2D2013();
		exit(0);

	} else  { // (read)


	for( int i1=0; i1 < n1; i1++ ) { // (a)

		for( int i2=0; i2 < n2; i2++ ) {
	
		int index = ((int) i2*n1 + i1); 

		fread(&readdata, sizeof(double),1, readOmega);

		*(OMEGAB + index) = (double)readdata;

		}

	} // (a)

	fclose(readOmega);

	} // (read)

	sprintf(inputFile, "%s","omegaAB.img");

	if ((readOmega = fopen(inputFile,"rb"))==NULL)
	{

		std::cout << "Cannot open input file, Now Exit..." << endl;
		fprintf(savedata,"%s\n", "Cannot open input file, Now Exit...");

	    free(OMEGAA);
		free(OMEGAB);
		free(OMEGAAB);
		SRI.~SRI2D2013();
		exit(0);

	} else  { // (read)


	for( int i1=0; i1 < n1; i1++ ) { // (a)

		for( int i2=0; i2 < n2; i2++ ) {
	
		int index = ((int) i2*n1 + i1); 

		fread(&readdata, sizeof(double),1, readOmega);

		*(OMEGAAB + index) = (double)readdata;

		}

	} // (a)

	fclose(readOmega);

	} // (read)
	/// read omega (end)
	/// allocate memory & store OMEGAs (end)

	std::cout << "Image data scaled" << endl;

    // shift the image of the given misplacement (begin)
	// calculate curvature of model interpolation function (begin)
	double omegaA, omegaB, omegaAB;
	double convA, convB, convAB;
	double TheShift, term;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 

		
		int index = ((int) i2*n1 + i1); 

		omegaA =  (double)*(OMEGAA + index);

		omegaB =  (double)*(OMEGAB + index);
		
		omegaAB = (double)*(OMEGAAB + index);


		convA = ((double)theAconstant * misplacement_X * misplacement_X) + ((double)misplacement_Y);

		convB = ((double)theBconstant * misplacement_Y * misplacement_Y) - ((double)misplacement_X);

		convAB = ((double)theAconstant * misplacement_X * misplacement_X * misplacement_Y) + 
			     ((double)theBconstant * misplacement_Y * misplacement_Y * misplacement_X);


		SRI.pointer->shifted_OriginalSignal[i1][i2] = 
			
		(double)SRI.pointer->OriginalSignal[i1][i2] + 
		
		((double)convA*omegaA + (double)convB*omegaB + (double)convAB*omegaAB);
		
		
		TheShift = (double)sqrt( ((double)misplacement_X * misplacement_X) + 
			                     ((double)misplacement_Y * misplacement_Y) );

		term = ((double)2.0 * theAconstant * misplacement_X + 2.0 * theBconstant * misplacement_Y);

		SRI.pointer->ClassicCurvature[i1][i2] = ((double)omegaA * 2.0 * theAconstant) +
		                                        ((double)omegaAB * 2.0 * theAconstant * misplacement_Y) + 
												((double)omegaB * 2.0 * theBconstant) +
												((double)omegaAB * 2.0 * theBconstant * misplacement_X) +
												((double)2.0 * omegaAB * term);

		} else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim


	std::cout << "Image shifted of the given misplacement" << endl;
	std::cout << "Curvature of model interpolation function calculated" << endl;
	// shift the image of the given misplacement (end)
	// calculate curvature of model interpolation function (end)



	// calculate the signal resilient to interpolation (begin)
	double Fxy, Hxy, Bxy, Cxy, Dxy;
	double termOne, termTwo, termThree, termFour, termFive;
	double b1, b2, b3, b4, b5;
	double c1, c2, c3, c4, c5;
	double d1, d2, d3, d4, d5;
	double x4, x3, x2, y4, y3, y2;
	double dBxydx, dBxydy;
	double e1, e2, e3, e4, e5;
	double f1, f2, f3, f4, f5;
	double dCxydx, dCxydy;
	double g1, g2, g3, g4, g5;
	double h1, h2, h3, h4, h5;
	double dDxydx, dDxydy;
	double j1, j2, j3, j4, j5;
	double k1, k2, k3, k4, k5;
	double dFxydx, dFxydy, d2Fxydxdy, d2Fxydydx;
	double L6, L7, L8, L9;
	double d2Fxydx2, d2Fxydy2;
	double L4, L5;
	double d2Bxydx2, d2Bxydy2;
	double m1, m2, m3, m4, m5;
	double o1, o2, o3, o4, o5;
	double d2Cxydx2, d2Cxydy2;
	double p1, p2, p3, p4, p5;
	double q1, q2, q3, q4, q5;
	double d2Dxydx2, d2Dxydy2;
	double r1, r2, r3, r4, r5;
	double s1, s2, s3, s4, s5;
	double d2Bxydxdy, d2Cxydxdy, d2Dxydxdy;
	double d2Bxydydx, d2Cxydydx, d2Dxydydx;
	double t1, t2, t3;
	double d2SRIdx2, d2SRIdy2, d2SRIdxdy, d2SRIdydx;
	double t4, t5, t6, t7, t8, t9; 
	double u1, u2, u3;


	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim


		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 


		int index = ((int) i2*n1 + i1); 

		omegaA =  (double)*(OMEGAA + index);

		omegaB =  (double)*(OMEGAB + index);
		
		omegaAB = (double)*(OMEGAAB + index);


		termOne =   ((double)2.0 * theAconstant * misplacement_X * misplacement_Y);
		termTwo =   ((double)theAconstant * misplacement_X * misplacement_Y * misplacement_Y);
		termThree = ((double)2.0 * theBconstant * misplacement_X * misplacement_Y);
		termFour =  ((double)theBconstant * misplacement_X * misplacement_X * misplacement_Y);
		termFive =  ((double)theAconstant * misplacement_X * misplacement_X * misplacement_Y) +
			        ((double)theBconstant * misplacement_X * misplacement_Y * misplacement_Y );

	    Hxy = ((double)omegaA*termOne + omegaAB*termTwo + omegaB*termThree + omegaAB*termFour + 2.0*omegaAB*termFive);

		Fxy = (((double)omegaA * 2.0 * theAconstant + omegaB * 2.0 * theBconstant) * 
			  ((double)misplacement_X * misplacement_Y)) - ((double)Hxy);

		x4 = (double) pow( (double)misplacement_X, (double)4.0 );
		x3 = (double) pow( (double)misplacement_X, (double)3.0 );
		x2 = (double) pow( (double)misplacement_X, (double)2.0 );

		y4 = (double) pow( (double)misplacement_Y, (double)4.0 );
		y3 = (double) pow( (double)misplacement_Y, (double)3.0 );
		y2 = (double) pow( (double)misplacement_Y, (double)2.0 );
	
		b1 = ((double)theAconstant * x3 * misplacement_Y / 3.0) + 
			 ((double)misplacement_X * y2 / 2.0);

		b2 = ((double)theAconstant * theAconstant * x3 * y2 / 3.0) + 
			 ((double)2.0 * theAconstant * misplacement_X * y3 / 3.0);

		b3 = ((double)theAconstant * misplacement_Y * x3 / 3.0) + 
			 ((double)misplacement_X * y2 / 2.0);

		b4 = ((double)theAconstant * theBconstant * misplacement_Y * x4 / 2.0) + 
			 ((double)theBconstant * x2 * y2 / 2.0);

		b5 = ((double)theAconstant * theAconstant * misplacement_Y * x4 / 2.0) + 
			 ((double)theAconstant * x2 * y2 / 2.0) + 
			 ((double)theAconstant * theBconstant * x3 * y2 / 3.0) + 
			 ((double)2.0 * theBconstant * misplacement_X * y3 / 3.0);

		Bxy = ((double)omegaA * b1 * omegaA * 2.0 * theAconstant) + 
			  ((double)omegaA * b2 * omegaAB) +
			  ((double)omegaA * b3 * omegaB * 2.0 * theBconstant) +
			  ((double)omegaA * b4 * omegaAB) +
			  ((double)omegaA * b5 * 2.0 * omegaAB);

		c1 = ((double)theBconstant * misplacement_X * y3 / 3.0) - 
			 ((double)misplacement_Y * x2 / 2.0);

		c2 = ((double)theAconstant * theBconstant * misplacement_X * y4 / 2.0) - 
			 ((double)theAconstant * x2 * y2 / 2.0);

		c3 = ((double)theBconstant * misplacement_X * y3 / 3.0) - 
			 ((double)x2 * misplacement_Y / 2.0);

		c4 = ((double)theBconstant * theBconstant * x2 * y3 / 3.0) - 
			 ((double)2.0 * theBconstant * misplacement_Y * x3 / 3.0);

		c5 = ((double)theAconstant * theBconstant * x2 * y3 / 3.0) - 
			 ((double)2.0 * theAconstant * misplacement_Y * x3 / 3.0) + 
			 ((double)theBconstant * theBconstant * misplacement_X * y4 / 2.0) - 
			 ((double)theBconstant * x2 * y2 / 2.0);

		Cxy = ((double)omegaB * c1 * omegaA * 2.0 * theAconstant) + 
			  ((double)omegaB * c2 * omegaAB) +
			  ((double)omegaB * c3 * omegaB * 2.0 * theBconstant) +
			  ((double)omegaB * c4 * omegaAB) +
			  ((double)omegaB * c5 * 2.0 * omegaAB);

		d1 = ((double)theAconstant * x3 * y2 / 6.0) + 
			 ((double)theBconstant * y3 * x2 / 6.0);

		d2 = ((double)2.0 * theAconstant * theAconstant * x3 * y3 / 9.0) + 
			 ((double)theAconstant * theBconstant * y4 * x2 / 4.0);

		d3 = ((double)theAconstant * x3 * y2 / 6.0) + 
			 ((double)theBconstant * y3 * x2 / 6.0);

		d4 = ((double)theAconstant * theBconstant * x4 * y2 / 4.0) + 
			 ((double)2.0 * theBconstant * theBconstant * y3 * x3 / 9.0);

		d5 = ((double)theAconstant * theAconstant * x4 * y2 / 4.0) + 
			 ((double)2.0 * theAconstant * theBconstant * y3 * x3 / 9.0) + 
			 ((double)2.0 * theAconstant * theBconstant * x3 * y3 / 9.0) + 
			 ((double)theBconstant * theBconstant * y4 * x2 / 4.0);

		Dxy = ((double)omegaAB * d1 * omegaA * 2.0 * theAconstant) + 
			  ((double)omegaAB * d2 * omegaAB) +
			  ((double)omegaAB * d3 * omegaB * 2.0 * theBconstant) +
			  ((double)omegaAB * d4 * omegaAB) +
			  ((double)omegaAB * d5 * 2.0 * omegaAB);

		if ( (double)Fxy != 0.0 )
			 SRI.pointer->ResilientSignal[i1][i2] = ((double)Bxy + Cxy + Dxy) / ((double)Fxy);
		else SRI.pointer->ResilientSignal[i1][i2] =  (double)SRI.pointer->OriginalSignal[i1][i2];
		// calculate the signal resilient to interpolation (end)

		// calculate the curvature of the Resilient Signal (begin) 
		e1 = ((double)theAconstant * x2 * misplacement_Y) + 
			 ((double)y2 / 2.0);

		e2 = ((double)theAconstant * theAconstant * x2 * y2) + 
			 ((double)2.0 * theAconstant * y3 / 3.0);

		e3 = ((double)theAconstant * misplacement_Y * x2) + 
			 ((double)y2 / 2.0);

		e4 = ((double)2.0 * theAconstant * theBconstant * misplacement_Y * x3) + 
			 ((double)theBconstant * misplacement_X * y2);

		e5 = ((double)2.0 * theAconstant * theAconstant * misplacement_Y * x3) + 
			 ((double)theAconstant * misplacement_X * y2) + 
			 ((double)theAconstant * theBconstant * x2 * y2) + 
			 ((double)2.0 * theBconstant * y3 / 3.0);

		dBxydx = ((double)omegaA * e1 * omegaA * 2.0 * theAconstant) + 
			     ((double)omegaA * e2 * omegaAB) +
			     ((double)omegaA * e3 * omegaB * 2.0 * theBconstant) +
			     ((double)omegaA * e4 * omegaAB) +
			     ((double)omegaA * e5 * 2.0 * omegaAB);
		
		f1 = ((double)theAconstant * x3 / 3.0) + 
			 ((double)misplacement_X * misplacement_Y);

		f2 = ((double)2.0 * theAconstant * theAconstant * x3 * misplacement_Y / 3.0) + 
			 ((double)2.0 * theAconstant * misplacement_X * y2);

		f3 = ((double)theAconstant * x3 / 3.0) +  
			 ((double)misplacement_X * misplacement_Y);

		f4 = ((double)theAconstant * theBconstant * x4 / 2.0) + 
			 ((double)theBconstant * x2 * misplacement_Y);

		f5 = ((double)theAconstant * theAconstant * x4 / 2.0) + 
			 ((double)theAconstant * x2 * misplacement_Y) + 
			 ((double)2.0 * theAconstant * theBconstant * x3 * misplacement_Y / 3.0) + 
			 ((double)2.0 * theBconstant * misplacement_X * y2);

		dBxydy = ((double)omegaA * f1 * omegaA * 2.0 * theAconstant) + 
			     ((double)omegaA * f2 * omegaAB) +
			     ((double)omegaA * f3 * omegaB * 2.0 * theBconstant) +
			     ((double)omegaA * f4 * omegaAB) +
			     ((double)omegaA * f5 * 2.0 * omegaAB);

		g1 = ((double)theBconstant * y3 / 3.0) - 
			 ((double)misplacement_X * misplacement_Y);

		g2 = ((double)theAconstant * theBconstant * y4 / 2.0) - 
			 ((double)theAconstant * misplacement_X * y2);

		g3 = ((double)theBconstant * y3 / 3.0) -  
			 ((double)misplacement_X * misplacement_Y);

		g4 = ((double)2.0 * theBconstant * theBconstant * misplacement_X * y3 / 3.0) - 
			 ((double)2.0 * theBconstant * misplacement_Y * x2);

		g5 = ((double)2.0 * theAconstant * theBconstant * misplacement_X * y3 / 3.0) - 
			 ((double)2.0 * theAconstant * misplacement_Y * x2) + 
			 ((double)theBconstant * theBconstant * y4 / 2.0) - 
			 ((double)theBconstant * misplacement_X * y2);

		dCxydx = ((double)omegaB * g1 * omegaA * 2.0 * theAconstant) + 
			     ((double)omegaB * g2 * omegaAB) +
			     ((double)omegaB * g3 * omegaB * 2.0 * theBconstant) +
			     ((double)omegaB * g4 * omegaAB) +
			     ((double)omegaB * g5 * 2.0 * omegaAB);

		h1 = ((double)theBconstant * misplacement_X * y2) - 
			 ((double)x2 / 2.0);

		h2 = ((double)2.0 * theAconstant * theBconstant * misplacement_X * y3) - 
			 ((double)theAconstant * x2 * misplacement_Y);

		h3 = ((double)theBconstant * misplacement_X * y2) -  
			 ((double)x2 / 2.0);

		h4 = ((double)theBconstant * theBconstant * x2 * y2) - 
			 ((double)2.0 * theBconstant * x3 / 3.0);

		h5 = ((double)theAconstant * theBconstant * x2 * y2) - 
			 ((double)2.0 * theAconstant * x3 / 3.0) + 
			 ((double)2.0 * theBconstant * theBconstant * misplacement_X * y3) - 
			 ((double)theBconstant * x2 * misplacement_Y);
		
		dCxydy = ((double)omegaB * h1 * omegaA * 2.0 * theAconstant) + 
			     ((double)omegaB * h2 * omegaAB) +
			     ((double)omegaB * h3 * omegaB * 2.0 * theBconstant) +
			     ((double)omegaB * h4 * omegaAB) +
			     ((double)omegaB * h5 * 2.0 * omegaAB);
	    
		j1 = ((double)theAconstant * x2 * y2 / 2.0) +
			 ((double)theBconstant * y3 * misplacement_X / 3.0);

		j2 = ((double)2.0 * theAconstant * theAconstant * x2 * y3 / 3.0) + 
			 ((double)theAconstant * theBconstant * y4 * misplacement_X / 2.0);

		j3 = ((double)theAconstant * x2 * y2 / 2.0) +  
			 ((double)theBconstant * y3 * misplacement_X / 3.0) ;

		j4 = ((double)theAconstant * theBconstant * x3 * y2) +
			 ((double)2.0 * theBconstant * theBconstant * y3 * x2 / 3.0);

		j5 = ((double)theAconstant * theAconstant * x3 * y2) + 
			 ((double)2.0 * theAconstant * theBconstant * y3 * x2 / 3.0) + 
			 ((double)2.0 * theAconstant * theBconstant * x2 * y3 / 3.0) + 
			 ((double)theBconstant * theBconstant * y4 * misplacement_X / 2.0);
	
		dDxydx = ((double)omegaAB * j1 * omegaA * 2.0 * theAconstant) + 
			     ((double)omegaAB * j2 * omegaAB) +
			     ((double)omegaAB * j3 * omegaB * 2.0 * theBconstant) +
			     ((double)omegaAB * j4 * omegaAB) +
			     ((double)omegaAB * j5 * 2.0 * omegaAB);
		
		k1 = ((double)theAconstant * x3 * misplacement_Y / 3.0) +
			 ((double)theBconstant * y2 * x2 / 2.0);

		k2 = ((double)2.0 * theAconstant * theAconstant * x3 * y2 / 3.0) + 
			 ((double)theAconstant * theBconstant * y3 * x2);

		k3 = ((double)theAconstant * x3 * misplacement_Y / 3.0) +  
			 ((double)theBconstant * y2 * x2 / 2.0) ;

		k4 = ((double)theAconstant * theBconstant * x4 * misplacement_Y / 2.0) +
			 ((double)2.0 * theBconstant * theBconstant * y2 * x3 / 3.0);

		k5 = ((double)theAconstant * theAconstant * x4 * misplacement_Y / 2.0) + 
			 ((double)2.0 * theAconstant * theBconstant * y2 * x3 / 3.0) + 
			 ((double)2.0 * theAconstant * theBconstant * x3 * y2 / 3.0) + 
			 ((double)theBconstant * theBconstant * y3 * x2);

		dDxydy = ((double)omegaAB * k1 * omegaA * 2.0 * theAconstant) + 
			     ((double)omegaAB * k2 * omegaAB) +
			     ((double)omegaAB * k3 * omegaB * 2.0 * theBconstant) +
			     ((double)omegaAB * k4 * omegaAB) +
			     ((double)omegaAB * k5 * 2.0 * omegaAB);
	
		L6 = ((double)omegaA * 2.0 * theAconstant) + ((double)omegaB * 2.0 * theBconstant);

		L7 = ((double)omegaA * 2.0 * theAconstant * misplacement_Y) + 
			 ((double)omegaAB * theAconstant * y2) +
			 ((double)omegaB * 2.0 * theBconstant * misplacement_Y) +
			 ((double)omegaAB * 2.0 * theBconstant * misplacement_X * misplacement_Y) +
			 ((double)2.0 * omegaAB * 
			 (((double)2.0 * theAconstant * misplacement_X * misplacement_Y) +
			 ((double)theBconstant * y2) ) ); 

		L8 = ((double)omegaA * 2.0 * theAconstant * misplacement_X) + 
			 ((double)omegaAB * 2.0 * theAconstant * misplacement_X * misplacement_Y) +
			 ((double)omegaB * 2.0 * theBconstant * misplacement_X) +
			 ((double)omegaAB * theBconstant * x2) +
			 ((double)2.0 * omegaAB * 
			 (((double)theAconstant * x2) +
			 ((double)2.0 * theBconstant * misplacement_X * misplacement_Y) ) ); 

		dFxydx = ((double)L6 * misplacement_Y) - (double)L7;
		dFxydy = ((double)L6 * misplacement_X) - (double)L8;

		m1 = ((double)2.0 * theAconstant * misplacement_X * misplacement_Y);

		m2 = ((double)2.0 * theAconstant * theAconstant * misplacement_X * y2);

		m3 = ((double)2.0 * theAconstant * misplacement_X * misplacement_Y);

		m4 = ((double)6.0 * theAconstant * theBconstant * misplacement_Y * x2) +
			 ((double)theBconstant * y2);

		m5 = ((double)6.0 * theAconstant * theAconstant * misplacement_Y * x2) + 
			 ((double)theAconstant * y2) + 
			 ((double)2.0 * theAconstant * theBconstant * misplacement_X * y2);

		d2Bxydx2 = ((double)omegaA * m1 * omegaA * 2.0 * theAconstant) + 
			       ((double)omegaA * m2 * omegaAB) +
			       ((double)omegaA * m3 * omegaB * 2.0 * theBconstant) +
			       ((double)omegaA * m4 * omegaAB) +
			       ((double)omegaA * m5 * 2.0 * omegaAB);

		o1 = ((double)omegaA * misplacement_X * omegaA * 2.0 * theAconstant);

		o2 = ((double)2.0 * theAconstant * theAconstant * x3 / 3.0) + 
			 ((double)4.0 * theAconstant * misplacement_X * misplacement_Y);

		o3 = ((double)omegaA * misplacement_X * omegaB * 2.0 * theBconstant);

		o4 = ((double)omegaA * theBconstant * x2 * omegaAB);

		o5 = ((double)theAconstant * x2) + 
			 ((double)2.0 * theAconstant * theBconstant * x3 / 3.0) + 
			 ((double)4.0 * theBconstant * misplacement_X * misplacement_Y);

		d2Bxydy2 = ((double) o1 + ((double)omegaA * o2 * omegaAB) + o3 + o4 + ((double)omegaA * o5 * 2.0 * omegaAB) );

		p1 = -((double)omegaB * misplacement_Y * omegaA * 2.0 * theAconstant);

		p2 = -((double)omegaB * theAconstant * y2 * omegaAB);

		p3 = -((double)omegaB * misplacement_Y * omegaB * 2.0 * theBconstant);

		p4 = ((double)2.0 * theBconstant * theBconstant * y3 / 3.0) - 
			 ((double)4.0 * theBconstant * misplacement_Y * misplacement_X);

		p5 = ((double)2.0 * theAconstant * theBconstant * y3 / 3.0) - 
			 ((double)4.0 * theAconstant * misplacement_Y * misplacement_X) -
			 ((double)theBconstant * y2);

		d2Cxydx2 = ((double) p1 + p2 + p3 + ((double)omegaB * p4 * omegaAB) + ((double)omegaB * p5 * 2.0 * omegaAB) );

		q1 = ((double)omegaB * 2.0 * theBconstant * misplacement_X * misplacement_Y * omegaA * 2.0 * theAconstant);

		q2 = ((double)6.0 * theAconstant * theBconstant * misplacement_X * y2) -
			 ((double)theAconstant * x2);

		q3 = ((double)omegaB * 2.0 * theBconstant * misplacement_X * misplacement_Y * omegaB * 2.0 * theBconstant);

		q4 = ((double)omegaB * 2.0 * theBconstant * theBconstant * x2 * misplacement_Y * omegaAB);

		q5 = ((double)2.0 * theAconstant * theBconstant * x2 * misplacement_Y) + 
			 ((double)6.0 * theBconstant * theBconstant * misplacement_X * y2) -
			 ((double)theBconstant * x2);

		d2Cxydy2 = ((double) q1 + ((double)omegaB * q2 * omegaAB) + q3 + q4 + ((double)omegaB * q5 * 2.0 * omegaAB) );

		r1 = ((double)theAconstant * misplacement_X * y2) +
			 ((double)theBconstant * y3 / 3.0);

		r2 = ((double)4.0 * theAconstant * theAconstant * misplacement_X * y3 / 3.0) + 
			 ((double)theAconstant * theBconstant * y4 / 2.0);

		r3 = ((double)theAconstant * misplacement_X * y2) +  
			 ((double)theBconstant * y3 / 3.0) ;

		r4 = ((double)3.0 * theAconstant * theBconstant * x2 * y2) +
			 ((double)4.0 * theBconstant * theBconstant * y3 * misplacement_X / 3.0);

		r5 = ((double)3.0 * theAconstant * theAconstant * x2 * y2) + 
			 ((double)4.0 * theAconstant * theBconstant * y3 * misplacement_X / 3.0) + 
			 ((double)4.0 * theAconstant * theBconstant * misplacement_X * y3 / 3.0) + 
			 ((double)theBconstant * theBconstant * y4 / 2.0);

		d2Dxydx2 = ((double)omegaAB * r1 * omegaA * 2.0 * theAconstant) + 
			       ((double)omegaAB * r2 * omegaAB) +
			       ((double)omegaAB * r3 * omegaB * 2.0 * theBconstant) +
			       ((double)omegaAB * r4 * omegaAB) +
			       ((double)omegaAB * r5 * 2.0 * omegaAB);
			
		s1 = ((double)theAconstant * x3 / 3.0) +
			 ((double)theBconstant * misplacement_Y * x2);

		s2 = ((double)4.0 * theAconstant * theAconstant * x3 * misplacement_Y / 3.0) + 
			 ((double)3.0 * theAconstant * theBconstant * y2 * x2);

		s3 = ((double)theAconstant * x3 / 3.0) +  
			 ((double)theBconstant * misplacement_Y * x2) ;

		s4 = ((double)theAconstant * theBconstant * x4 / 2.0) +
			 ((double)4.0 * theBconstant * theBconstant * misplacement_Y * x3 / 3.0);

		s5 = ((double)theAconstant * theAconstant * x4 / 2.0) + 
			 ((double)4.0 * theAconstant * theBconstant * misplacement_Y * x3 / 3.0) + 
			 ((double)4.0 * theAconstant * theBconstant * x3 * misplacement_Y / 3.0) + 
			 ((double)3.0 * theBconstant * theBconstant * y2 * x2);

		d2Dxydy2 = ((double)omegaAB * s1 * omegaA * 2.0 * theAconstant) + 
			       ((double)omegaAB * s2 * omegaAB) +
			       ((double)omegaAB * s3 * omegaB * 2.0 * theBconstant) +
			       ((double)omegaAB * s4 * omegaAB) +
			       ((double)omegaAB * s5 * 2.0 * omegaAB);

		L9 = ((double)omegaA * 2.0 * theAconstant) + ((double)omegaB * 2.0 * theBconstant); 
		
		d2Fxydxdy = ((double) L9 - SRI.pointer->ClassicCurvature[i1][i2]);

		d2Fxydydx = (double)d2Fxydxdy;

		L4 = ((double)omegaAB * ((double)2.0 * theBconstant * misplacement_Y)) + 
			 ((double)2.0 * omegaAB * ((double)2.0 * theAconstant * misplacement_Y));

		d2Fxydx2 = -(double)L4;

		L5 = ((double)omegaAB * ((double)2.0 * theAconstant * misplacement_X)) + 
			 ((double)2.0 * omegaAB * ((double)2.0 * theBconstant * misplacement_X));

		d2Fxydy2 = -(double)L5;


		t1 = ((double)omegaA * (theAconstant * x2 + misplacement_Y));

		d2Bxydxdy = ((double)t1 * SRI.pointer->ClassicCurvature[i1][i2]);

		t2 = ((double)omegaB * (theBconstant * y2 - misplacement_X));

		d2Cxydxdy = ((double)t2 * SRI.pointer->ClassicCurvature[i1][i2]);

		t3 = ((double)omegaAB * (((double)theAconstant * x2 * misplacement_Y) + 
			                     ((double)theBconstant * y2 * misplacement_X)) );

		d2Dxydxdy = ((double)t3 * SRI.pointer->ClassicCurvature[i1][i2]);

	 
		d2Bxydydx = ((double)d2Bxydxdy); 
		d2Cxydydx = ((double)d2Cxydxdy); 
		d2Dxydydx = ((double)d2Dxydxdy);

		t4 = ((double)dBxydx + dCxydx + dDxydx);
		t5 = ((double)Bxy + Cxy + Dxy);
		t6 = ((double)d2Bxydx2 + d2Cxydx2 + d2Dxydx2);
		t7 = ((double)2.0 * Fxy * dFxydx);

		d2SRIdx2 =  ((double)dFxydx * t4 - d2Fxydx2 * t5 + Fxy * t6 - dFxydx * t4) * ((double)Fxy * Fxy);
		d2SRIdx2 -= (((double)Fxy * t4 - dFxydx * t5) * (double)t7);

		if ( (double)Fxy != 0.0 ) 
		d2SRIdx2 /= (double)pow((double)Fxy, (double)4.0);
		else d2SRIdx2 = (double)0.0;

		t8  = ((double)dBxydy + dCxydy + dDxydy);
		t9  = ((double)d2Bxydy2 + d2Cxydy2 + d2Dxydy2);
		u1 = ((double)2.0 * Fxy * dFxydy);

		d2SRIdy2 =  ((double)dFxydy * t8 - d2Fxydy2 * t5 + Fxy * t9 - dFxydy * t8) * ((double)Fxy * Fxy);
		d2SRIdy2 -= (((double)Fxy * t8 - dFxydy * t5) * (double)u1);

		if ( (double)Fxy != 0.0 ) 
		d2SRIdy2 /= (double)pow((double)Fxy, (double)4.0);
		else d2SRIdy2 = (double)0.0;

		u2 = ((double)d2Bxydxdy + d2Cxydxdy + d2Dxydxdy);


		d2SRIdxdy =  ((double)dFxydy * t4 - d2Fxydxdy * t5 + Fxy * u2 - dFxydx * t8) * ((double)Fxy * Fxy);
		d2SRIdxdy -= (((double)Fxy * t4 - dFxydx * t5) * (double)u1);

		if ( (double)Fxy != 0.0 ) 
		d2SRIdxdy /= (double)pow((double)Fxy, (double)4.0);
		else d2SRIdxdy = (double)0.0;

		u3 = ((double)d2Bxydydx + d2Cxydydx + d2Dxydydx);

		d2SRIdydx =  ((double)dFxydx * t8 - d2Fxydydx * t5 + Fxy * u3 - dFxydy * t4) * ((double)Fxy * Fxy);
		d2SRIdydx -= (((double)Fxy * t8 - dFxydy * t5) * (double)t7);

		if ( (double)Fxy != 0.0 ) 
		d2SRIdydx /= (double)pow((double)Fxy, (double)4.0);
		else d2SRIdydx = (double)0.0;


		SRI.pointer->ResilientCurvature[i1][i2] = ((double)d2SRIdx2 + d2SRIdy2 + d2SRIdxdy + d2SRIdydx);

		// calculate the curvature of the Resilient Signal (end) 

		/// Calculation of the Intensity-Curvature Functional (begin) ///
		double E0 = (double)SRI.pointer->OriginalSignal[i1][i2] * 
			        (2.0 * theAconstant * omegaA + 2.0 * theBconstant * omegaB);

		double EIN = (double)SRI.pointer->OriginalSignal[i1][i2] * 
			         (omegaA * (2.0 * theAconstant * misplacement_X * misplacement_Y) + 
					  omegaAB * (theAconstant * misplacement_X * misplacement_Y * misplacement_Y) +
					  omegaB * (2.0 * theBconstant * misplacement_X * misplacement_Y) +
					  omegaAB * (theBconstant * misplacement_X * misplacement_X * misplacement_Y) +
					  2.0 * omegaAB * (theAconstant * misplacement_X * misplacement_X * misplacement_Y +
                                       theBconstant * misplacement_X * misplacement_Y * misplacement_Y)

					 ) + (omegaA * (theAconstant * pow(misplacement_X, 3.0) * misplacement_Y/3.0 + 
					                misplacement_X * misplacement_Y * misplacement_Y/2.0) * omegaA * 2.0 * theAconstant)
						 + 
						 (omegaA * (theAconstant * theAconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/3.0 + 
					                theAconstant * misplacement_X * pow(misplacement_Y, 3.0) * 2.0/3.0) * omegaAB
						 ) +
						 (omegaA * (theAconstant * pow(misplacement_X, 3.0) * misplacement_Y/3.0 + 
					                misplacement_X * pow(misplacement_Y, 2.0)/2.0) * omegaB * 2.0 * theBconstant
						 ) +
						 (omegaA * (theAconstant * theBconstant * pow(misplacement_X, 4.0) * misplacement_Y/2.0 + 
					                theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0) * omegaAB
						 ) +
						 (omegaA * (theAconstant * theAconstant * pow(misplacement_X, 4.0) * misplacement_Y/2.0 + 
					                theAconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0 +
									theAconstant * theBconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/3.0 + 
									theBconstant * misplacement_X * pow(misplacement_Y, 3.0) * 2.0/3.0 ) * 2.0 * omegaAB
						 ) + 

						 (omegaB * (theBconstant * pow(misplacement_Y, 3.0) * misplacement_X/3.0 - 
					                pow(misplacement_X, 2.0) * misplacement_Y/2.0) * omegaA * 2.0 * theAconstant)
						 + 
						 (omegaB * (theAconstant * theBconstant * pow(misplacement_Y, 4.0) * misplacement_X/2.0 - 
					                theAconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0) * omegaAB)
						 + 
						 (omegaB * (theBconstant * pow(misplacement_Y, 3.0) * misplacement_X/3.0 - 
					                pow(misplacement_X, 2.0) * misplacement_Y/2.0) * omegaB * 2.0 * theBconstant)
						 + 
						 (omegaB * (theBconstant * theBconstant * pow(misplacement_Y, 3.0) * pow(misplacement_X, 2.0)/3.0 - 
					                theBconstant * pow(misplacement_X, 3.0) * misplacement_Y * 2.0/3.0) * omegaAB)
						 +
						 (omegaB * (theAconstant * theBconstant * pow(misplacement_Y, 3.0) * pow(misplacement_X, 2.0)/3.0 - 
					                theAconstant * pow(misplacement_X, 3.0) * misplacement_Y * 2.0/3.0 + 
									theBconstant * theBconstant * misplacement_X * pow(misplacement_X, 4.0)/2.0 -
									theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 2.0)/2.0) * 2.0 * omegaAB)
						 +

						 (omegaAB * (theAconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/6.0 + 
					                 theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 3.0)/6.0 ) * 
									 omegaA * 2.0 * theAconstant 
						 ) + 
						 (omegaAB * (theAconstant * theAconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0) * 2.0/9.0 + 
					                 theAconstant * theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 4.0)/4.0 ) * 
									 omegaAB 
						 ) + 
						 (omegaAB * (theAconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 2.0)/6.0 + 
					                 theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 3.0)/6.0 ) * 
									 omegaB * 2.0 * theBconstant
						 ) + 
						 (omegaAB * (theAconstant * theBconstant * pow(misplacement_X, 4.0) * pow(misplacement_Y, 2.0)/4.0 + 
					                 theBconstant * theBconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0) * 2.0/9.0 ) * 
									 omegaAB
						 ) + 
						 (omegaAB * (theAconstant * theAconstant * pow(misplacement_X, 4.0) * pow(misplacement_Y, 2.0)/4.0 + 
					                 theAconstant * theBconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0) * 2.0/9.0 +
									 theAconstant * theBconstant * pow(misplacement_X, 3.0) * pow(misplacement_Y, 3.0) * 2.0/9.0 + 
									 theBconstant * theBconstant * pow(misplacement_X, 2.0) * pow(misplacement_Y, 4.0)) * 2.0 * omegaAB 
						 ); 

		SRI.pointer->ictBI[i1][i2] = (double)E0;

		SRI.pointer->ictAI[i1][i2] = (double)EIN;

		if ( EIN != 0.0 ) SRI.pointer->ICF[i1][i2] = ((double)E0/EIN);
		else			  SRI.pointer->ICF[i1][i2] = (double)0.0;
		/// Calculation of the Intensity-Curvature Functional (end) ///

		/// Calculation of the Intensity-Curvature Functional (end) ///

		} else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim	
	
	std::cout << "Signal Resilient to Interpolation (SRI) calculated" << endl;
	std::cout << "Curvature of the SRI calculated" << endl;

	double A, A1, A2;
	// calculation of PIC classic (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 

		A1 = ((double)(n1-(2.0*PAD))*XPixelSize);

		A2 = ((double)(n2-(2.0*PAD))*YPixelSize);

		A = (double)sqrt( pow(A1, (double)2.0) + pow(A2, (double)2.0) );
        
		TheShift = (double) sqrt( ((double)misplacement_X * misplacement_X) + 
                                  ((double)misplacement_Y * misplacement_Y) );
           
		SRI.pointer->ThePixICClassic[i1][i2] = (double)TheShift * 
	    ((double)tan( (double)SRI.pointer->ClassicCurvature[i1][i2] * 2.0 * 3.1415926535897 ) / A); 
	

		} else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim	
	// calculation of PIC classic (end)


	std::cout << "PIC Classic calculated" << endl;

	// calculation of PIC resilient (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 

		A1 = ((double)(n1-(2.0*PAD))*XPixelSize);

		A2 = ((double)(n2-(2.0*PAD))*YPixelSize);

		A = (double)sqrt( pow(A1, (double)2.0) + pow(A2, (double)2.0) );
        
		TheShift = (double) sqrt( ((double)misplacement_X * misplacement_X) + 
                                  ((double)misplacement_Y * misplacement_Y) );
           
		SRI.pointer->ThePixICResilient[i1][i2] = (double)TheShift * 
	    ((double)tan( (double)SRI.pointer->ResilientCurvature[i1][i2] * 2.0 * 3.1415926535897 ) / A); 


		} else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim	
	// calculation of PIC resilient (end)

	std::cout << "PIC Resilient calculated" << endl;

	// reconstruct the signal -> Classic (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim

		for (int i2=0; i2 < n2; i2++) {// y dim

		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 

		SRI.pointer->shifted_OriginalSignal[i1][i2] = 
			
		(double)SRI.pointer->OriginalSignal[i1][i2] + 

		(double)SRI.pointer->ThePixICClassic[i1][i2];
		
		} else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim	
	// reconstruct the signal -> Classic (end)


	// reconstruct the signal -> Resilient (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 

		SRI.pointer->the_SRI_shifted_OriginalSignal[i1][i2] = 
			
		(double)SRI.pointer->OriginalSignal[i1][i2] + 

		(double)SRI.pointer->ThePixICResilient[i1][i2];
		
		} else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim	
	// reconstruct the signal -> Resilient (end)

	// reconstruct Hybrid Signal (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 


		if ( (double)abs( (double)SRI.pointer->ThePixICClassic[i1][i2] ) <= 
			
			 (double)abs( (double)SRI.pointer->ThePixICResilient[i1][i2] ) )
		{ 
    
		SRI.pointer->HybridSignal[i1][i2] = (double)SRI.pointer->shifted_OriginalSignal[i1][i2];
			
		} else 
		if ( (double)abs( (double)SRI.pointer->ThePixICClassic[i1][i2] ) > 
			
			 (double)abs( (double)SRI.pointer->ThePixICResilient[i1][i2] ) )
		{ 

		SRI.pointer->HybridSignal[i1][i2] = (double)SRI.pointer->the_SRI_shifted_OriginalSignal[i1][i2];

		}

		} else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim	
   // reconstruct Hybrid Signal (end)


	std::cout << "Signals reconstructed" << endl;

	// compute SE (begin)
	double SE_sum = 0.0;
	long int counter = 0;
	
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 

		counter++;

		SRI.pointer->SE[i1][i2] = (double) sqrt( ( SRI.pointer->OriginalSignal[i1][i2] - SRI.pointer->shifted_OriginalSignal[i1][i2] ) * 
			                                     ( SRI.pointer->OriginalSignal[i1][i2] - SRI.pointer->shifted_OriginalSignal[i1][i2] ) );  

		SE_sum = SRI.pointer->SE[i1][i2] + SE_sum;

		}

		} // y dim
        
	}  // x dim	
	// compute SE (end)

	std::cout << "SE with no SRI is: " << ((double)SE_sum / (counter))  << endl;
	std::cout << endl;

	fprintf(savedata,"%s%e\n", "SE with no SRI is: ", ((double)SE_sum / (counter)) );
	fprintf(savedata,"\n");


	// compute SE_SRI (begin)
	double SE_SRI_sum = 0.0;
	counter = 0;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 

		counter++;

		SRI.pointer->SE_SRI[i1][i2] = (double) sqrt( ( SRI.pointer->OriginalSignal[i1][i2] - SRI.pointer->the_SRI_shifted_OriginalSignal[i1][i2] ) * 
			                                         ( SRI.pointer->OriginalSignal[i1][i2] - SRI.pointer->the_SRI_shifted_OriginalSignal[i1][i2] ) );  

		SE_SRI_sum = SRI.pointer->SE_SRI[i1][i2] + SE_SRI_sum;

		}

		} // y dim
        
	}  // x dim
	// compute SE_SRI (end)

	std::cout << "SE with SRI is: " << ((double)SE_SRI_sum / (counter))  << endl;
	std::cout << endl;

	fprintf(savedata,"%s%e\n", "SE with SRI is: ", ((double)SE_SRI_sum / (counter)) );
	fprintf(savedata,"\n");

	SRI.save();

	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	fclose(savedata);

	free(OMEGAA);
	free(OMEGAB);
	free(OMEGAAB);

	delete SRI.pointer;
	SRI.~SRI2D2013();
	} // processing (end)

	} // run the program (end)

	
	return 0;
} // end of main 