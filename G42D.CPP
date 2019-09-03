//  This file contains sourcecode distributed as freeware. 
//  The intellectual property of the sourcecode is shown 
//  here to belong to Carlo Ciulla.

// Reference (book):

// [2] Carlo Ciulla "Signal Resilient to Interpolation: An Exploration on the 
// Approximation Properties of the Mathematical Functions" - CreateSpace Publisher - 
// June 2012; pp. 1 - 336, ISBN: 978-1477567487.

// Disclaimer:

// The website here named www.sourcecodewebsiteCarloCiulla.com 2013 [2] does not intend 
// to convey the meaning of profit making for what pertains to the content
// provided. --->>> Instead, when the content is downloaded, the user(s) are
// kindly invited to donate money to charity organizations involved in 
// helping people in need of food and water. <<<---


//  Project Title: Bivariate Cubic SRI Interpolation
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

if (argc < 12) { std::cout << endl;
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
				 std::cout << "Please enter that 'n >= 0' number of pixels along X, Y and Z which will pad the image" << endl;
				 std::cout << "Some data will be lost if you enter a value of n which is greater than" << endl;
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

	double the_A_constant = atof(argv[9]);

	char imageFileName[128];

	int PAD = atoi(argv[10]);

	unsigned char Fyn = 'n';

	Fyn = *argv[11];

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
	std::cout << "The 'n' number of pixels along X and Y which you entered to pad the image is: " << atoi(argv[10]) << endl;
	std::cout << "Some data will be lost if you enter a value of n which is greater than" << endl;
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
	fprintf(savedata,"%s%lf\n", "The a constant is: ", the_A_constant);
	fprintf(savedata,"%s%d\n", "The 'n' number of pixels along X and Y which you entered to pad the image is: ", PAD);
	fprintf(savedata,"%s\n", "Some data will be lost if you enter a value of n which is greater than");
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
	sprintf(command, "%s%s%s%d%s%d%s%c%s%s%s%s", "java -jar Convolutions2013.jar ", 
		                                          imageFileName, " ", n1, " ", n2, 
												  " ", Fyn, " ", "alpha3.img", " ", 
												  "alpha2.img"); 
	system(command);
	/// run JAVA (end)

	
	/// allocate memory & store ALPHAs (begin)
	double * ALPHA3 = 0;
	double * ALPHA2 = 0;

	/// allocate ALPHAs (begin)
	if ((ALPHA3 = (double *) calloc( n1*n2, sizeof(double)) ) == NULL)
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
	SRI.~SRI2D2013();
	exit(0);

	} else { // else


	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
		
		int index = ((int) i2*n1 + i1); 

		*(ALPHA2 + index) = (double)0.0;


		} // y dim
        
	}  // x dim 

	}//else
	/// allocate ALPHAs (end)

	/// read alpha (begin)
	FILE * readAlpha;
	char inputFile[128];
	double readdata;
	
	sprintf(inputFile, "%s","alpha3.img");

	if ((readAlpha = fopen(inputFile,"rb"))==NULL)
	{

		std::cout << "Cannot open input file, Now Exit..." << endl;
		fprintf(savedata,"%s\n", "Cannot open input file, Now Exit...");

		free(ALPHA3);
		free(ALPHA2);
		SRI.~SRI2D2013();
		exit(0);
		
	} else  { // (read)


	for( int i1=0; i1 < n1; i1++ ) { // (a)

		for( int i2=0; i2 < n2; i2++ ) {
	
		int index = ((int) i2*n1 + i1); 

		fread(&readdata, sizeof(double),1, readAlpha);

		*(ALPHA3 + index) = (double)readdata;

		}

	} // (a)

	fclose(readAlpha);

	} // (read)

	sprintf(inputFile, "%s","alpha2.img");

	if ((readAlpha = fopen(inputFile,"rb"))==NULL)
	{

		std::cout << "Cannot open input file, Now Exit..." << endl;
		fprintf(savedata,"%s\n", "Cannot open input file, Now Exit...");

		free(ALPHA3);
		free(ALPHA2);
		SRI.~SRI2D2013();
		exit(0);

	} else  { // (read)


	for( int i1=0; i1 < n1; i1++ ) { // (a)

		for( int i2=0; i2 < n2; i2++ ) {
	
		int index = ((int) i2*n1 + i1); 

		fread(&readdata, sizeof(double),1, readAlpha);

		*(ALPHA2 + index) = (double)readdata;

		}

	} // (a)

	fclose(readAlpha);

	} // (read)
	/// read alpha (end)
	/// allocate memory & store ALPHAs (end)


	std::cout << "Image data scaled" << endl;

    // shift the image of the given misplacement (begin)
	// calculate curvature of model interpolation function (begin)
	double alpha3, alpha2;
	double a, b, c, d, e;
	double G1, G2;
	double TheShift;
	
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 


		int index = ((int) i2*n1 + i1); 

		alpha3 =  (double)*(ALPHA3 + index);

		alpha2 =  (double)*(ALPHA2 + index);

		// compute polynomial coefficients (begin)
		a = ((double) (misplacement_X + misplacement_Y) * 
			          (misplacement_X + misplacement_Y) * 
					  (misplacement_X + misplacement_Y) );

		b = ((double) 0.5 * (misplacement_X + misplacement_Y) * (misplacement_X + misplacement_Y));

		c = ((double) 0.25 * (misplacement_X + misplacement_Y));

		d = ((double) (misplacement_X + misplacement_Y) * (misplacement_X + misplacement_Y));

		e = ((double) 2.0 * (misplacement_X + misplacement_Y));

		// compute polynomial coefficients (end)

		// compute the h functions without transfer function (begin) 
		G1 = ((double) the_A_constant * (a + b + c + 1.0) * alpha2 );

		G2 = ((double) the_A_constant * (d + e + 1.0) * alpha3);
		// compute the h functions without transfer function (end) 

		SRI.pointer->shifted_OriginalSignal[i1][i2] = 
			
		(double)SRI.pointer->OriginalSignal[i1][i2] + ((double)G1 + (double)G2);
		
		
		TheShift = (double)sqrt( ((double)misplacement_X * misplacement_X) + 
			                     ((double)misplacement_Y * misplacement_Y) );

		// calculate curvature of model interpolation function (Classic)
		SRI.pointer->ClassicCurvature[i1][i2] = 
		
		((double) 4.0 * ( alpha2 * the_A_constant * (6.0 * (misplacement_X + misplacement_Y) + 1.0 ) + 
		                                            (2.0 * alpha3 * the_A_constant) ) );

		} else if ( (double)SRI.pointer->OriginalSignal[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim
	
	std::cout << "Image shifted of the given misplacement" << endl;
	std::cout << "Curvature of model interpolation function calculated" << endl;
	// shift the image of the given misplacement (end)
	// calculate curvature of model interpolation function (end)



	// calculate the signal resilient to interpolation (begin)
	double K1, K2, K3, K4, q1;
	double W1, W2, W3, W4, W5;
	double x1, x2, x3, x4, x5;
	double y1, y2, y3, y4, y5;
	double dero2x2, y2dero2; 
	double de2xydero2, yxde2dero2;
	double the_sum, singleTerm; 
	double yNUM, xythe_sum, deDENde;
	double yxdeK4de, yxdeK3de; 
	double yxdeK2de, yxdeK1de, NUM;
	double ydeDENde, theTerm;
	double de2K4dex2, de2K3dex2;
	double de2K2dex2, de2K1dex2; 
	double deK4dex, deK3dex; 
	double deK2dex, deK1dex;
	double y2de2K4de, y2de2K3de; 
	double y2de2K2de, y2de2K1de;
	double ydeK4de, ydeK3de;
	double ydeK2de, ydeK1de;
	double xydeK4de, xydeK3de;
	double xydeK2de, xydeK1de;
	
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim


		if ( (double)SRI.pointer->OriginalSignal[i1][i2] > TH ) { 


		int index = ((int) i2*n1 + i1); 

		alpha3 =  (double)*(ALPHA3 + index);

		alpha2 =  (double)*(ALPHA2 + index);	


		x5 = (double)pow((double)misplacement_X, 5.0);

		x4 = (double)pow((double)misplacement_X, 4.0);

		x3 = (double)pow((double)misplacement_X, 3.0);

		x2 = (double)pow((double)misplacement_X, 2.0);

		x1 = (double)misplacement_X;


		y5 = (double)pow((double)misplacement_Y, 5.0);

		y4 = (double)pow((double)misplacement_Y, 4.0);

		y3 = (double)pow((double)misplacement_Y, 3.0);

		y2 = (double)pow((double)misplacement_Y, 2.0);

		y1 = (double)misplacement_Y;


		K4 = ((double) (1.0/5.0) * y1 * x5 + 
			           (1.0/2.0) * x4 * y2 + 
					   (2.0/3.0) * x3 * y3 + 
					   (1.0/2.0) * x2 * y4 + 
					   (1.0/5.0) * x1 * y5);

		K3 = ((double) (1.0/4.0) * y1 * x4 + 
			           (1.0/6.0) * x3 * y2 + 
					   (1.0/3.0) * x3 * y2 + 
					   (1.0/3.0) * x2 * y3 + 
					   (1.0/6.0) * x2 * y3 + 
					   (1.0/4.0) * x1 * y4);

		K2 = ((double) (1.0/3.0) * y1 * x3 + 
			           (1.0/2.0) * x2 * y2 + 
					   (1.0/3.0) * x1 * y3);

		K1 = ((double) (1.0/2.0) * y1 * x2 + 
					   (1.0/2.0) * x1 * y2);


		W1 = (double) 6.0 * pow( ((double)the_A_constant * alpha2), 2.0); // lambda1

		// lambda2
		W2 = ((double) 4.0 * (double)pow(the_A_constant * alpha2, 2.0) + 
			                 (2.0 * alpha2 * alpha3 * the_A_constant * the_A_constant));

		W2 = (double)W2 + ((double) 6.0 * alpha2 * alpha3 * the_A_constant * the_A_constant);
		// lambda2

		// lambda3
		W3 = ((double) 2.0 * (double)pow( ((double)the_A_constant * alpha2), 2.0) + 
			                 (alpha2 * alpha3 * the_A_constant * the_A_constant));

		W3 = (double)W3 + ((double) (13.0 * alpha2 * alpha3 * the_A_constant * the_A_constant) + 
			                        2.0 * pow( ((double)the_A_constant * alpha3), 2.0) );
		// lambda3

		// lambda4
		W4 = ((double) (25.0/4.0) * (double)pow( ((double)the_A_constant * alpha2), 2.0) + 
			            (1.0/2.0) * (alpha2 * alpha3 * the_A_constant * the_A_constant) );

		W4 = (double)W4 + ((double) (8.0 * alpha2 * alpha3 * the_A_constant * the_A_constant) + 
			                         4.0 * (double)pow( ((double)the_A_constant * alpha3), 2.0) );
		// lambda4


		// lambda5
		W5 = (double)pow( ((double)the_A_constant * alpha2), 2.0) + 
		 	((double) 2.0 * (alpha2 * alpha3 * the_A_constant * the_A_constant));

		W5 = (double)W5 + ((double) (alpha2 * alpha3 * the_A_constant * the_A_constant) + 
			                         2.0 * (double)pow( ((double)the_A_constant * alpha3), 2.0) );
		// lambda5


		// denominator of Eq. (94) (begin)
		q1 = ((double) (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant) * x1 * y1); // lambda6  * x1 * y1;

		q1 = (double)q1 - ((double) (alpha2 * the_A_constant * (6.0 * ((y1 * x2 / 2.0) + (x1 * y2 / 2.0))) + (x1 * y1) ) + 
			                        (2.0 * alpha3 * the_A_constant) * x1 * y1);
		// denominator of Eq. (94) (end)


		if ( ((double)q1) != 0.0 )
			 SRI.pointer->ResilientSignal[i1][i2] = ((double) ( (K4 * W1) + 
			                                                    (K3 * W2) + 
																(K2 * W3) + 
																(K1 * W4) + 
																(W5 * x1 * y1) ) / q1 );

		else SRI.pointer->ResilientSignal[i1][i2] = (double)SRI.pointer->OriginalSignal[i1][i2];
		// calculate the signal resilient to interpolation (end)


		// calculate the curvature of the Resilient Signal (begin) 
		
		deK1dex = ((double) y1 * x1 + (1.0/2.0) * y2); // Eq.(70)

		deK2dex = ((double) y1 * x2 + x1 * y2 + (1.0/3.0) * y3); // Eq. (71)

		deK3dex = ((double) y1 * x3 + (1.0/2.0) * x2 * y2 + x2 * y2 + 
			               (2.0/3.0) * x1 * y3 + (1.0/3.0) * x1 * y3 + (1.0/4.0) * y4); // Eq. (72)

		deK4dex = ((double) y1 * x4 + (2.0 * x3 * y2) + (2.0 * x2 * y3) + 
			                x1 * y4 + (1.0/5.0) * y5); // Eq. (73)


		xydeK1de = ((double) y1 + x1); // Eq.(74)

		xydeK2de = ((double) x2 + 2.0 * x1 * y1 + y2); // Eq.(75)

		xydeK3de = ((double) x3 + x2 * y1 + (2.0 * x2 * y1) + 
			                 (2.0 * x1 * y2) + x1 * y2 + y3); // Eq.(76)

		xydeK4de = ((double) x4 + (4.0 * x3 * y1) + (6.0 * x2 * y2) + 
			                 (4.0 * x1 * y3) + y4); // Eq.(77)

		de2K1dex2 = (double)y1; // Eq.(78)

		de2K2dex2 = ((double) (2.0 * y1 * x1) + y2); // Eq.(79)

		de2K3dex2 = ((double) (3.0 * y1 * x2) + x1 * y2 + 
			                  (2.0 * x1 * y2) + (2.0/3.0) * y3 + (1.0/3.0) * y3); // Eq.(80)

		de2K4dex2 = ((double) (4.0 * y1 * x3) + (6.0 * x2 * y2) + (4.0 * x1 * y3) + y4); // Eq.(81)


		ydeK1de = ((double)y1 * x1 + (1.0/2.0) * x2); // Eq.(82)

		ydeK2de = ((double) (1.0/3.0) * x3 + x2 * y1 + x1 * y2); // Eq.(83)

		ydeK3de = ((double) (1.0/4.0) * x4 + (1.0/3.0) * x3 * y1 + 
			                (2.0/3.0) * x3 * y1 + x2 * y2 + 
				            (1.0/2.0) * x2 * y2 + (x1 * y3)); // Eq.(84)

		ydeK4de = ((double) (1.0/5.0) * x5 + x4 * y1 + (2.0 * x3 * y2) + 
			                (2.0 * x2 * y3) + (x1 * y4)); // Eq.(85)


		yxdeK1de = ((double) y1 + x1); // Eq.(86)

		yxdeK2de = ((double) x2 + 2.0 * x1 * y1 + y2); // Eq.(87)

		yxdeK3de = ((double) x3 + x2 * y1 + (2.0 * x2 * y1) + 
			                (2.0 * x1 * y2) + x1 * y2 + y3); // Eq.(88)

		yxdeK4de = ((double) x4 + (4.0 * x3 * y1) + (6.0 * x2 * y2) + 
			                (4.0 * x1 * y3) + y4); // Eq.(89)


		y2de2K1de = (double)x1; // Eq.(90)

		y2de2K2de = ((double) x2 + (2.0 * y1 * x1)); // Eq.(91)

		y2de2K3de = ((double) (1.0/3.0) * x3 + (2.0/3.0) * x3 + 
			                  (2.0 * x2 * y1) + x2 * y1 + (3.0 * x1 * y2)); // Eq.(92)

		y2de2K4de = ((double) x4 + (4.0 * x3 * y1) + (6.0 * x2 * y2) + (4.0 * x1 * y3)); // Eq.(93)


		// lambda6 = (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant)

		// Eq. (94.a)(begin) -> except the common denominator (the fourth power)
		the_sum = ((double) de2K4dex2 * W1 + de2K3dex2 * W2 + de2K2dex2 * W3 + de2K1dex2 * W4 );
		
		xythe_sum = ((double) (K4 * W1) + (K3 * W2) + (K2 * W3) + (K1 * W4) + (W5 * x1 * y1) );


		// q1 is the denominator to raise to the fourth power
		q1 = ((double) (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant) * x1 * y1); // lambda6  * x1 * y1;

		q1 = (double)q1 - ((double) (alpha2 * the_A_constant * (6.0 * ((y1 * x2 / 2.0) + (x1 * y2 / 2.0))) + (x1 * y1) ) + 
			                        (2.0 * alpha3 * the_A_constant) * x1 * y1);
		// q1 is the denominator to raise to the fourth power

		NUM = ((double) deK4dex * W1 + deK3dex * W2 + deK2dex * W3 + deK1dex * W4 + W5 * y1);

		deDENde = ((double) (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant) * y1); // lambda6 * y1;

		deDENde = (double)deDENde - ((double) (alpha2 * the_A_constant * (6.0 * ((y1 * x1) + (y2 / 2.0))) + y1 ) + 
			                                  (2.0 * alpha3 * the_A_constant) * y1);

		theTerm = ((double) 6.0 * y1 * alpha2 * the_A_constant);

		singleTerm = ((double) ( NUM * q1 ) - ( xythe_sum * deDENde )); // 2011

		dero2x2 = ((double) (( the_sum * q1 + NUM * deDENde - NUM * deDENde - xythe_sum * theTerm ) * q1 * q1) - 
			                ((singleTerm * 2.0 * deDENde) * q1)); //2011
		// Eq. (94.a)(end) 


		// lambda6 = (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant)

		// Eq. (95.a)(begin) -> except the common denominator (the fourth power)                                                        
		the_sum = ((double) y2de2K4de * W1 + y2de2K3de * W2 + y2de2K2de * W3 + y2de2K1de * W4 );

		xythe_sum = ((double) (K4 * W1) + (K3 * W2) + (K2 * W3) + (K1 * W4) + (W5 * x1 * y1) );

		// q1 is the denominator to raise to the fourth power

		q1 = ((double) (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant) * x1 * y1); // lambda6  * x1 * y1;

		q1 = (double)q1 - ((double) (alpha2 * the_A_constant * (6.0 * ((y1 * x2 / 2.0) + (x1 * y2 / 2.0))) + (x1 * y1) ) + 
			                        (2.0 * alpha3 * the_A_constant) * x1 * y1);

		// q1 is the denominator to raise to the fourth power


		NUM = ((double) ydeK4de * W1 + ydeK3de * W2 + ydeK2de * W3 + ydeK1de * W4 + W5 * x1);

		deDENde = ((double) (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant) * x1); // lambda6 * x1;

		deDENde = (double)deDENde - ((double) (alpha2 * the_A_constant * (6.0 * ((x2 / 2.0) + ( x1 * y1 ))) + x1 ) + 
			                                  (2.0 * alpha3 * the_A_constant) * x1);

		theTerm = ((double) 6.0 * x1 * alpha2 * the_A_constant);

		singleTerm = ((double) ( NUM * q1 ) - ( xythe_sum * deDENde )); // 2011


		y2dero2 = ((double) (( the_sum * q1 + NUM * deDENde - NUM * deDENde - xythe_sum * theTerm ) * q1 * q1) - 
			                ((singleTerm * 2.0 * deDENde) * q1)); //2011
		// Eq. (95.a)(end) 


		// lambda6 = (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant)

		// Eq. (96)(begin) -> except the common denominator (the fourth power)
		the_sum = ((double) xydeK4de * W1 + xydeK3de * W2 + xydeK2de * W3 + xydeK1de * W4 + W5);

		xythe_sum = ((double) (K4 * W1) + (K3 * W2) + (K2 * W3) + (K1 * W4) + (W5 * x1 * y1) );


		// q1 is the denominator to raise to the fourth power

		q1 = ((double) (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant) * x1 * y1); // lambda6  * x1 * y1;

		q1 = (double)q1 - ((double) (alpha2 * the_A_constant * (6.0 * ((y1 * x2 / 2.0) + (x1 * y2 / 2.0))) + (x1 * y1) ) + 
			                        (2.0 * alpha3 * the_A_constant) * x1 * y1);

		// q1 is the denominator to raise to the fourth power


		NUM = ((double) deK4dex * W1 + deK3dex * W2 + deK2dex * W3 + deK1dex * W4 + W5 * y1);

		yNUM = ((double) ydeK4de * W1 + ydeK3de * W2 + ydeK2de * W3 + ydeK1de * W4 + W5 * x1);


		deDENde = ((double) (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant) * x1); // lambda6 * x1;

		deDENde = (double)deDENde - ((double) (alpha2 * the_A_constant * (6.0 * ((x2 / 2.0) + ( x1 * y1))) + x1 ) + 
			                                  (2.0 * alpha3 * the_A_constant) * x1);


		ydeDENde = ((double) (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant) * y1); // lambda6 * y1;

		ydeDENde = (double)ydeDENde - ((double) (alpha2 * the_A_constant * (6.0 * ((y1 * x1) + (y2 / 2.0))) + y1 ) + 
												(2.0 * alpha3 * the_A_constant) * y1);


		theTerm = ((double) (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant)); // lambda6;

		theTerm = (double)theTerm - ((double) (alpha2 * the_A_constant * (6.0 * (y1 + x1) + 1.0 )) + 
											  (2.0 * alpha3 * the_A_constant));


		singleTerm = ((double) ( NUM * q1 ) - ( xythe_sum * ydeDENde )); // 2011


		de2xydero2 = ((double) (( the_sum * q1 + NUM * deDENde - yNUM * ydeDENde - xythe_sum * theTerm ) * q1 * q1) - 
							   ((singleTerm * 2.0 * deDENde) * q1)); // 2011

		// Eq. (96)(end)


		// lambda6 = (alpha2 * the_A_constant + 2.0 * alpha3 * the_A_constant)

		// Eq. (97)(begin) -> except the common denominator (the fourth power)

		the_sum = ((double) yxdeK4de * W1 + yxdeK3de * W2 + yxdeK2de * W3 + yxdeK1de * W4 + W5);


		singleTerm = ((double) ( yNUM * q1 ) - ( xythe_sum * deDENde )); // 2011

		yxde2dero2 = ((double) (( the_sum * q1 + yNUM * ydeDENde - NUM * deDENde - xythe_sum * theTerm ) * q1 * q1) - 
			                   ((singleTerm * 2.0 * ydeDENde) * q1)); //2011
	
		// Eq. (97)(end)


		if ( (double)pow( (double)q1 , 4.0 ) != 0 ) 

			 SRI.pointer->ResilientCurvature[i1][i2] = ((double) ( dero2x2 + y2dero2 + de2xydero2 + yxde2dero2 ) / 
			                                                     pow( (double)q1, 4.0) ); 

		else SRI.pointer->ResilientCurvature[i1][i2] = (double)0.0; 

		
		// calculate the curvature of the Resilient Signal (end) 

		/// Calculation of the Intensity-curvature Functional (begin) ///
		double E0 = (double) 4.0 * the_A_constant * SRI.pointer->OriginalSignal[i1][i2] * x1 * y1 * (alpha2 + 2.0 * alpha3);

		double M4 = ((double) 1.0/5.0 * y1 * x5 + 1.0/2.0 * x4 * y2 + 2.0/3.0 * x3 * y3 + 1.0/2.0 * x2 * y4 + 1.0/5.0 * x1 * y5);

		double M3 = ((double)(1.0/4.0 * y1 * x4 + 1.0/6.0 * x3 * y2 + 1.0/3.0 * x3 * y2 + 
			                  1.0/3.0 * x2 * y3 + 1.0/6.0 * x2 * y3 + 1.0/4.0 * x1 * y4));

		double M2 = ((double)(1.0/3.0 * y1 * x3 + 1.0/2.0 * x2 * y2 + 1.0/3.0 * x1 * y3));

		double M1 = ((double) 1.0/2.0 * y1 * x2 + 1.0/2.0 * x1 * y2);

		double EIN = (double) 4.0 * (SRI.pointer->OriginalSignal[i1][i2] * (alpha2 * the_A_constant * (6.0 * (y1 * x2/2.0 + x1 * y2/2.0) + x1 * y1) +
			                  2.0 * alpha3 * the_A_constant * x1 * y1) + 6.0 * M4 * (alpha2 * the_A_constant * alpha2 * the_A_constant) + M3 *
							  (4.0 * (alpha2 * the_A_constant * alpha2 * the_A_constant) + 2.0 * alpha3 * alpha2 * the_A_constant * the_A_constant +
							  6.0 * alpha3 * alpha2 * the_A_constant * the_A_constant) + M2 * (2.0 * (alpha2 * the_A_constant * alpha2 * the_A_constant) + 
							  alpha3 * alpha2 * the_A_constant * the_A_constant + 13.0 * alpha3 * alpha2 * the_A_constant * the_A_constant + 2.0 * 
							  (alpha3 * the_A_constant * alpha3 * the_A_constant)) + M1 * (25.0/4.0 * (alpha2 * the_A_constant * alpha2 * the_A_constant) +
							  1.0/2.0 * alpha3 * alpha2 * the_A_constant * the_A_constant + 8.0 * alpha3 * alpha2 * the_A_constant * the_A_constant + 
							  4.0 * (alpha3 * the_A_constant * alpha3 * the_A_constant)) + ((alpha2 * the_A_constant * alpha2 * the_A_constant) + 
							  2.0 * alpha3 * alpha2 * the_A_constant * the_A_constant + alpha3 * alpha2 * the_A_constant * the_A_constant + 2.0 * 
							  (alpha3 * the_A_constant * alpha3 * the_A_constant) ) * x1 * y1);

		SRI.pointer->ictBI[i1][i2] = (double)E0;

		SRI.pointer->ictAI[i1][i2] = (double)EIN;
		
		if ( EIN != 0.0 ) SRI.pointer->ICF[i1][i2] = ((double)E0/EIN);
		else			  SRI.pointer->ICF[i1][i2] = (double)0.0;
		/// Calculation of the Intensity-curvature Functional (end) ///

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

	free(ALPHA3);
	free(ALPHA2);
	delete SRI.pointer;
	SRI.~SRI2D2013();
	} // processing (end)

	} // run the program (end)

	
	return 0;
} // end of main 