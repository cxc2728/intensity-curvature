#define _CRT_SECURE_NO_WARNINGS
// name of the project: FTinvFT2018
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

using namespace std;

#define ft_SCALE 255

void OnFourierTransform(char imageFilename[], int rcxres, int rcyres);
void OnInverseFourierTransform(char filename[], int rcyres, int rcxres);

// declare a class by the name 'FTinvFT2018'
class FTinvFT2018 {

	// the point is to assign to 'n1' and 'n2' 
	// the correct values from the command console
	int n1; // matrix size x (number of pixels along the x direction)
	int n2; // matrix size y (number of pixels along the y direction)

// declare the class methods
// and embed the methods body
// in the class specification
public:

	// declare a method that returns the number
	// of pixels of the image along the x direction
	int getNofPixelsX(void) { return this->n1; };

	// declare a method that returns the number
	// of pixels of the image along the y direction
	int getNofPixelsY(void) { return this->n2; };

	// declare a method that sets the number
	// of pixels of the image along the x direction
	void setNofPixelsX(int x) { this->n1 = x; };

	// declare a method that sets the number 
	// of pixels of the image along the y direction
	void setNofPixelsY(int y) { this->n2 = y; };

public:
	// declare a structure ' data ' that defines the
	// pointers to pointers to the image
	struct data {

		double **Signal; // declare the pointer to pointer to the matrix (image) 

		}*pointer; // pointer to the element of the structure 'data'
	               // the pointer points to the memory address of the
	               // pointers to pointers 

public:

	// constructor of the class 'FTinvFT2018'
	FTinvFT2018(int x, int y) : n1(x), n2(y) { }; 
		
	// declare the prototype of the 
	// function 'allocateData'
	// the function belongs to the 
	// sets of methods of the class 'FTinvFT2018'
	void allocateData();

	// declare the prototype of the 
	// function 'save'
	// the function belongs to the 
	// sets of methods of the class 'FTinvFT2018'
	void save();

	// destructor of the class 'FTinvFT2018'
	~FTinvFT2018() { } 

};

void FTinvFT2018::allocateData() { // allocate data


	// (1) allocate struct 'data' (begin)
	// define a pointer by the name 'pointer'
	// and assign to it the instance of the
	// structure 'data' (the instance is a
	// memory address
	 pointer = new data;
	
	 // assign to the pointer to a pointer 'Signal' (pointer->Signal)
	 // the instance of the memory address of the pointer to pointer
	 // *[this->n2]. Practically define the memory address of the 
	 // rows of the matrix containing the image 'Signal'.
	 pointer->Signal = new double*[this->n2];

	 
	 for( int v=0; v < this->n2; v++ ) { // (1)
	 // at each iteration of the for loop
	 // assign to the pointer 'Signal[v]' the instance 
	 // of the memory address of the pointer [this-n1].
	 // Practically define the memory address of the 
	 // columns of the matrices containing the image:
	 // 'Signal.

		 pointer->Signal[v] = new double[this->n1];

		  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for(int v=0; v < this->n2; v++ ) { // (a)

			for( int f=0; f < this->n1 ; f++ ) { // (b)
			
			// at each iteration of the two for loops
			// initializes the value of the pixel of
			// the image to zero. This is done for the
			// image 'Signal'.
			pointer->Signal[v][f] = (double)0.0;

			} //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


void FTinvFT2018::save() { // saveImages

	// declare a pointer to file
	// to be used to read the image
	FILE * savedata;
	// declare a string which contains
	// the file name of the image
	char outputFile[128];
	
	// assign the name string "Signal.img"
	// to the string 'outputFile'
	sprintf(outputFile, "%s","Signal.img");

	// open the image file in write-binary mode
	if ((savedata = fopen(outputFile,"wb"))==NULL)
	{
		// alert the user of possible failure in opening the file
		std::cout << "Cannot open output file, Now Exit..." << endl;
		exit(0);

	} else  { // (save)


	for( int v=0; v < this->n2; v++ ) { // (a)

		for( int f=0; f < this->n1; f++ ) 
	
		// at each iteration of the for loop saves the pixel
		// value contained at the memory address: ' &pointer->Signal[v][f] '
		fwrite(&pointer->Signal[v][f],sizeof(double),1,savedata);

	} // (a)

	// close the file after saving
	fclose(savedata);

	} // (save)

	} // saveImages

// read the input parameters from the
// console command line
int main ( int argc, char * argv[] ) { 

	// assign to the char string 'outputFile'
	// the value "FTinvFt2018.log"
	char outputFile[128] = "FTinvFT2018.log";

	// declare a pointer to a file by the
	// name 'savedata'
	FILE * savedata;

// tell the user of the list of input parameters necessary to tun the program:
if (argc < 4) { std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << endl;
				 exit(0); }

else { // run the program (begin)
	
	// opens the log file which name is 
	// contained in 'outputFile', opens 
	// in write mode
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{
		// alert the user of possible failure in opening the file
		std::cout << "Cannot open output file, Now Exit..." << endl;
		exit(0);

	} else  { // processing (begin)

	// declare an array of char to contain a string
	char imageFileName[128];
	
	// transfer into the array 'imageFileName'
	// the image file name as per input from
	// the command console. The image file name
	// is 'argv[1]'
	sprintf(imageFileName, "%s", argv[1]);

	// reads from the command console the
	// value of the image size (number of rows
	// and number of columns)
	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);

	// inform the user of the image size 
	// (number of rows and number of columns
	// of the matrix containing the image)
	std::cout << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;

	// save into the log file the image size 
	// (number of rows and number of columns
	// of the matrix containing the image)
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);

	// call to the constructor 'FT' so to create
	// an object of type 'FT'. The data type of 'FT'
	// is 'FTinvFT2018'
	FTinvFT2018 FT(n1,n2);

	// the object of type 'FT' 
	// sends a message (invokes)
	// to the method 'allocateData()'
	FT.allocateData();

	/// read image file (begin)
	// declare a file pointer
	FILE * pf;

	// open the file containing the image to
	// process. The image to process is by the
	// name of the string contained by 'imageFileName'
	if ((pf = fopen(imageFileName,"rb"))==NULL)
	{
		// alert the user and save into the log file 
		// the event that the program cannot open the 
		// file containing the image to process
		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		exit(0);

	} else { // else

	double number;

	for (int i1=0; i1 < n2; i1++) {// x dim
       	
		for (int i2=0; i2 < n1; i2++) { // y dim
		
		// at each iteration of the two for loops
		// the program reads the pixel value from the
		// file containing the image and 
		fread(&number,sizeof(double),1,pf);
		
		// assigns the pixel value 'number' to the
		// pixel value 'FT.pointer->Signal[i1][i2]'
		// where 'FT' invokes the pointer to the 
		// image 'Signal' at the matrix location '[i1][i2]'
		FT.pointer->Signal[i1][i2] = (double)number;
                          
		} // y dim
        
	}  // x dim 

    // close the file containg the image to process 	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl;  

	// scale Signal.img (begin) 
	double max=-FT.pointer->Signal[0][0];
	double min=FT.pointer->Signal[0][0];

	for (int i1=0; i1<n2; i1++)
	{ 
			for (int i2=0; i2<n1; i2++)
			{

				if( (double)FT.pointer->Signal[i1][i2] > (double)max ) 
			
					max = (double)FT.pointer->Signal[i1][i2];
              
				if( (double)FT.pointer->Signal[i1][i2] < (double)min ) 
			
					min = (double)FT.pointer->Signal[i1][i2];
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1<n2; i1++)
	{ 
			for (int i2=0; i2<n1; i2++)
			{

				if ( max == min ) (double)FT.pointer->Signal[i1][i2];

				else (double)FT.pointer->Signal[i1][i2] = (double) ft_SCALE * (min - (double)FT.pointer->Signal[i1][i2]) / (min - max) ;
			
		} // y dim
        
	}  // x dim
	// scale Signal.img (end)

	// save all of the images
	// the object 'FT' invokes (sends a message)
	// the method by the name of 'save()'
	FT.save(); 

    OnFourierTransform("Signal.img", n1, n2);
	OnInverseFourierTransform("Signal.img", n1, n2);

	// alert the user of the end of the program
	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	// save to log file
	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	fclose(savedata);
    // delete the memorry address
	// allocated to the images
	// do this by the command 'delete'
	// applied to the object 'FT' which
	// invokes the pointer to the data
	// structure 'data' containing the image
	// 'Signal'
	delete FT.pointer;

	// the object 'FT' invokes the
	// class destructor
	FT.~FTinvFT2018();
	} // processing (end)

	} // run the program (end)

	// ANSI C requires the 'main'
	// function returning a value: 
	// zero in this case
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
	system( logfilename );

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
	system( logfilename );
		
	// FIFO memory deallocation method
	free(kSpaceR);
	free(kSpaceI);
	free(reconSignal);

}