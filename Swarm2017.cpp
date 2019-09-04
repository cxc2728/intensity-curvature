//  Project Title: Swarm2017 (2D Filter)
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

using namespace std;

#define ft_SCALE 255 
// this constant is necessary in order to control the 
// contrast brightness of the Fourier Transformations 
// (direct and inverse)

void OnFourierTransform(char imageFilename[], int rcxres, int rcyres);
void OnInverseFourierTransform(char filename[], int rcyres, int rcxres);

class Swarm2017 {

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

		double **SwarmH; // pointer to the matrix entry

		double **SwarmHPF; // pointer to the particle velocity

		double **particleBestPosition; // pointer to the particle best position

	    double **ParticlePosition; // pointer to the particle position

		double **g; // swarm's best known position

		double **store_particleBestPosition; // storage 

		double **store_g; // storage
		
	}*pointer; // pointer to the matrices

public:

	Swarm2017(int x, int y) : n1(x), n2(y) { }; // constructor 

	void allocateData();

	void save();

	~Swarm2017() { } // destructor

};

void Swarm2017::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;
			
	 pointer->Signal = new double*[this->n1];

	 pointer->SwarmH = new double*[this->n1]; 

	 pointer->SwarmHPF = new double*[this->n1];

	 pointer->particleBestPosition = new double*[this->n1];

	 pointer->ParticlePosition = new double*[this->n1];
	
	 pointer->g = new double*[this->n1];

	 pointer->store_particleBestPosition = new double*[this->n1];

	 pointer->store_g = new double*[this->n1];


	 for( int v=0; v < this->n1; v++ ) { // (1)

		 pointer->Signal[v] = new double[this->n2];

		 pointer->SwarmH[v] = new double[this->n2];

		 pointer->SwarmHPF[v] = new double[this->n2];

		 pointer->particleBestPosition[v] = new double[this->n2];

		 pointer->ParticlePosition[v] = new double[this->n2];

		 pointer->g[v] = new double[this->n2];

		 pointer->store_particleBestPosition[v] = new double[this->n2];

		 pointer->store_g[v] = new double[this->n2];

	  } // (1) allocate struct 'data' (end)

					
		for( int v=0; v < this->n1; v++ ) { // (a)

			for( int f=0; f < this->n2 ; f++ ) { // (b)
				
			 		pointer->Signal[v][f] = (double)0.0;

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


void Swarm2017::save() { // saveImages

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

	char outputFile[128]="Swarm-2017.log";

	FILE * savedata;

if (argc < 9) {  std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the cut-off frequency (double) > 0 and < 1" << endl;
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
		exit(0);

	} else  { // processing (begin)

	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);

	char imageFileName[128];

	sprintf(imageFileName, "%s", argv[1]);

	double cutOFF_frequency = atof(argv[4]);

	double tolerance = atof(argv[5]);
	double omega = atof(argv[6]);
	double Phi_P = atof(argv[7]);
	double Phi_G = atof(argv[8]);
	double error;

	
	std::cout << "The image file name is: " << imageFileName << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[3]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[2]) << endl;
	std::cout << "The cut-off frequency is: " << atof(argv[4]) << endl;
	std::cout << "The value of the tolerance is: " << atof(argv[5]) << endl;
	std::cout << "The value of omega is: " << atof(argv[6]) << endl;
	std::cout << "The value of Phi_P is: " << atof(argv[7]) << endl;
	std::cout << "The value of Phi_G is: " << atof(argv[8]) << endl;


	fprintf(savedata,"%s%s\n", "The image file name is: " , imageFileName);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n2);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n1);
	fprintf(savedata,"%s%lf\n", "The cut-off frequency is: ", atof(argv[4]));
	fprintf(savedata,"%s%lf\n", "The value of the tolerance is: ", atof(argv[5]));
	fprintf(savedata,"%s%lf\n", "The value of omega is: ", atof(argv[6]));
	fprintf(savedata,"%s%lf\n", "The value of Phi_P is: ", atof(argv[7]));
	fprintf(savedata,"%s%lf\n", "The value of Phi_G is: ", atof(argv[8]));
	
	Swarm2017 Swarm(n1, n2); 

	Swarm.allocateData();

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
		
		Swarm.pointer->Signal[i1][i2] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	OnFourierTransform(imageFileName, n2, n1);
	OnInverseFourierTransform(imageFileName, n2, n1);

    double highPassFilter, RC;
	double pi = 3.141592;
	double x = 0.0, y = 0.0, w = 0.0;
	double deltaT = 1.0;
	double condition_fpi_fg;
	
	// Particle swarm optimization (begin)
	// www.en.wikipedia.org/wiki/Particle_swarm_optimization (begins)

	// initialize (begins)
	// Initialize the particle's position with the signal + random number:
	for (int i1=0; i1 < n1; i1++) { // x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
		
			Swarm.pointer->ParticlePosition[i1][i2] = (double)Swarm.pointer->Signal[i1][i2] + ((double)rand() / (double) RAND_MAX); // xi
			

			} // for each dimension (y)

		} // for each dimension (x)

	    // Initialize the particle's best known position to its initial position (the signal + random number):
    	for (int i1=0; i1 < n1; i1++) {// x dim
		   	
			for (int i2=0; i2 < n2; i2++) { // y dim
		
				Swarm.pointer->particleBestPosition[i1][i2] = (double) Swarm.pointer->ParticlePosition[i1][i2]; // pi <- xi
	
			   // initialize 'g': the swarm's best known position (the signal + random number):
			   Swarm.pointer->g[i1][i2] = (double) Swarm.pointer->Signal[i1][i2] + ((double) rand() / (double) RAND_MAX);

			  
			   // condition calculated for if statement (1)
			   condition_fpi_fg = (double) Swarm.pointer->particleBestPosition[i1][i2] - (double) Swarm.pointer->g[i1][i2]; 
			   // expected to be negative
			
			   	if( condition_fpi_fg <= (double) 0.0 ) { // if (1)

			    // Update the swarm's best known position:
				Swarm.pointer->g[i1][i2] = Swarm.pointer->particleBestPosition[i1][i2]; // g <- pi

	         } // if (1)

			} // for each dimension (y)

		} // for each dimension (x)
		// initialize (ends)

			
    // Initialize the particle's velocity (SwarmHPF) with the signal plus a random number
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		 Swarm.pointer->SwarmHPF[i1][i2] = (double) Swarm.pointer->Signal[i1][i2] + ((double) rand() / (double) RAND_MAX); // vi
	
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

					Swarm.pointer->SwarmH[i1][i2] = 
						
				    ((double)highPassFilter * Swarm.pointer->SwarmH[i1-1][i2]) + 
				   (((double) Phi_P * random_p) * 
				    ((double)Swarm.pointer->particleBestPosition[i1][i2] - (double)Swarm.pointer->ParticlePosition[i1-1][i2])) +
				   (((double) Phi_G * random_g) * 
				   ((double)Swarm.pointer->g[i1][i2] - (double)Swarm.pointer->ParticlePosition[i1-1][i2])); // vi <-
					
					/* equation of the HPF 
					//Swarm.pointer->SwarmH[i1][i2] = ((double)highPassFilter * y ) + 
				                                     (((double)highPassFilter) * ((double)x - w )); */
					
					//Update the particle's position:
					Swarm.pointer->ParticlePosition[i1][i2] += (double) Swarm.pointer->SwarmH[i1][i2]; // xi <- xi + vi
														   	
					// expected to be negative:
					double condition_fxi_fpi = (double)Swarm.pointer->ParticlePosition[i1][i2] - 
						                       (double)Swarm.pointer->particleBestPosition[i1][i2];

					// expected to be negative:
					double condition_fpi_fg  = (double)Swarm.pointer->particleBestPosition[i1][i2] - 
						                       (double)Swarm.pointer->g[i1][i2];
 
					if( condition_fxi_fpi <= (double)0.0 ) { // if (2)

						// Update the particle's best known position using the Delta Rule (first order partial derivative of the cost function (**)
						Swarm.pointer->particleBestPosition[i1][i2] = Swarm.pointer->store_particleBestPosition[i1][i2] -
							                                          (double) 2.0 * ((double) Swarm.pointer->SwarmH[i1][i2] - 
						                                              (double) Swarm.pointer->SwarmHPF[i1][i2]) *
																     ((double) Phi_P * random_p);
					}

					if( condition_fpi_fg <= (double)0.0) { // if (3)

						// Update the swarm's best known position using the Delta Rule (first order partial derivative of the cost function (**)   
						Swarm.pointer->g[i1][i2] = Swarm.pointer->store_g[i1][i2] - 
							                       (double) 2.0 * ((double) Swarm.pointer->SwarmH[i1][i2] - 
						                           (double) Swarm.pointer->SwarmHPF[i1][i2]) *
									              ((double) Phi_G * random_g);
					}

				} // for each dimension (y)

			} // for each dimension (x)

				
					// calculate the error of fit (cost function (**) (begins)
					error = 0.0;
					for (int i1=0; i1 < n1; i1++) {// x dim
       	
						for (int i2=0; i2 < n2; i2++) { // y dim

					error += pow( ((double) Swarm.pointer->SwarmH[i1][i2] - 
						           (double) Swarm.pointer->SwarmHPF[i1][i2]), 2.0); 
					
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
						
							Swarm.pointer->SwarmHPF[i1][i2] = (double)Swarm.pointer->SwarmH[i1][i2];

							} // for each dimension (y)

					} // for each dimension (x)
					// update the value of particleVelocity (ends)

					// store particle Best Position and store g; of the current iteration so that can be used at the next iteration
			    	for (int i1=0; i1 < n1; i1++) {// x dim
       	
							for (int i2=0; i2 < n2; i2++) { // y dim
					
							Swarm.pointer->store_particleBestPosition[i1][i2] = (double)Swarm.pointer->particleBestPosition[i1][i2];
							
							Swarm.pointer->store_g[i1][i2] = (double)Swarm.pointer->g[i1][i2];

							} // for each dimension (y)

					} // for each dimension (x)
		
		} while (error > tolerance && Iterations < 100); // do while loop ends here
	   // www.en.wikipedia.org/wiki/Particle_swarm_optimization (ends)
	   // Particle swarm optimization (ends)

	std::cout << "Filter calculated" << endl;

	char FTfilename[128];

	// scale data (begins)
	double MAX = 5000000000000000000.0;
	// SwarmHPF
	double max=-MAX;
	double min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
	
		if( Swarm.pointer->SwarmHPF[i1][i2] > (double)max ) 
			
			max = (double)Swarm.pointer->SwarmHPF[i1][i2];
              
		if( Swarm.pointer->SwarmHPF[i1][i2] < (double)min ) 
			
			min = (double)Swarm.pointer->SwarmHPF[i1][i2];
		

		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           if ( max == min ) Swarm.pointer->SwarmHPF[i1][i2] = (double)0.0;

           else Swarm.pointer->SwarmHPF[i1][i2] = (double) ft_SCALE * (min - Swarm.pointer->SwarmHPF[i1][i2]) / (min - max) ;
			
		} // y dim
        
	}  // x dim 	
	///  scale data (ends)
	
	Swarm.save();

	/// calculate the K-space and the magnitude of the K-Space after filtering (begins)
    sprintf(FTfilename, "%s", "SwarmHPF.img");	
	OnFourierTransform(FTfilename, n2, n1);
	OnInverseFourierTransform(FTfilename, n2, n1);
	/// calculate the K-space and the magnitude of the K-Space after filtering (ends)

	

	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	fclose(savedata);
	
	delete Swarm.pointer;
	Swarm.~Swarm2017();
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
