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

int OnChebyshev(int npoles, char imageFilename[], int m_rcxresFilter, int m_rcyresFilter, double m_eps, double m_cof, double Gain);

int main ( int argc, char * argv[] ) {

	char outputFile[128]="Chebyshev2017.log";

if (argc < 8) { std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the number of Poles of the Chebyshev Filter (integer):" << endl; 
				 std::cout << "Please use a value between 2 and 8" << endl;
				 std::cout << "Please enter the ripple factor (double)" << endl; 
				 std::cout << "Please enter the Gain of the filter (double)" << endl; 
				 std::cout << "Please enter the cut-off frequency (double) ] 0, 1 [" << endl; 
				 std::cout << endl;
				 std::cout << endl;
				 exit(0); }

else { // run the program (begin)

	char imageFileName[128];
	sprintf(imageFileName, "%s", argv[1]);

	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);

	int npoles = atoi(argv[4]);

	double m_eps = atof(argv[5]);
	double Gain = atof(argv[6]);
	double m_cof = atof(argv[7]);

	int r = 1;
	char DataType[128];
	sprintf(DataType, "%s", "DOUBLE");
	FILE * savedata;
	char outputFile[128]="DATA.log";
	

	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		printf("Cannot open output file, Now Exit...");

	} else  { // processing (begin)

	printf("\n");
	printf("The image file name is: %s \n", imageFileName);
	printf("The DataType is: %s \n" , DataType);
	printf("The number of pixels along the X direction is: %d \n" , n1);
	printf("The number of pixels along the Y direction is: %d \n" , n2);
	
	fprintf(savedata,"%s%s\n", "The image file name is: \n" , imageFileName);
	fprintf(savedata,"%s%s\n", "The DataType is: \n" , DataType);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: \n", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: \n", n2);
	
	fprintf(savedata,"\n");

	if ( (r = OnChebyshev(npoles, imageFileName, n2, n1, m_eps, m_cof, Gain)) == 0 )
	{

		cout << "The Chebyshev filter was applied to the input image" << endl;
		
	} else {

		cout << "The Chebyshev filter was not applied to the input image" << endl;
		
	}

	} // processing begin (end)

	return 0;
	} // run the program (end)

}// main

int OnChebyshev(int npoles, char imageFilename[], int m_rcxresFilter, int m_rcyresFilter, double m_eps, double m_cof, double Gain) 
{

	int NofXpixels = npoles, NofYpixels = 1;
	double pi = 3.141592;
	double savedata = 0;
	int i, j, s, p;
	double savedata2 = 0;

	FILE * logfile;
	char logfilename[128]="Chebyshev.log";

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	fprintf(logfile,"%s\n", "Unable to open log file, now exit...");
	exit(0);
	
	} else { // run the Filtering Process (begin) 

	if ( npoles < 2 || npoles > 8 )
	{

		fprintf(logfile,"%s\n", "Unable to Process");
		fprintf(logfile,"%s\n", "Please use a Number of Poles of the Chebyshev Filter between 2 and 8");

		std::cout << "Unable to Process: " << endl; 
		std::cout << "Please use a Number of Poles of the Chebyshev Filter between 2 and 8" << endl;
		exit(0);
	
	} 

		struct data {

		double **ChebyshevFilteredSignal; // pointer to the matrix entry

		double **Signal; // pointer to the matrix entry

	}*data_pointer; // pointer to the matrices

	struct ChebyshevPoles {

		double **spmR; // pointer to the matrix entry

		double **spmI; // pointer to the matrix entry

		double **spm; // pointer to the matrix entry

	}*Chebyshev_pointer; // pointer to the matrices
	

	 // (1) allocate struct 'data' (begin)
	 data_pointer = new data;

	 Chebyshev_pointer = new ChebyshevPoles;


	 data_pointer->ChebyshevFilteredSignal = new double*[m_rcxresFilter];

	 data_pointer->Signal = new double*[m_rcxresFilter];


	 Chebyshev_pointer->spmR = new double*[npoles];

	 Chebyshev_pointer->spmI = new double*[npoles];

	 Chebyshev_pointer->spm = new double*[npoles];


	 for( int v=0; v <= npoles; v++ ) { // (1)

		 Chebyshev_pointer->spmR[v] = new double[1];

		 Chebyshev_pointer->spmI[v] = new double[1];

		 Chebyshev_pointer->spm[v] = new double[1];

	  } // (1) allocate struct 'data' (end)


	  for( int v=0; v < m_rcxresFilter; v++ ) { // (1)
		 
		 data_pointer->ChebyshevFilteredSignal[v] = new double[m_rcyresFilter];

		 data_pointer->Signal[v] = new double[m_rcyresFilter];

	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for( int v=0; v <= npoles; v++ ) { // (a)

			for( int f=0; f < 1 ; f++ ) { // (b)

			Chebyshev_pointer->spmR[v][f] = (double)1.0;

			Chebyshev_pointer->spmI[v][f] = (double)1.0;

			Chebyshev_pointer->spm[v][f] = (double)1.0;

			} //(b)

		 } //(a)
		// (2) initialize (end)

		// (2) initialize (begin)
		for( int v=0; v < m_rcxresFilter; v++ ) { // (a)

			for( int f=0; f < m_rcyresFilter ; f++ ) { // (b)
		 
			data_pointer->ChebyshevFilteredSignal[v][f] = (double)0.0;

			data_pointer->Signal[v][f] = (double)0.0;

			} //(b)

		 } //(a)
		// (2) initialize (end)

	

	FILE * pf;
	double savedata = 0.0;
	char FTfilename[128];
	sprintf(FTfilename, "%s", imageFilename);

    fprintf(logfile, "%s\t%s\n", "Now Reading Signal in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 delete data_pointer;
	 delete Chebyshev_pointer;

	 exit(0);
	
	} else { // read data


	for (int i=0; i<m_rcxresFilter; i++)
	{ ///read image space data
		for (int j=0; j<m_rcyresFilter; j++)
		{

            fread(&savedata,sizeof(double),1,pf);

		    data_pointer->Signal[i][j] = (double)savedata;

		}
	} ///read image space data

	fprintf(logfile,"%s\n", "Signal Read in");

	fclose (pf);
	} // read data

		double eps = 1.0/m_eps; // ripple factor is m_eps
		// calculate Chebyshev Space
		for (i=0; i<=NofXpixels; i++)
		{ 
			
				double thetam = (pi / 2.0) * (2 * i - 1.0) / ( NofXpixels * NofYpixels );
				
				// calculate negative pole only to use in the transfer function
				double argument = ((double)1.0 / ( NofXpixels * NofYpixels )) * 
					              
					   log( (double) eps  + sqrt( ((double) eps*eps ) + 1.0 ) );
				
	
				double sinh = ( (double) exp ( (double) argument ) - (double) exp ( (double) -argument ) ) / 2.0;

				double cosh = ( (double) exp ( (double) argument ) + (double) exp ( (double) -argument ) ) / 2.0;

				Chebyshev_pointer->spmR[i][0] = -(double)sinh * (double) sin ( (double) thetam );

				Chebyshev_pointer->spmI[i][0] = (double)cosh * (double) cos ( (double) thetam ); 

			}
			

	FILE * opf;
	char ofilename[328];

	sprintf(ofilename, "%s", "ChebyshevR.img");

    fprintf(logfile, "%s\t%s\n", "Now Saving Chebyshev Real Space in File: ", ofilename);

    if ((opf = fopen(ofilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save Chebyshev Real Space");

	 delete data_pointer;
	 delete Chebyshev_pointer;

	 exit(0);
	
	} else { // save data


	for (i=0; i<=NofXpixels; i++)
	{ ///save Chebyshev-space data
		for (j=0; j<NofYpixels; j++)
		{

			savedata = (double)Chebyshev_pointer->spmR[i][0];
          
            fwrite(&savedata,sizeof(double),1,opf);

		}
	} ///save Chebyshev-space data

	fprintf(logfile,"%s\n", "Chebyshev Real Space Saved");

	fclose (opf);
	} // save data


	sprintf(ofilename, "%s", "ChebyshevI.img");

    fprintf(logfile, "%s\t%s\n", "Now Saving Chebyshev Imaginary Space in File: ", ofilename);

    if ((opf = fopen(ofilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save Chebyshev Imaginary Space");

	 delete data_pointer;
	 delete Chebyshev_pointer;

	 exit(0);
	
	} else { // save data


	for (i=0; i<=NofXpixels; i++)
	{ ///save Chebyshev-space data
		for (j=0; j<NofYpixels; j++)
		{

			savedata = (double)Chebyshev_pointer->spmI[i][0];
          
            fwrite(&savedata,sizeof(double),1,opf);

		}
	} ///save Chebyshev-space data

	fprintf(logfile,"%s\n", "Chebyshev Imaginary Space Saved");

	fclose (opf);
	} // save data

		/// calculate magnitude of the Poles
		for (s=0; s<=NofXpixels; s++)
		{ 
			for (p=0; p<NofYpixels; p++)
			{
				
				// magnitude of Pole
				Chebyshev_pointer->spm[s][0] = (double) sqrt ( ((double)Chebyshev_pointer->spmR[s][0] * 
																	    Chebyshev_pointer->spmR[s][0] ) + 
				                  		                       ((double)Chebyshev_pointer->spmI[s][0] * 
																	    Chebyshev_pointer->spmI[s][0] ) ); 
		}
	} /// calculate magnitude of the Poles


	
	fprintf(logfile,"%s\n", "Now Filtering Data");

		
		int k = 0;	    
		double highPassFilter;
		double pi = 3.141592;
		double x = 0.0, y = 0.0, w = 0.0;
		double deltaT = 1.0;
		double s = 0.0, dx = 0.0, dy = 0.0;

		int XNEI = (int)2;
        int n7 = ( (int)floor( (double)m_rcxresFilter/2.0) );  
	    int n8 = ( (int)floor( (double)m_rcyresFilter/2.0) );
		
		/// Calculate Chebyshev polynomials & Filter (begin)
		for (int i =-n7+XNEI; i < n7-XNEI; i++) {

            for (int j =-n8; j < n8; j++) {

				dx = (i - m_rcxresFilter/2);
				dy = (j - m_rcyresFilter/2);

				s = ((double)sqrt( (double)dx*dx / (m_rcxresFilter*m_rcyresFilter) + 
					               (double)dy*dy / (m_rcxresFilter*m_rcyresFilter) )  ); 

		        double SCALE = (m_rcxresFilter*m_rcyresFilter);

				s /= (double) SCALE;

				if ( s == 0.0 ) s = 1.0;

		        double product = 1.0;
				
				for (int u=0; u<=NofXpixels; u++)
				{ 
					for (int p=0; p<NofYpixels; p++)
					{
						if ( (double)Chebyshev_pointer->spmR[u][0] < 0 )
						{

							
						product *= ( (double)1.0 / (  (double) s - 
						
						                             ((double)Chebyshev_pointer->spmR[u][0] + 
												   
												      (double)Chebyshev_pointer->spmI[u][0])) );
						}
			
					}
				}
	


				double tf = ((double) Gain / (( pow(2.0, (npoles - 1) ) ) * m_eps)) * ((double)(product*m_cof));

				double encode = (double) exp(-(double)tf) / ( (double) 1.0 + (double) exp(-(double)tf) );

				x = data_pointer->Signal[i + n7][j + n8];
				
				w = data_pointer->Signal[i + n7 - 1][j + n8];

				y = data_pointer->ChebyshevFilteredSignal[i + n7 - 1][j + n8];

				double RC = (double) 1.0 / ((double) 2.0 * pi * encode);

				highPassFilter = ((double) RC) / ((double) deltaT + RC); 
				
				/// Filter (begin)	
				data_pointer->ChebyshevFilteredSignal[i + n7][j + n8] = ((double) highPassFilter * y ) + 
					                                                      (((double)highPassFilter) * ((double)x - w ));
			    /// Filter (ends)
			}
		} //Calculate Chebyshev polynomials & Filter  (end)
		

	char filename[328];

	fprintf(logfile,"%s\n", "Now Saving Filtered Data");
    

	sprintf(filename, "%s%s", "ChebyshevFiltered-", imageFilename);

    fprintf(logfile, "%s\t%s\n", "Now Saving Chebyshev Filtered Data in File: ", filename);

    if ((pf = fopen(filename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save Chebyshev Filtered Data");

	 delete data_pointer;
	 delete Chebyshev_pointer;

	 exit(0);
	
	} else { // save data


	for (int i=0; i<m_rcxresFilter; i++)
	{ ///save FilteredSignal data
		for (int j=0; j<m_rcyresFilter; j++)
		{

			savedata2 = (double)data_pointer->ChebyshevFilteredSignal[i][j];
          
            fwrite(&savedata2,sizeof(double),1,pf);

		}
	} ///save FilteredSignal data

	fprintf(logfile,"%s\n", "Chebyshev Filtered Data Saved");

	fclose (pf);
	} // save data

	delete data_pointer;
	delete Chebyshev_pointer;

	fclose(logfile);
	system( logfilename );
	
	} // run the Filtering Process (begin) 

	return 0;
}// end of function