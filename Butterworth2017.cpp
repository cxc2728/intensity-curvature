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

void OnButterworth(int onbp, int m_rcxresFilter, int m_rcyresFilter, char* imageFileName, double m_cof, double Gain);

int main ( int argc, char * argv[] ) {

	char outputFile[128]="Butterworth2017.log";

if (argc < 7) { std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the order of the Normalized Butterworth Polynomials (integer):" << endl; 
				 std::cout << "Please use a value between 2 and 8" << endl;
				 std::cout << "Please enter the Gain of the Bessel Filter (double)" << endl; 
				 std::cout << "Please enter the cut-off frequency (double) ] 0, 1 [" << endl; 
				 std::cout << endl;
				 exit(0); }

else { // run the program (begin)

	char imageFileName[128];
	sprintf(imageFileName, "%s", argv[1]);

	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);
	int onbp = atoi(argv[4]);
	double Gain = atof(argv[5]);
	double m_cof = atof(argv[6]);

	if ( onbp < 2 || onbp > 8 )
	{

	    printf("%s\n", "Unable to Process");
		printf("%s\n", "Please use a value of the order of the Normalized Butterworth Polynomials between 2 and 8");
		exit(0);
	
	} 

	OnButterworth(onbp, n2, n1, imageFileName, m_cof, Gain);
	
	cout << "The Butterworth filter was calculated" << endl;

	return 0;
	} // run the program (end)

}// main

void OnButterworth(int onbp, int m_rcxresFilter, int m_rcyresFilter, char* imageFileName, double m_cof, double Gain) 
{

	int Order = onbp;
	double pi = 3.141592;
	double savedata = 0;
	double savedata2 = 0;

	FILE * logfile;
	char logfilename[128]="Butterworth.log";

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	fprintf(logfile,"%s\n", "Unable to open log file, now exit...");
	exit(0);
	
	} else { // run the Filtering Process (begin) 

	if ( Order < 2 || Order > 8 )
	{

		fprintf(logfile,"%s\n", "Unable to Process");
		fprintf(logfile,"%s\n", "Please use a value of order of the Normalized Butterworth Polynomials between 2 and 8");
		exit(0);
	
	} 

	struct data {

		double **ButterworthFilteredSignal; // pointer to the matrix entry

		double **Signal; // pointer to the matrix entry

	}*data_pointer; // pointer to the matrices

	struct ButterworthPolynomial {

		double **Bw_Odd; // pointer to the matrix entry

		double **Bw_Even; // pointer to the matrix entry

		double **BW; // pointer to the matrix entry

	}*Butterworth_Polynomial_pointer; // pointer to the matrices
	

	 // (1) allocate struct 'data' (begin)
	 data_pointer = new data;

	 Butterworth_Polynomial_pointer = new ButterworthPolynomial;
		
	 data_pointer->ButterworthFilteredSignal = new double*[m_rcxresFilter];

	 data_pointer->Signal = new double*[m_rcxresFilter];

	 Butterworth_Polynomial_pointer->Bw_Odd = new double*[Order+1];

	 Butterworth_Polynomial_pointer->Bw_Even = new double*[Order+1];

	 Butterworth_Polynomial_pointer->BW = new double*[Order+1];


	 for( int v=0; v <= Order; v++ ) { // (1)

		 Butterworth_Polynomial_pointer->Bw_Odd[v] = new double[1];

		 Butterworth_Polynomial_pointer->Bw_Even[v] = new double[1];

		 Butterworth_Polynomial_pointer->BW[v] = new double[1];

	  } // (1) allocate struct 'data' (end)


	  for( int v=0; v < m_rcxresFilter; v++ ) { // (1)
		 
		 data_pointer->ButterworthFilteredSignal[v] = new double[m_rcyresFilter];

		 data_pointer->Signal[v] = new double[m_rcyresFilter];

	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for( int v=0; v <= Order; v++ ) { // (a)

			for( int f=0; f < 1 ; f++ ) { // (b)

			Butterworth_Polynomial_pointer->Bw_Odd[v][f] = (double)1.0;

			Butterworth_Polynomial_pointer->Bw_Even[v][f] = (double)1.0;

			Butterworth_Polynomial_pointer->BW[v][f] = (double)1.0;

			} //(b)

		 } //(a)
		// (2) initialize (end)

		// (2) initialize (begin)
		for( int v=0; v < m_rcxresFilter; v++ ) { // (a)

			for( int f=0; f < m_rcyresFilter ; f++ ) { // (b)
		 
			data_pointer->ButterworthFilteredSignal[v][f] = (double)0.0;

			data_pointer->Signal[v][f] = (double)0.0;

			} //(b)

		 } //(a)
		// (2) initialize (end)
	

	FILE * pf;
	double savedata = 0.0;
	char FTfilename[128];
	sprintf(FTfilename, "%s", imageFileName);

    fprintf(logfile, "%s\t%s\n", "Now Reading Signal in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 delete data_pointer;
	 delete Butterworth_Polynomial_pointer;

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

		
		double Bw_term = 1.0;
		int k = 0;	    
		double highPassFilter;
		double pi = 3.141592;
		double x = 0.0, y = 0.0, w = 0.0;
		double deltaT = 1.0;
		double s = 0.0, dx = 0.0, dy = 0.0;

		int XNEI = (int)2;
        int n7 = ( (int)floor( (double)m_rcxresFilter/2.0) );  
	    int n8 = ( (int)floor( (double)m_rcyresFilter/2.0) );
		
		/// Calculate Butterworth polynomials (begin)
		for (int i =-n7+XNEI; i < n7-XNEI; i++) {

            for (int j =-n8; j < n8; j++) {

				dx = (i - m_rcxresFilter/2);
				dy = (j - m_rcyresFilter/2);

				s = ((double)sqrt( (double)dx*dx / (m_rcxresFilter*m_rcyresFilter) + 
					               (double)dy*dy / (m_rcxresFilter*m_rcyresFilter) )  ); 

		        double SCALE = (m_rcxresFilter*m_rcyresFilter);

				s /= (double) SCALE;

				if ( s == 0.0 ) s = 1.0;

				// Butterworth (even)
				for (int o=1; o<=Order; o++)
				{ 

					if ( o % 2 == 0 )
					{ // if even

					k++;

					Bw_term = (double)s*s - ((double)2.0 * s * cos ( pi * ((double)(2.0 * k + Order - 1.0) / (2.0 * Order)) ) + 1.0); 

					Butterworth_Polynomial_pointer->Bw_Even[o][0] *= ((double)s - Bw_term);

					} // if even
		
				}// Butterworth (even)
		
				k = 0;
				
				

			  	// Butterworth (odd)
				for (int o=1; o<=Order; o++)
				{ 

					if ( o % 2 != 0 )
					{ // if odd

					k++;

					double Bw_term = ((double)s + 1.0) * ((double)s*s - 2.0 * s * cos ( pi * ((double)(2.0 * k + Order - 1.0) / (2.0 * Order)) ) + 1.0); 

					Butterworth_Polynomial_pointer->Bw_Odd[o][0] *= ((double)s - Bw_term);

					} // if odd
				
				} // Butterworth Odd
				
				// calculate Butterworth polynomial (begin)
				double Butterworth_Polynomial = 1.0;

				for (int p=1; p<=Order; p++)
				{ 
				if ( p % 2 == 0 )
				{ 

				Butterworth_Polynomial *= ((double)Butterworth_Polynomial_pointer->Bw_Even[p][0]); 
			
				} else if ( p % 2 != 0 )
				{ 
		
			    Butterworth_Polynomial *= ((double)Butterworth_Polynomial_pointer->Bw_Odd[p][0]);
					         			
				}
				
				}// calculate Butterwoth polynomial (end)	


				double tf = ((double)Gain/(Butterworth_Polynomial/m_cof));

				double encode = (double) exp(-(double)tf) / ( (double) 1.0 + (double) exp(-(double)tf) );

				x = data_pointer->Signal[i + n7][j + n8];
				
				w = data_pointer->Signal[i + n7 - 1][j + n8];

				y = data_pointer->ButterworthFilteredSignal[i + n7 - 1][j + n8];

				double RC = (double) 1.0 / ((double) 2.0 * pi * encode);

				highPassFilter = ((double) RC) / ((double) deltaT + RC); 
				
				/// Filter (begin)	
				data_pointer->ButterworthFilteredSignal[i + n7][j + n8] = ((double) highPassFilter * y ) + 
					                                                      (((double)highPassFilter) * ((double)x - w ));
			    /// Filter (ends)
			}
		}//Calculate Butterwoth polynomials & filter (end)
		

	char filename[328];

	fprintf(logfile,"%s\n", "Now Saving Filtered Data");
    

	sprintf(filename, "%s%s", "ButterworthFiltered-", imageFileName);

    fprintf(logfile, "%s\t%s\n", "Now Saving Butterworth Filtered Data in File: ", filename);

    if ((pf = fopen(filename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save Butterworth Filtered Data");

	 delete data_pointer;
	 delete Butterworth_Polynomial_pointer;

	 exit(0);
	
	} else { // save data


	for (int i=0; i<m_rcxresFilter; i++)
	{ ///save BWFilteredSignal data
		for (int j=0; j<m_rcyresFilter; j++)
		{

			savedata2 = (double)data_pointer->ButterworthFilteredSignal[i][j];
          
            fwrite(&savedata2,sizeof(double),1,pf);

		}
	} ///save BWFilteredSignal data

	fprintf(logfile,"%s\n", "Butterworth Filtered Data Saved");

	fclose (pf);
	} // save data

	delete data_pointer;
	delete Butterworth_Polynomial_pointer;

	fclose(logfile);
	system( logfilename );
	
	} // run the Filtering Process (begin) 

}// end of function