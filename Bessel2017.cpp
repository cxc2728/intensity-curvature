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

void OnBessel(int n_poles, double Gain, char * m_FileName, int m_rcxresFilter, int m_rcyresFilter, double m_cof);

int main ( int argc, char * argv[] ) {

	char outputFile[128]="Bessel2017.log";

if (argc < 7) {  std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the number of poles of the Bessel Filter (integer):" << endl; 
				 std::cout << "Please use a value between 1 and 5" << endl;
				 std::cout << "Please enter the Gain of the Bessel Filter (double)" << endl; 
				 std::cout << "Please enter the cut-off frequency (double) ] 0, 1 [" << endl; 
				 std::cout << endl;
				 exit(0); }

else { // run the program (begin)
	
	char imageFileName[128];
	sprintf(imageFileName, "%s", argv[1]);

	int m_rcxres = atoi(argv[2]);
	int m_rcyres = atoi(argv[3]);
	int n_poles = atoi(argv[4]);
	double Gain = atof(argv[5]);
	double m_cof = atof(argv[6]);

	if ( n_poles < 1 || n_poles > 5 )
	{

	    printf("%s\n", "Unable to Process");
		printf("%s\n", "Please use a value of the number of poles between 1 and 5");
		exit(0);
	
	} 

	OnBessel(n_poles, Gain, imageFileName, m_rcyres, m_rcxres, m_cof);
	
	cout << "The Bessel filter was calculated" << endl;

	} // run the program (end)

	
	return 0;
} // end of main 

void OnBessel(int n_poles, double Gain, char * m_FileName, int m_rcxresFilter, int m_rcyresFilter, double m_cof)  
{

	FILE * logfile;
	char logfilename[128]="Bessel.log";

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	fprintf(logfile,"%s\n", "Unable to open log file, now exit...");
	exit(0);
	
	} else { // run the Filtering Process (begin) 

		
		struct data {

		double **BesselFilteredSignal; // pointer to the matrix entry

		double **Signal; // pointer to the matrix entry

	}*data_pointer; // pointer to the matrices

	struct BesselPolynomial {

		double **ak; // pointer to the matrix entry

		double **bk; // pointer to the matrix entry

	}*BesselPolynomial_pointer; // pointer to the matrices


	 // (1) allocate struct 'data' (begin)
	 data_pointer = new data;

	 BesselPolynomial_pointer = new BesselPolynomial;
			
	 data_pointer->BesselFilteredSignal = new double*[m_rcxresFilter];

	 data_pointer->Signal = new double*[m_rcxresFilter];

	 BesselPolynomial_pointer->ak = new double*[n_poles+1];

	 BesselPolynomial_pointer->bk = new double*[n_poles+1];


	 for( int v=0; v <= n_poles; v++ ) { // (1)

		 BesselPolynomial_pointer->ak[v] = new double[1];

		 BesselPolynomial_pointer->bk[v] = new double[1];

	  } // (1) allocate struct 'data' (end)

	  for( int v=0; v < m_rcxresFilter; v++ ) { // (1)
		 
		 data_pointer->BesselFilteredSignal[v] = new double[m_rcyresFilter];

		 data_pointer->Signal[v] = new double[m_rcyresFilter];

	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for( int v=0; v <= n_poles; v++ ) { // (a)

			for( int f=0; f < 1 ; f++ ) { // (b)

			BesselPolynomial_pointer->ak[v][f] = (double)0.0;

			BesselPolynomial_pointer->bk[v][f] = (double)0.0;

			} //(b)

		 } //(a)
		// (2) initialize (end)

		// (2) initialize (begin)
		for( int v=0; v < m_rcxresFilter; v++ ) { // (a)

			for( int f=0; f < m_rcyresFilter ; f++ ) { // (b)
		 
			data_pointer->BesselFilteredSignal[v][f] = (double)0.0;

			data_pointer->Signal[v][f] = (double)0.0;

			} //(b)

		 } //(a)
		// (2) initialize (end)

	int NofXpixels = n_poles;
	int NofYpixels = 1;
	double savedata = 0;
	int i, j;
	double savedata2 = 0;

	FILE * logfile;
	char logfilename[128]="BesselFilter.log";
		 

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

		cout << "Unable to open log File Now Exit" << endl;
		delete data_pointer;
	    delete BesselPolynomial_pointer;
		exit(0);
	
	} 

		// initialize pointers
		for (i=0; i<=NofXpixels; i++)
		{ 

				BesselPolynomial_pointer->ak[i][0] = (double)1.0;

				BesselPolynomial_pointer->bk[i][0] = (double)1.0;
				
		}

		int k = 0;
		double ck;
		int n;

		for (i=0; i<=NofXpixels; i++)
		{ 
				k++;

				n = (int)NofXpixels*NofYpixels;

				BesselPolynomial_pointer->ak[i][0] *= (double) 2.0 * n - (double)k;

				BesselPolynomial_pointer->bk[i][0] *= (double) k * ((double) n - (double)k);	
	
		}

		k=0;

		for (i=0; i<=NofXpixels; i++)
		{ 
				k++;
	
				n = (int)NofXpixels*NofYpixels;

				ck = (double) pow ( (double) 2.0 , ((double)n - k) ); 
		
				if ( BesselPolynomial_pointer->bk[i][0] != 0.0 && ck != 0.0 )
				
					 BesselPolynomial_pointer->ak[i][0] /= ((double)BesselPolynomial_pointer->bk[i][0] * ck);

				else BesselPolynomial_pointer->ak[i][0] = ((double)1.0);

		}

	FILE * pf;
	char filename[328];

	sprintf(filename, "%s", "Bessel.img");

    fprintf_s(logfile, "%s\t%s\n", "Now Saving Bessel Space in File: ", filename);

    if ((pf = fopen(filename,"wb+"))==NULL)
	{

	 fprintf_s(logfile, "%s\n", "Cannot open file to save Bessel Space Signal");
	 delete data_pointer;
	 delete BesselPolynomial_pointer;
	 exit(0);
	
	} else { // save data


	for (i=0; i<=NofXpixels; i++)
	{ ///save Bessel-Polynomial coefficients data
		for (j=0; j<NofYpixels; j++)
		{

			savedata = (double)BesselPolynomial_pointer->ak[i][j];
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	}///save Bessel-Polynomial coefficients data

	fprintf_s(logfile,"%s\n", "Bessel-Polynomial coefficients Saved");

	fclose (pf);
	} // save data

	
	char FTfilename[128];
	sprintf(FTfilename, "%s", m_FileName);

    fprintf(logfile, "%s\t%s\n", "Now Reading Signal in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 delete data_pointer;
	 delete BesselPolynomial_pointer;

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

	fprintf_s(logfile,"%s\n", "Now Filtering Data");

		double Bessel_Polynomial;	    
		double highPassFilter;
		double pi = 3.141592;
		double x = 0.0, y = 0.0, w = 0.0;
		double deltaT = 1.0;
		double s = 0.0, dx = 0.0, dy = 0.0;

		int XNEI = (int)2;
        int n7 = ( (int)floor( (double)m_rcxresFilter/2.0) );  
	    int n8 = ( (int)floor( (double)m_rcyresFilter/2.0) );
		
		/// Calculate Bessel polynomials (begin)
		for (int i =-n7+XNEI; i < n7-XNEI; i++) {

            for (int j =-n8; j < n8; j++) {

				dx = (i - m_rcxresFilter/2);
				dy = (j - m_rcyresFilter/2);

				s = ((double)sqrt( (double)dx*dx / (m_rcxresFilter*m_rcyresFilter) + 
					               (double)dy*dy / (m_rcxresFilter*m_rcyresFilter) )  ); 

				Bessel_Polynomial = 0.0;
			
				for (int k=0; k<=n_poles; k++)
				{ 

					 // Bessel
						Bessel_Polynomial += (double)BesselPolynomial_pointer->ak[k][0] * 
											 (double)pow ( (double)s , (double)k);
			
				}

				double tf = ((double)Gain/(Bessel_Polynomial/m_cof));

				double encode = (double) exp(-(double)tf) / ( (double) 1.0 + (double) exp(-(double)tf) );

				x = data_pointer->Signal[i + n7][j + n8];
				
				w = data_pointer->Signal[i + n7 - 1][j + n8];

				y = data_pointer->BesselFilteredSignal[i + n7 - 1][j + n8];

				double RC = (double) 1.0 / ((double) 2.0 * pi * encode);

				highPassFilter = ((double) RC) / ((double) deltaT + RC); 
				
				/// Filter (begin)	
				data_pointer->BesselFilteredSignal[i + n7][j + n8] = ((double) highPassFilter * y ) + 
					                                                 (((double)highPassFilter) * ((double)x - w ));
			    /// Filter (ends)
			}
		}//Calculate Bessel polynomials & filter (end)
		

	fprintf(logfile,"%s\n", "Now Saving Filtered Data");
    

	sprintf(filename, "%s%s", "BesselFiltered-", m_FileName);

    fprintf(logfile, "%s\t%s\n", "Now Saving Bessel Filtered Data in File: ", filename);

    if ((pf = fopen(filename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save Bessel Filtered Data");

	 delete data_pointer;
	 delete BesselPolynomial_pointer;

	 exit(0);
	
	} else { // save data


	for (int i=0; i<m_rcxresFilter; i++)
	{ ///save BesselFilteredSignal data
		for (int j=0; j<m_rcyresFilter; j++)
		{

			double savedata2 = (double)data_pointer->BesselFilteredSignal[i][j];
          
            fwrite(&savedata2,sizeof(double),1,pf);

		}
	} ///save BesselFilteredSignal data

	fprintf(logfile,"%s\n", "Bessel Filtered Data Saved");

	fclose (pf);
	} // save data

	delete data_pointer;
	delete BesselPolynomial_pointer;

	fclose(logfile);
	system( logfilename );
	
	} // run the Filtering Process (begin) 

}// end of function