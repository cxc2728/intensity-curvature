#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <process.h>
#include <malloc.h>
#include <conio.h>
#include <string.h>

void OnInverseFourierTransform(char filename1[], char filename2[], char filename3[], char filename4[], int rcyres, int rcxres);

double MAX=-1000000000000000.0;
double MIN=1000000000000000.0;

int main ( int argc, char * argv[] ) {

	char outputFile[128]="INVFT2016.log";

	FILE * savedata;

if (argc < 7) { printf("\n");
				printf("Please type the k-spaceR image file name of the reference\n");
				printf("Please type the k-spaceR image file name to subtract \n");
				printf("Please type the k-spaceI image file name of the reference\n");
				printf("Please type the k-spaceI image file name to subtract \n");
				printf("Please make sure that the image format is in double precision Analyze\n");
				printf("Please enter the number of pixels along the X direction (integer) \n");
				printf("Please enter the number of pixels along the Y direction (integer) \n");
				printf("\n");
				printf("The program will calculate the INVFT of the Difference between the k-spaces\n");
				printf("\n");
				exit(0); }

else { // run the program (begins)

	char imageFileName_R_reference[128];
	char imageFileName_I_reference[128];
	char imageFileName_R_to_subtract[128];
	char imageFileName_I_to_subtract[128];

	sprintf(imageFileName_R_reference, "%s", argv[1]);
	sprintf(imageFileName_R_to_subtract, "%s", argv[2]);
	sprintf(imageFileName_I_reference, "%s", argv[3]);
	sprintf(imageFileName_I_to_subtract, "%s", argv[4]);

	int rcxres = atoi(argv[5]);
	int rcyres = atoi(argv[6]);
	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		printf("Cannot open output file, Now Exit...");

	} else  { // processing (begin)

	printf("\n");
	printf("The image k-spaceR reference file name is: %s \n", imageFileName_R_reference);
	printf("The image k-spaceR file name to subtract is: %s \n", imageFileName_R_to_subtract);
	printf("The image k-spaceI reference file name is: %s \n", imageFileName_I_reference);
	printf("The image k-spaceI file name to subtract is: %s \n", imageFileName_I_to_subtract);
	printf("The number of pixels along the X direction is: %d \n" , atoi(argv[5]));
	printf("The number of pixels along the Y direction is: %d \n" , atoi(argv[6]));
	
	
	fprintf(savedata,"%s%s\n", "The image k-spaceR reference file name is: \n", imageFileName_R_reference);
	fprintf(savedata,"%s%s\n", "The image k-spaceR file name to subtract is: \n", imageFileName_R_to_subtract);
	fprintf(savedata,"%s%s\n", "The image k-spaceI reference file name is: \n", imageFileName_I_reference);
	fprintf(savedata,"%s%s\n", "The image k-spaceI file name to subtract is: \n", imageFileName_I_to_subtract);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: \n", rcxres);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: \n", rcyres);
	
	fprintf(savedata,"\n");

		
		OnInverseFourierTransform(imageFileName_R_reference,
			                      imageFileName_R_to_subtract,
			                      imageFileName_I_reference,
			                      imageFileName_I_to_subtract, rcyres, rcxres);
		

		} // processing begin (end)

	 } // run the program (end)


return 0;
} // main function

void OnInverseFourierTransform(char filename1[], char filename2[], char filename3[], char filename4[], int rcyres, int rcxres)
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

	double * kSpaceR_Subtract = 0;
	double * kSpaceI_Subtract = 0;
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

	if ((kSpaceR_Subtract = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
        free(kSpaceR);
		free(kSpaceI);
		exit(0);

	}

	if ((kSpaceI_Subtract = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);
		free(kSpaceR_Subtract);
		exit(0);

	}

	if ((reconSignal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
	
		fprintf(logfile,"%s\n", "Not enough memory to allocate Imaginary Image data: Exit");
	
		// FIFO memory deallocation method
	    free(kSpaceR);
		free(kSpaceI);
		free(kSpaceR_Subtract);
		free(kSpaceI_Subtract);
		exit(0);

	}

	} // allocate memory

	
	//// read image data and initialize pointers (begins)
    sprintf(imageFilename, "%s", filename1);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofXpixels; i++)
		{ 
			for (j=0; j<NofYpixels; j++)
			{

				index = ((j*NofXpixels)+i);

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
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofXpixels; i++)
		{ 
			for (j=0; j<NofYpixels; j++)
			{

				index = ((j*NofXpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceR_Subtract+index) = (double) number;

			}

		}

		fclose(image);


		for (i=0; i<NofXpixels; i++)
		{ 
			for (j=0; j<NofYpixels; j++)
			{

				index = ((j*NofXpixels)+i);

				*(reconSignal+index) = (double)0.0;
					
			}

		}


	}// read data and initialize pointers

	char imageFilename3[128];

	sprintf(imageFilename3, "%s", filename3);

    if ((image = fopen(imageFilename3,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofXpixels; i++)
		{ 
			for (j=0; j<NofYpixels; j++)
			{

				index = ((j*NofXpixels)+i);

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
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	
	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofXpixels; i++)
		{ 
			for (j=0; j<NofYpixels; j++)
			{

				index = ((j*NofXpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(kSpaceI_Subtract+index) = (double) number;

		
			}

		}

		fclose(image);

	} //// read image data and initialize pointers (ends)

	double real = 0.0, imaginary = 0.0;
	
	///// INV Fourier Transform //////
	for (i=0; i<NofXpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofYpixels; j++)
		{
		
	    	dx = ((int) i - NofXpixels/2);
		    dy = ((int) j - NofYpixels/2);
		
	  	    k2 = ((int)(dx*NofYpixels)+dy);

			w = ((j*NofXpixels)+i);

			real = (double)0.0;
			imaginary = (double)0.0;

			
			for (int s=0; s<NofXpixels; s++)
			{ ///process k-space data

				for (int p=0; p<NofYpixels; p++)
				{ 

					ds = ((int) s - NofXpixels/2);
		            dp = ((int) p - NofYpixels/2);

					k3 = ((int)(dp*NofXpixels)+ds);  
				
					t = ((p*NofXpixels)+s);
					
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

	sprintf(reconFilename, "%s", "RecontructedSignal.img");


    fprintf(logfile, "%s\t%s\n", "Now Saving Reconstructed Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);
	 exit(0);
	
	} else { // save data


	for (i=0; i<NofXpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofYpixels; j++)
		{

			index = ((j*NofXpixels)+i);

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
	 free(kSpaceR_Subtract);
	 free(kSpaceI_Subtract);
	 free(reconSignal);

}