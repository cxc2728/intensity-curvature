#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <process.h>
#include <malloc.h>
#include <conio.h>
#include <string.h>

void On_Intensity_Curvature_Process(char filename1[], char filename2[], char filename3[], char filename4[], char filename5[], int rcyres, int rcxres, double Ncc, double Nicf, double Nsri, double Nrc);

int main ( int argc, char * argv[] ) {

	char outputFile[128]="IC2016.log";

	FILE * savedata;

if (argc < 12) { printf("\n");
				printf("Please type the MRI image file name \n");
				printf("Please type the Classic-Curvature (CC) image file name \n");
				printf("Please type the Intensity-Curvature Functional (ICF) image file name \n");
				printf("Please type the Signal Resilient to Interpolation (SRI) image file name \n");
				printf("Please type the Resilient-Curvature (RC) image file name \n");
				printf("Please make sure that the images format is in double precision Analyze\n");
				printf("Please enter the number of pixels along the X direction (integer) \n");
				printf("Please enter the number of pixels along the Y direction (integer) \n");
				printf("Please enter the number of multiplications of the CC image (double) \n");
				printf("Please enter the number of multiplications of the ICF image (double) \n");
				printf("Please enter the number of multiplications of the SRI image (double) \n");
				printf("Please enter the number of multiplications of the RC image (double) \n");
				printf("\n");
				printf("The program will marge the MRI with the Intensity-Curvature images\n");
				printf("\n");
				exit(0); }

else { // run the program (begins)

	char imageFileName_MRI[128];
	char imageFileName_Classic_Curvature[128];
	char imageFileName_Intensity_Curvature_Functional[128];
	char imageFileName_Signal_Resilient_to_Interpolation[128];
	char imageFileName_Resilient_Curvature[128];

	
	sprintf(imageFileName_MRI, "%s", argv[1]);
	sprintf(imageFileName_Classic_Curvature, "%s", argv[2]);
	sprintf(imageFileName_Intensity_Curvature_Functional, "%s", argv[3]);
	sprintf(imageFileName_Signal_Resilient_to_Interpolation, "%s", argv[4]);
	sprintf(imageFileName_Resilient_Curvature, "%s", argv[5]);


	int rcxres = atoi(argv[6]);
	int rcyres = atoi(argv[7]);

	double Ncc = atof(argv[8]);
	double Nicf = atof(argv[9]);
	double Nsri = atof(argv[10]);
	double Nrc = atof(argv[11]);
	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		printf("Cannot open output file, Now Exit...");

	} else  { // processing (begin)

	printf("\n");
	printf("The MRI file name is: %s \n", imageFileName_MRI);
	printf("The Classic-Curvature file name is: %s \n", imageFileName_Classic_Curvature);
	printf("The Intensity-Curvature Functional file name is: %s \n", imageFileName_Intensity_Curvature_Functional);
	printf("The Signal Resilient to Interpolation file name is: %s \n", imageFileName_Signal_Resilient_to_Interpolation);
	printf("The Resilient-Curvature file name is: %s \n", imageFileName_Resilient_Curvature);
	printf("The number of pixels along the X direction is: %d \n" , atoi(argv[6]));
	printf("The number of pixels along the Y direction is: %d \n" , atoi(argv[7]));
	printf("The number of multiplications of the CC image is: %lf \n" , atof(argv[8]));
	printf("The number of multiplications of the ICF image is: %lf \n" , atof(argv[9]));
	printf("The number of multiplications of the SRI image is: %lf \n" , atof(argv[10]));
	printf("The number of multiplications of the RC image is: %lf \n" , atof(argv[11]));
	
	
	fprintf(savedata,"%s%s\n", "The MRI file name is: ", imageFileName_MRI);
	fprintf(savedata,"%s%s\n", "The Classic-Curvature file name is: ", imageFileName_Classic_Curvature);
	fprintf(savedata,"%s%s\n", "Intensity-Curvature Functional file name is: ", imageFileName_Intensity_Curvature_Functional);
	fprintf(savedata,"%s%s\n", "Signal Resilient to Interpolation file name is: ", imageFileName_Signal_Resilient_to_Interpolation);
	fprintf(savedata,"%s%s\n", "Resilient-Curvature file name is: ", imageFileName_Resilient_Curvature);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", rcxres);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", rcyres);
	fprintf(savedata,"%s%lf\n", "The number of multiplications of the CC image is: ", Ncc);
	fprintf(savedata,"%s%lf\n", "The number of multiplications of the ICF image is: ", Nicf);
	fprintf(savedata,"%s%lf\n", "The number of multiplications of the SRI image is: ", Nsri);
	fprintf(savedata,"%s%lf\n", "The number of multiplications of the RC image is: ", Nrc);
	
	fprintf(savedata,"\n");

	
		On_Intensity_Curvature_Process(imageFileName_MRI, imageFileName_Classic_Curvature,
			                           imageFileName_Intensity_Curvature_Functional, 
									   imageFileName_Signal_Resilient_to_Interpolation, 
									   imageFileName_Resilient_Curvature, rcyres, rcxres, Ncc, Nicf, Nsri, Nrc);

		} // processing begin (end)

	 } // run the program (end)


return 0;
} // main function

void On_Intensity_Curvature_Process(char filename1[], char filename2[], char filename3[], char filename4[], char filename5[], int rcyres, int rcxres, double Ncc, double Nicf, double Nsri, double Nrc)
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;
	
	int i, j, index;
	
	FILE * logfile;
	char logfilename[128]="IC2016Process.log";

	FILE *image;
	char imageFilename[256];

	double * MRI = 0;
	double * CC = 0;
	double * ICF = 0;
	double * SRI = 0;
	double * RC = 0;
	double * reconSignal = 0;

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	 exit(0);
	
	} else { // allocate memory

		
	printf("%s\n", "Now Intensity-Curvature Processing...");
    fprintf(logfile,"%s\n", "Now Intensity-Curvature Processing...");

	if ((MRI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Image data: Exit");
   
		exit(0);

	}

	if ((CC = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to Image data: Exit");
   
		// FIFO memory deallocation method
		free(MRI);
		exit(0);

	}

	if ((ICF = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to Image data: Exit");
   
		// FIFO memory deallocation method
		free(MRI);
		free(CC);
		exit(0);

	}

	if ((SRI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to Image data: Exit");
   
		// FIFO memory deallocation method
		free(MRI);
		free(CC);
		free(ICF);
		exit(0);

	}

	if ((RC = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to Image data: Exit");
   
		// FIFO memory deallocation method
		free(MRI);
		free(CC);
		free(ICF);
		free(SRI);
		exit(0);

	}

	if ((reconSignal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
	
		fprintf(logfile,"%s\n", "Not enough memory to allocate Image data: Exit");
	
		// FIFO memory deallocation method
	    free(MRI);
		free(CC);
		free(ICF);
		free(SRI);
		free(RC);
		exit(0);

	}

	} // allocate memory

	
	//// read image data and initialize pointers (begins)
    sprintf(imageFilename, "%s", filename1);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	  free(MRI);
	  free(CC);
	  free(ICF);
	  free(SRI);
	  free(RC);
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
				
				*(MRI+index) = (double) number;

		
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
	free(MRI);
	free(CC);
	free(ICF);
	free(SRI);
	free(RC);
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
				
				*(CC+index) = (double) number;

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
	
	sprintf(imageFilename, "%s", filename3);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	  free(MRI);
	  free(CC);
	  free(ICF);
	  free(SRI);
	  free(RC);
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
				
				*(ICF+index) = (double) number;

		
			}

		}

		fclose(image);

	}// read data and initialize pointers


	sprintf(imageFilename, "%s", filename4);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	  free(MRI);
	  free(CC);
	  free(ICF);
	  free(SRI);
	  free(RC);
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
				
				*(SRI+index) = (double) number;

		
			}

		}

		fclose(image);

	}// read data and initialize pointers

	sprintf(imageFilename, "%s", filename5);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	  free(MRI);
	  free(CC);
	  free(ICF);
	  free(SRI);
	  free(RC);
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
				
				*(RC+index) = (double) number;

		
			}

		}

		fclose(image);

	}// read data and initialize pointers

	///process begins
	for (i=0; i<NofXpixels; i++)
	{ 

		for (j=0; j<NofYpixels; j++)
		{
		
	    	index = ((j*NofXpixels)+i);

	
			*(reconSignal + index) = ( (double)*(MRI+index) + (double)pow( (double)*(CC+index) , (double)Ncc ) 
				
															+ (double)pow( (double)*(ICF+index) , (double)Nicf )  
															
															+ (double)pow( (double)*(SRI+index) , (double)Nsri ) 
															
															+ (double)pow( (double)*(RC+index) , (double)Nrc ) );

			} 
			
	}///process ends

	double savedata = 0.0;
	FILE * pf;
	char reconFilename[128];

	sprintf(reconFilename, "%s", "Intensity-Curvature.img");


    fprintf(logfile, "%s\t%s\n", "Now Saving Intensity-Curvature Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save Signal");

	  // FIFO memory deallocation method
	  free(MRI);
	  free(CC);
	  free(ICF);
	  free(SRI);
	  free(RC);
	  free(reconSignal);
	  exit(0);
	
	} else { // save data


	for (i=0; i<NofXpixels; i++)
	{ ///save IC data
		for (j=0; j<NofYpixels; j++)
		{

			index = ((j*NofXpixels)+i);

			savedata = (double)*(reconSignal+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save IC data

	fprintf(logfile,"%s\n", "Intensity-Curvature Signal Saved");

	fclose (pf);
	} // save data

	fclose(logfile);
	system( logfilename );
		
	  // FIFO memory deallocation method
	  free(MRI);
	  free(CC);
	  free(ICF);
	  free(SRI);
	  free(RC);
	  free(reconSignal);

}