#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include "bmplib.h"
using namespace std;

//============================Add function prototypes here======================

void gaussian(double k[][11], int N, double sigma);
void gaussian_filter(unsigned char output[][SIZE][3], unsigned char input[][SIZE][3], int N, double sigma);
void unsharp(unsigned char output[][SIZE][3], unsigned char input[][SIZE][3], int N, double sigma, double alpha);
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB]);
void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], 
int N, double kernel[][11]);
void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB]);
//============================Do not change code in main()======================



#ifndef AUTOTEST

int main(int argc, char* argv[])
{
   //First check argc
  if(argc < 3)
    {
      //we need at least ./filter <input file> <filter name> to continue
      cout << "usage: ./filter <input file> <filter name> <filter parameters>";
      cout << " <output file name>" << endl;
      return -1;
    }
   //then check to see if we can open the input file
   unsigned char input[SIZE][SIZE][RGB];
   unsigned char output[SIZE][SIZE][RGB];
   char* outfile;
   int N;
   double sigma, alpha;
   //double kernel[11][11];

   // read file contents into input array
   int status = readRGBBMP(argv[1], input); 
   if(status != 0)
   {
      cout << "unable to open " << argv[1] << " for input." << endl;
      return -1;
   }
   //Input file is good, now look at next argument
   if( strcmp("sobel", argv[2]) == 0)
   {
     sobel(output, input);
     outfile = argv[3];
   }
   else if( strcmp("blur", argv[2]) == 0)
   {
     if(argc < 6)
       {
	 cout << "not enough arguments for blur." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     outfile = argv[5];
     gaussian_filter(output, input, N, sigma);
   }
   else if( strcmp("unsharp", argv[2]) == 0)
   {
     if(argc < 7)
       {
	 cout << "not enough arguments for unsharp." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     alpha = atof(argv[5]);
     outfile = argv[6];
     unsharp(output, input, N, sigma, alpha);

   }
   else if( strcmp("dummy", argv[2]) == 0)
   {
     //do dummy
     dummy(output, input);
     outfile = argv[3];
   }
   else
   {
      cout << "unknown filter type." << endl;
      return -1;
   }

   if(writeRGBBMP(outfile, output) != 0)
   {
      cout << "error writing file " << argv[3] << endl;
   }

}   

#endif 

//=========================End Do not change code in main()=====================


// Creates an identity kernel (dummy kernel) that will simply
// copy input to output image and applies it via convolve()
//
// ** This function is complete and need not be changed.
// Use this as an example of how to create a kernel array, fill it in
// appropriately and then use it in a call to convolve.
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   for (int i = 0; i < 3; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         k[i][j] = 0;
      }
   }
   k[1][1] = 1;
   convolve(out, in, 3, k);
}


// Convolves an input image with an NxN kernel to produce the output kernel
// You will need to complete this function by following the 
//  instructions in the comments
void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], 
	      int N, double kernel[][11])
{
 
   int padded[SIZE+10][SIZE+10][RGB];  // Use for input image with appropriate 
                                       // padding
   int temp[SIZE][SIZE][RGB];          // Use for the unclamped output pixel 
                                       // values then copy from temp to out, 
                                       // applying clamping 

   //first set all of padded to 0 (black)
   for(int d = 0; d < SIZE+10; d++){
    for(int e = 0; e < SIZE+10; e++){
      for(int f = 0; f < RGB; f++){
        
        padded[d][e][f] = 0;

      }
    }
   }

 //now copy input into padding to appropriate locations
    for(int g = 0; g < SIZE; g++){ 
      for(int h = 0; h < SIZE; h++){
        for(int o = 0; o < RGB; o++){
        
          padded[g +(N/2)][h + (N /2)][o] = in[g][h][o];

        }
      }
   }



   //initialize temp pixels to 0 (black)
   for(int a = 0; a < SIZE; a++){
    for(int b = 0; b < SIZE; b++){
      for(int c = 0; c < RGB; c++){
     
        temp[a][b][c] = 0; // Initialize values to 0

      }
    }
   }


  //now perform convolve (using convolution equation on each pixel of the 
  // actual image) placing the results in temp (i.e. unclamped result)
  //Here we give you the structure of the convolve for-loops, you need
  //to figure out the loop limits
  for(int y=0 ; y < SIZE; y++){
    for(int x=0; x < SIZE; x++){
      for(int k=0; k < RGB; k++){
         for(int i = -N / 2 ; i <= N / 2; i++){
            for(int j= -N / 2; j<= N / 2; j++){
                temp[y][x][k] += (padded[y + i + (N /2)][x + j +(N / 2)][k])*(kernel[i + (N / 2)][j + (N / 2)]); // Unclamped result produced
   //now clamp and copy to output
   // You may need to cast to avoid warnings from the compiler:
   // (i.e. out[i][j][k] = (unsigned char) temp[i][j][k];)
             }
           }
         }
       }
     }
      
      for(int t = 0; t < SIZE; t++){
        for(int f = 0; f < SIZE; f++){
          for(int e = 0; e < RGB; e++){
            
            if(temp[t][f][e] > 255){ // Clamping
              temp[t][f][e] = 255;
            }

            else if(temp[t][f][e] < 0){ 
              temp[t][f][e] = 0;
            }

            out[t][f][e] = (unsigned char) temp[t][f][e]; // Set value of temp to out
          }
        }
      }

 
}

// You will need to complete this function by following the 
//  instructions in the comments
void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   //double s_h1[3][3] = { {-1, 0, 1}, 
                         //{-2, 0, 2}, 
                         //{-1, 0, 1} };
   //double s_h2[3][3] = { {1, 0, -1}, 
                         //{2, 0, -2}, 
                         //{1, 0, -1} };
   
   unsigned char h1_soble[SIZE][SIZE][RGB]; //hold intemediate images
   unsigned char h2_soble[SIZE][SIZE][RGB]; 

   for (int i = 0; i < 11; i++)
   {
      for(int j=0; j < 11; j++)
      {
         k[i][j] = 0;
      }
   }


   // Copy in 1st 3x3 horizontal sobel kernel (i.e. copy s_h1 into k)
    k[0][0] = -1; // Values copied
    k[0][1] = 0; 
    k[0][2] = 1;
    k[1][0] = -2;
    k[1][1] = 0;
    k[1][2] = 2;
    k[2][0] = -1;
    k[2][1] = 0;
    k[2][2] = 1;

   // Call convolve to apply horizontal sobel kernel with result in h1_soble
    convolve(h1_soble, in, 3, k);


   // Copy in 2nd 3x3 horizontal sobel kernel (i.e. copy s_h2 into k)
    k[0][0] = 1; // Values copied 
    k[0][1] = 0;
    k[0][2] = -1;
    k[1][0] = 2;
    k[1][1] = 0;
    k[1][2] = -2;
    k[2][0] = 1;
    k[2][1] = 0;
    k[2][2] = -1;

   // Call convolve to apply horizontal sobel kernel with result in h2_soble
    convolve(h2_soble, in, 3, k);


   // Add the two results (applying clamping) to produce the final output
    for(int e = 0; e < SIZE; e++){
      for(int f = 0; f < SIZE; f++){
        for(int g = 0; g < RGB; g++){

          if((int)h1_soble[e][f][g] + h2_soble[e][f][g] > 255){ // If the results are greater than 255 Clamp
            out[e][f][g] = 255;
          }

          else if((int)h1_soble[e][f][g] + h2_soble[e][f][g] < 0){
            out[e][f][g] = 0;
          }

          else{
            out[e][f][g] = h1_soble[e][f][g] + h2_soble[e][f][g];
          }

          }
        }
      }
    }




// Add the rest of your functions here (i.e. gaussian, gaussian_filter, unsharp)

void gaussian(double k[][11], int N, double sigma){
  double total = 0.0;
 

  for(int d = 0; d < N; d++){
    for(int e = 0; e < N; e++){

      //Formula Integration
      k[d][e] = (double) exp(-((pow(d - N / 2, 2) / (2 * pow(sigma, 2))) + (pow(e - N / 2, 2) / (2 *pow(sigma,2)))));
      // Total up kernel values
      total+= k[d][e];
    }
  }

  for(int w = 0; w < N; w++){
    for(int x = 0; x < N; x++){
      k[w][x] /= total;
    }
  }

  for(int s = 0; s < N; s++){
    for(int v = 0; v < N; v++){
      cout << fixed << setprecision(4) << k[s][v] << " "; // Output values - Fix the values to 4
    }
    cout << endl;
  }

}

void gaussian_filter(unsigned char output[][SIZE][3], unsigned char input[][SIZE][3], int N, double sigma){

  double kernel[11][11];
  gaussian(kernel, N, sigma);
  convolve(output, input, N, kernel);

}

void unsharp(unsigned char output[][SIZE][3], unsigned char input[][SIZE][3], int N, double sigma, double alpha){

  unsigned char distort[SIZE][SIZE][3]; // Variable decleration
  gaussian_filter(distort, input, N, sigma);
  double detailed_image[SIZE][SIZE][3];

  for(int a = 0; a < SIZE; a++){
    for(int b = 0; b < SIZE; b++){
      for(int c = 0; c < 3; c++){
        detailed_image[a][b][c] = (double)input[a][b][c] - (double)distort[a][b][c]; // Formula implementation 
        //Perform the subtraction operation for each pixel
      }
    }
  }

  for(int r = 0; r < SIZE; r++){
    for(int t = 0; t < SIZE; t++){
      for(int s = 0; s < 3; s++){
        double val = (double)input[r][t][s] + alpha * detailed_image[r][t][s]; // Formula implementation
        if(val > 255){ // Clamp values
          val = 255;
        }
        else if(val < 0){
          val = 0;
        } 
        output[r][t][s] = (unsigned char)val;
      }
    }
  }


}


