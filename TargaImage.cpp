///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if


    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }

    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
	
	for(int i=0;i<height*width*4; i+=4)
	{
		unsigned char rgb1[3];
		RGBA_To_RGB(data + i, rgb1);
		data[i]		=0.299*rgb1[0] + 0.587*rgb1[1]+0.114*rgb1[2];
		data[i+1]	=0.299*rgb1[0] + 0.587*rgb1[1]+0.114*rgb1[2];
		data[i+2]	=0.299*rgb1[0] + 0.587*rgb1[1]+0.114*rgb1[2];
	}

	//ClearToBlack();
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
	for(int i=0;i<height*width*4; i+=4)
	{
		unsigned char rgb1[3];
		RGBA_To_RGB(data + i, rgb1);
		
		data[i]		=floor((double)rgb1[0]/32)*32;
		data[i+1]	=floor((double)rgb1[1]/32)*32;
		data[i+2]	=floor((double)rgb1[2]/64)*64;
	}

    //ClearToBlack();
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
	//cout<<"This algorithm may take upto 5-10 minutes to run.But gives a perfect output.Please be patient.Thanks"<<endl;
	int* temp_data=new int[height*width*4];
	for(int i=0;i<height*width*4; i+=4)
	{
		unsigned char rgb1[3];
		RGBA_To_RGB(data + i, rgb1);
		
		temp_data[i]	=floor((double)rgb1[0]/8)*8;
		temp_data[i+1]	=floor((double)rgb1[1]/8)*8;
		temp_data[i+2]	=floor((double)rgb1[2]/8)*8;
		temp_data[i+3]=1;
	}
	
	for(int i=0;i<height*width*4;i+=4)
	{

		for(int j=i+4;j<height*width*4;j+=4)
		{
			if((temp_data[i]==temp_data[j]) && (temp_data[i+1]==temp_data[j+1]) && (temp_data[i+2]==temp_data[j+2]) && (temp_data[i+3]>0))
			{
				temp_data[i+3]=temp_data[i+3]+1;
				temp_data[j+3]=0;
				
			}
		}

	}
	
	int* temp=new int[4];
    for(int i=0;i<height*width*4;i+=4)
    {
        for(int j=0;j<i;j+=4)
        {
            if(temp_data[i+3]>temp_data[j+3])
            {
                temp[0]=temp_data[i];
				temp[1]=temp_data[i+1];
				temp[2]=temp_data[i+2];
				temp[3]=temp_data[i+3];
				//swap 
                temp_data[i]=temp_data[j];
				temp_data[i+1]=temp_data[j+1];
				temp_data[i+2]=temp_data[j+2];
				temp_data[i+3]=temp_data[j+3];

				temp_data[j]=temp[0];
                temp_data[j+1]=temp[1];
				temp_data[j+2]=temp[2];
				temp_data[j+3]=temp[3];
            }

        }

    }
/*
	for(int i=0;i<10*4;i+=4)
	{
		cout<<temp_data[i]<<","<<temp_data[i+1]<<","<<temp_data[i+2]<<","<<temp_data[i+3]<<endl;
	}
	*/

	int* pop_col=new int[256*4];
	for(int i=0;i<256*4;i+=4)
	{
		pop_col[i]=temp_data[i];
		pop_col[i+1]=temp_data[i+1];
		pop_col[i+2]=temp_data[i+2];
		pop_col[i+3]=temp_data[i+3];
		//cout<<pop_col[i]<<","<<pop_col[i+1]<<","<<pop_col[i+2]<<","<<pop_col[i+3]<<endl;
	}

	int* distance=new int[height*width];
	for(int i=0;i<height*width*4;i+=4)
	{
		int k=0;
		for(int j=0;j<256*4;j+=4)
		{
			distance[k]=(data[i]-pop_col[j])*(data[i]-pop_col[j])+(data[i+1]-pop_col[j+1])*(data[i+1]-pop_col[j+1])+(data[i+2]-pop_col[j+2])*(data[i+2]-pop_col[j+2]);
			k=k+1;
		}
		int x = 0;
		int min = distance[x];
		int pos=x;
		while ( x < 256 )
		{
			if ( distance[x] < min )
			{
                min = distance[x];
				pos=x;
			}
			x++;
		}
		data[i]=pop_col[pos*4];
		data[i+1]=pop_col[pos*4+1];
		data[i+2]=pop_col[pos*4+2];
	}


	delete[] temp_data;
	delete[] temp;
	delete[] pop_col;
	delete[] distance;
    //ClearToBlack();
    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
	To_Grayscale();
	for(int i=0;i<height*width*4;i+=4)
	{
		if(data[i] < 128)
		{
			data[i]=0;
			data[i+1]=0;
			data[i+2]=0;
		}
		else
		{
			data[i]=255;
			data[i+1]=255;
			data[i+2]=255;
		}
	}

    //ClearToBlack();
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
	double* random_number=new double[height*width];
	for(int i =0;i<height*width;i++)
	{
		random_number[i]=rand()%41;
		random_number[i]=random_number[i]-20;
		random_number[i]=random_number[i]/100;
		
	}

	To_Grayscale();
	for(int i=0;i<height*width*4;i+=4)
	{
		
		data[i]=data[i]+random_number[i/4]*255;
		
		if(data[i] < 128)
		{
			data[i]=0;
			data[i+1]=0;
			data[i+2]=0;
		}
		else
		{
			data[i]=255;
			data[i+1]=255;
			data[i+2]=255;
		}
	}

	delete[] random_number;
    //ClearToBlack();
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{

	int error=0;
	int new_val=0;
	int k,i,j;
	
	
	To_Grayscale();
	for(j=0;j<(height/2)-1;j++)
	{
		for(i=2*width*4*j+4;i<width*4-4+2*width*4*j;i+=4)
		{
		
			
			if(data[i]<128)
			new_val=0;
			else
			new_val=255;

			error= (int)data[i]-new_val;

			data[i]=new_val;
			data[i+1]=new_val;
			data[i+2]=new_val;

			if(((16*data[i+4]+7*error)/16) <0)
				data[i+4]=0;
			else if(((16*data[i+4]+7*error)/16)>255)
				data[i+4]=255;
			else
				data[i+4]=(16*data[i+4]+7*error)/16;

			if(((16*data[i+4+1]+7*error)/16) <0)
				data[i+1+4]=0;
			else if(((16*data[i+4+1]+7*error)/16)>255)
				data[i+1+4]=255;
			else
				data[i+1+4]=(16*data[i+4+1]+7*error)/16;

			if(((16*data[i+4+2]+7*error)/16) <0)
				data[i+2+4]=0;
			else if(((16*data[i+4+2]+7*error)/16)>255)
				data[i+2+4]=255;
			else
				data[i+2+4]=(16*data[i+4+2]+7*error)/16;

			

			if(((16*data[i+width*4+4]+1*error)/16)<0)
				data[i+width*4+4]=0;
			else if(((16*data[i+width*4+4]+1*error)/16)>255)
				data[i+width*4+4]=255;
			else
				data[i+width*4+4]=(16*data[i+width*4+4]+1*error)/16;

			if(((16*data[i+width*4+4+1]+1*error)/16)<0)
				data[i+1+width*4+4]=0;
			else if(((16*data[i+width*4+4+1]+1*error)/16)>255)
				data[i+1+width*4+4]=255;
			else
				data[i+1+width*4+4]=(16*data[i+width*4+4+1]+1*error)/16;

			if(((16*data[i+width*4+4+2]+1*error)/16)<0)
				data[i+2+width*4+4]=0;
			else if(((16*data[i+width*4+4+2]+1*error)/16)>255)
				data[i+2+width*4+4]=255;
			else
				data[i+2+width*4+4]=(16*data[i+width*4+4+2]+1*error)/16;



			if(((16*data[i+width*4]+5*error)/16)<0)
				data[i+width*4]=0;
			else if(((16*data[i+width*4]+5*error)/16)>255)
				data[i+width*4]=255;
			else
				data[i+width*4]=(16*data[i+width*4]+5*error)/16;

			
			if(((16*data[i+width*4+1]+5*error)/16)<0)
				data[i+1+width*4]=0;
			else if(((16*data[i+width*4+1]+5*error)/16)>255)
				data[i+1+width*4]=255;
			else
				data[i+1+width*4]=(16*data[i+width*4+1]+5*error)/16;

			
			if(((16*data[i+width*4+2]+5*error)/16)<0)
				data[i+2+width*4]=0;
			else if(((16*data[i+width*4+2]+5*error)/16)>255)
				data[i+2+width*4]=255;
			else
				data[i+2+width*4]=(16*data[i+width*4+2]+5*error)/16;



		
			if(((16*data[i+width*4-4]+3*error)/16)<0)
				data[i+width*4-4]=0;
			else if(((16*data[i+width*4-4]+3*error)/16)>255)
				data[i+width*4-4]=255;
			else
				data[i+width*4-4]=(16*data[i+width*4-4]+3*error)/16;

			if(((16*data[i+width*4-4+1]+3*error)/16)<0)
				data[i+1+width*4-4]=0;
			else if(((16*data[i+width*4-4+1]+3*error)/16)>255)
				data[i+1+width*4-4]=255;
			else
				data[i+1+width*4-4]=(16*data[i+width*4-4+1]+3*error)/16;

			if(((16*data[i+width*4-4+2]+3*error)/16)<0)
				data[i+2+width*4-4]=0;
			else if(((16*data[i+width*4-4+2]+3*error)/16)>255)
				data[i+2+width*4-4]=255;
			else
				data[i+2+width*4-4]=(16*data[i+width*4-4+2]+3*error)/16;

			
		}

		for(int k=2*width*4-8+2*width*4*j;k>width*4+2*width*4*j;k-=4)
		{
			
			if(data[k]<128)
			new_val=0;
			else
			new_val=255;

			error= data[k]-new_val;

			data[k]=new_val;
			data[k+1]=new_val;
			data[k+2]=new_val;

			if(((16*data[k-4]+7*error)/16)<0)
				data[k-4]=0;
			else if(((16*data[k-4]+7*error)/16)>255)
				data[k-4]=255;
			else
				data[k-4]=(16*data[k-4]+7*error)/16;

			if(((16*data[k-4+1]+7*error)/16)<0)
				data[k+1-4]=0;
			else if(((16*data[k-4+1]+7*error)/16)>255)
				data[k+1-4]=255;
			else
				data[k+1-4]=(16*data[k-4+1]+7*error)/16;

			if(((16*data[k-4+2]+7*error)/16)<0)
				data[k+2-4]=0;
			else if(((16*data[k-4+2]+7*error)/16)>255)
				data[k+2-4]=255;
			else
				data[k+2-4]=(16*data[k-4+2]+7*error)/16;




			if(((16*data[k+width*4-4]+1*error)/16)<0)
				data[k+width*4-4]=0;
			else if(((16*data[k+width*4-4]+1*error)/16)>255)
				data[k+width*4-4]=255;
			else
				data[k+width*4-4]=(16*data[k+width*4-4]+1*error)/16;

			if(((16*data[k+width*4-4+1]+1*error)/16)<0)
				data[k+1+width*4-4]=0;
			else if(((16*data[k+width*4-4+1]+1*error)/16)>255)
				data[k+1+width*4-4]=255;
			else
				data[k+1+width*4-4]=(16*data[k+width*4-4+1]+1*error)/16;

			if(((16*data[k+width*4-4+2]+1*error)/16)<0)
				data[k+2+width*4-4]=0;
			else if(((16*data[k+width*4-4+2]+1*error)/16)>255)
				data[k+2+width*4-4]=255;
			else
				data[k+2+width*4-4]=(16*data[k+width*4-4+2]+1*error)/16;


	
			if(((16*data[k+width*4]+5*error)/16)<0)
				data[k+width*4]=0;
			else if(((16*data[k+width*4]+5*error)/16)>255)
				data[k+width*4]=255;
			else
				data[k+width*4]=(16*data[k+width*4]+5*error)/16;

			if(((16*data[k+width*4+1]+5*error)/16)<0)
				data[k+1+width*4]=0;
			else if(((16*data[k+width*4+1]+5*error)/16)>255)
				data[k+1+width*4]=255;
			else
				data[k+1+width*4]=(16*data[k+width*4+1]+5*error)/16;

			if(((16*data[k+width*4+2]+5*error)/16)<0)
				data[k+2+width*4]=0;
			else if(((16*data[k+width*4+2]+5*error)/16)>255)
				data[k+2+width*4]=255;
			else
				data[k+2+width*4]=(16*data[k+width*4+2]+5*error)/16;



		
			if(((16*data[k+width*4+4]+3*error)/16)<0)
				data[k+width*4+4]=0;
			else if(((16*data[k+width*4+4]+3*error)/16)>255)
				data[k+width*4+4]=255;
			else
				data[k+width*4+4]=(16*data[k+width*4+4]+3*error)/16;

			if(((16*data[k+width*4+4+1]+3*error)/16)<0)
				data[k+1+width*4+4]=0;
			else if(((16*data[k+width*4+4+1]+3*error)/16)>255)
				data[k+1+width*4+4]=255;
			else
				data[k+1+width*4+4]=(16*data[k+width*4+4+1]+3*error)/16;

			if(((16*data[k+width*4+4+2]+3*error)/16)<0)
				data[k+2+width*4+4]=0;
			else if(((16*data[k+width*4+4+2]+3*error)/16)>255)
				data[k+2+width*4+4]=255;
			else
				data[k+2+width*4+4]=(16*data[k+width*4+4+2]+3*error)/16;


			
		}
		
	}
	

    //ClearToBlack();
    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
	To_Grayscale();
	double sum_intensity=0;
	for(int i=0;i<height*width*4;i+=4)
	{
		sum_intensity=sum_intensity + data[i];
	}
	double avg_intensity = (sum_intensity)/(height*width);
	avg_intensity=avg_intensity/255;
	cout<<avg_intensity<<endl;
	
	
	int* temp=new int[width * height];
	for(int i =0;i<height*width*4;i+=4)
	{
		temp[i/4]=(int)data[i];
	}
	
	//sorting
	int gap , item , i , j ;
	for ( gap = (( height*width)-1 ) / 2 ; gap>0; gap /= 2 )
	{
		for ( i = gap ;i < height*width;i ++ )
		{
			item = temp [i ];
			j = i-gap ;
			while ( j >=0 && item < temp [j ] )
			{
				temp [ j + gap ] = temp [j ] ;
				j = j - gap ;
			}
			temp [j+ gap ] = item ;
		}
	}

	 int loc_threshold = (1- avg_intensity) * height * width;
	 int threshold = temp[loc_threshold];
	

	 for(int i=0;i<height*width*4;i+=4)
	{
		if(data[i] < threshold)
		{
			data[i]=0;
			data[i+1]=0;
			data[i+2]=0;
		}
		else
		{
			data[i]=255;
			data[i+1]=255;
			data[i+2]=255;
		}
	}
	
	 delete[] temp;
    //ClearToBlack();
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
	To_Grayscale();
	float t_m[4][4]={{0.7500*255,0.3750*255,0.6250*255,0.2500*255},
	{0.0625*255,1.0000*255,0.8750*255,0.4375*255},{0.5000*255,0.8125*255,0.9375*255,0.1250*255},
	{0.1875*255,0.5625*255,0.3125*255,0.6875*255}};
	for(int j=0;j<height/4;j++)
	{

	for(int i=0+4*j*width*4;i<width*4+4*j*width*4;i+=4)
	{
		if(data[i]<t_m[0][(i/4)%4])
		{
			data[i]=0;
			data[i+1]=0;
			data[i+2]=0;
		}
		else
		{
			data[i]=255;
			data[i+1]=255;
			data[i+2]=255;
		}
	}
	for(int i=width*4+4*j*width*4;i<2*width*4+4*j*width*4;i+=4)
	{
		if(data[i]<t_m[1][(i/4)%4])
		{
			data[i]=0;
			data[i+1]=0;
			data[i+2]=0;
		}
		else
		{
			data[i]=255;
			data[i+1]=255;
			data[i+2]=255;
		}
	}
	for(int i=2*width*4+4*j*width*4;i<3*width*4+4*j*width*4;i+=4)
	{
		if(data[i]<t_m[2][(i/4)%4])
		{
			data[i]=0;
			data[i+1]=0;
			data[i+2]=0;
		}
		else
		{
			data[i]=255;
			data[i+1]=255;
			data[i+2]=255;
		}
	}
	for(int i=3*width*4+4*j*width*4;i<4*width*4+4*j*width*4;i+=4)
	{
		if(data[i]<t_m[3][(i/4)%4])
		{
			data[i]=0;
			data[i+1]=0;
			data[i+2]=0;
		}
		else
		{
			data[i]=255;
			data[i+1]=255;
			data[i+2]=255;
		}
	}
	}
    //ClearToBlack();
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
	int error=0;
	int new_val=0;
	int k,i,j,av;
	
	for(j=0;j<(height/2)-1;j++)
	{
		for(i=2*width*4*j+4;i<width*4-4+2*width*4*j;i+=4)
		{
			//new_val=ceil((double)data[i]/32)*32;
			new_val=((floor(data[i]/32.0) +floor(data[i]/32.0) + 1)/2)*32;
			error= (int)data[i]-new_val;
			data[i]=new_val;

			if(((16*data[i+4]+7*error)/16) <0)
				data[i+4]=0;
			else if(((16*data[i+4]+7*error)/16)>255)
				data[i+4]=255;
			else
				data[i+4]=(16*data[i+4]+7*error)/16;

			if(((16*data[i+width*4+4]+1*error)/16)<0)
				data[i+width*4+4]=0;
			else if(((16*data[i+width*4+4]+1*error)/16)>255)
				data[i+width*4+4]=255;
			else
				data[i+width*4+4]=(16*data[i+width*4+4]+1*error)/16;

			if(((16*data[i+width*4]+5*error)/16)<0)
				data[i+width*4]=0;
			else if(((16*data[i+width*4]+5*error)/16)>255)
				data[i+width*4]=255;
			else
				data[i+width*4]=(16*data[i+width*4]+5*error)/16;

			if(((16*data[i+width*4-4]+3*error)/16)<0)
				data[i+width*4-4]=0;
			else if(((16*data[i+width*4-4]+3*error)/16)>255)
				data[i+width*4-4]=255;
			else
				data[i+width*4-4]=(16*data[i+width*4-4]+3*error)/16;


			
			//new_val=ceil((double)data[i+1]/32)*32;
			new_val=((floor(data[i+1]/32.0) +floor(data[i+1]/32.0) + 1)/2)*32;
			error= (int)data[i+1]-new_val;

			data[i+1]=new_val;
			
			if(((16*data[i+4+1]+7*error)/16) <0)
				data[i+1+4]=0;
			else if(((16*data[i+4+1]+7*error)/16)>255)
				data[i+1+4]=255;
			else
				data[i+1+4]=(16*data[i+4+1]+7*error)/16;

			if(((16*data[i+width*4+4+1]+1*error)/16)<0)
				data[i+1+width*4+4]=0;
			else if(((16*data[i+width*4+4+1]+1*error)/16)>255)
				data[i+1+width*4+4]=255;
			else
				data[i+1+width*4+4]=(16*data[i+width*4+4+1]+1*error)/16;

				if(((16*data[i+width*4+1]+5*error)/16)<0)
				data[i+1+width*4]=0;
			else if(((16*data[i+width*4+1]+5*error)/16)>255)
				data[i+1+width*4]=255;
			else
				data[i+1+width*4]=(16*data[i+width*4+1]+5*error)/16;

				if(((16*data[i+width*4-4+1]+3*error)/16)<0)
				data[i+1+width*4-4]=0;
			else if(((16*data[i+width*4-4+1]+3*error)/16)>255)
				data[i+1+width*4-4]=255;
			else
				data[i+1+width*4-4]=(16*data[i+width*4-4+1]+3*error)/16;



				
			//new_val=ceil((double)data[i+2]/64)*64;
			new_val=((floor(data[i+2]/64.0) +floor(data[i+2]/64.0) + 1)/2)*64;
			error= (int)data[i+2]-new_val;

			data[i+2]=new_val;

			if(((16*data[i+4+2]+7*error)/16) <0)
				data[i+2+4]=0;
			else if(((16*data[i+4+2]+7*error)/16)>255)
				data[i+2+4]=255;
			else
				data[i+2+4]=(16*data[i+4+2]+7*error)/16;

			if(((16*data[i+width*4+4+2]+1*error)/16)<0)
				data[i+2+width*4+4]=0;
			else if(((16*data[i+width*4+4+2]+1*error)/16)>255)
				data[i+2+width*4+4]=255;
			else
				data[i+2+width*4+4]=(16*data[i+width*4+4+2]+1*error)/16;

			
			if(((16*data[i+width*4+2]+5*error)/16)<0)
				data[i+2+width*4]=0;
			else if(((16*data[i+width*4+2]+5*error)/16)>255)
				data[i+2+width*4]=255;
			else
				data[i+2+width*4]=(16*data[i+width*4+2]+5*error)/16;


			if(((16*data[i+width*4-4+2]+3*error)/16)<0)
				data[i+2+width*4-4]=0;
			else if(((16*data[i+width*4-4+2]+3*error)/16)>255)
				data[i+2+width*4-4]=255;
			else
				data[i+2+width*4-4]=(16*data[i+width*4-4+2]+3*error)/16;
			
		}

		for(int k=2*width*4-8+2*width*4*j;k>width*4+2*width*4*j;k-=4)
		{
			
			//new_val=ceil((double)data[k]/32)*32;
			new_val=((floor(data[k]/32.0) +floor(data[k]/32.0) + 1)/2)*32;
			error= data[k]-new_val;

			data[k]=new_val;

			if(((16*data[k-4]+7*error)/16)<0)
				data[k-4]=0;
			else if(((16*data[k-4]+7*error)/16)>255)
				data[k-4]=255;
			else
				data[k-4]=(16*data[k-4]+7*error)/16;

			if(((16*data[k+width*4-4]+1*error)/16)<0)
				data[k+width*4-4]=0;
			else if(((16*data[k+width*4-4]+1*error)/16)>255)
				data[k+width*4-4]=255;
			else
				data[k+width*4-4]=(16*data[k+width*4-4]+1*error)/16;

			if(((16*data[k+width*4]+5*error)/16)<0)
				data[k+width*4]=0;
			else if(((16*data[k+width*4]+5*error)/16)>255)
				data[k+width*4]=255;
			else
				data[k+width*4]=(16*data[k+width*4]+5*error)/16;

			if(((16*data[k+width*4+4]+3*error)/16)<0)
				data[k+width*4+4]=0;
			else if(((16*data[k+width*4+4]+3*error)/16)>255)
				data[k+width*4+4]=255;
			else
				data[k+width*4+4]=(16*data[k+width*4+4]+3*error)/16;


			//new_val=ceil((double)data[k+1]/32)*32;
			new_val=((floor(data[k+1]/32.0) +floor(data[k+1]/32.0) + 1)/2)*32;
			error= data[k+1]-new_val;

			data[k+1]=new_val;
			
			if(((16*data[k-4+1]+7*error)/16)<0)
				data[k+1-4]=0;
			else if(((16*data[k-4+1]+7*error)/16)>255)
				data[k+1-4]=255;
			else
				data[k+1-4]=(16*data[k-4+1]+7*error)/16;

			if(((16*data[k+width*4-4+1]+1*error)/16)<0)
				data[k+1+width*4-4]=0;
			else if(((16*data[k+width*4-4+1]+1*error)/16)>255)
				data[k+1+width*4-4]=255;
			else
				data[k+1+width*4-4]=(16*data[k+width*4-4+1]+1*error)/16;

			
			if(((16*data[k+width*4+1]+5*error)/16)<0)
				data[k+1+width*4]=0;
			else if(((16*data[k+width*4+1]+5*error)/16)>255)
				data[k+1+width*4]=255;
			else
				data[k+1+width*4]=(16*data[k+width*4+1]+5*error)/16;

			
			if(((16*data[k+width*4+4+1]+3*error)/16)<0)
				data[k+1+width*4+4]=0;
			else if(((16*data[k+width*4+4+1]+3*error)/16)>255)
				data[k+1+width*4+4]=255;
			else
				data[k+1+width*4+4]=(16*data[k+width*4+4+1]+3*error)/16;

			//new_val=ceil((double)data[k+2]/64)*64;
			new_val=((floor(data[k+2]/64.0) +floor(data[k+2]/64.0) + 1)/2)*64;
			error= data[k+2]-new_val;

			data[k+2]=new_val;

			if(((16*data[k-4+2]+7*error)/16)<0)
				data[k+2-4]=0;
			else if(((16*data[k-4+2]+7*error)/16)>255)
				data[k+2-4]=255;
			else
				data[k+2-4]=(16*data[k-4+2]+7*error)/16;

			if(((16*data[k+width*4-4+2]+1*error)/16)<0)
				data[k+2+width*4-4]=0;
			else if(((16*data[k+width*4-4+2]+1*error)/16)>255)
				data[k+2+width*4-4]=255;
			else
				data[k+2+width*4-4]=(16*data[k+width*4-4+2]+1*error)/16;

			if(((16*data[k+width*4+2]+5*error)/16)<0)
				data[k+2+width*4]=0;
			else if(((16*data[k+width*4+2]+5*error)/16)>255)
				data[k+2+width*4]=255;
			else
				data[k+2+width*4]=(16*data[k+width*4+2]+5*error)/16;

			if(((16*data[k+width*4+4+2]+3*error)/16)<0)
				data[k+2+width*4+4]=0;
			else if(((16*data[k+width*4+4+2]+3*error)/16)>255)
				data[k+2+width*4+4]=255;
			else
				data[k+2+width*4+4]=(16*data[k+width*4+4+2]+3*error)/16;


			
		}
		
	}
	

   // ClearToBlack();
    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }
	
    int* alpha_f=new int[height*width];
	int* alpha_g=new int[height*width];
	for(int i =0;i<height*width*4;i+=4)
	{
		alpha_f[i/4]=data[i+3];
		alpha_g[i/4]=pImage->data[i+3];
	}

	//unsigned char* temp=new unsigned char[height*width*4];

	for(int i=0;i<height*width*4;i+=4)
	{
		data[i]=(data[i]+(1-(alpha_f[i/4]/255))*(pImage->data[i]));
		data[i+1]=(data[i+1]+(1-(alpha_f[i/4]/255))*(pImage->data[i+1]));
		data[i+2]=(data[i+2]+(1-(alpha_f[i/4]/255))*(pImage->data[i+2]));
		data[i+3]=(alpha_f[i/4]+(1-(alpha_f[i/4]/255))*alpha_g[i/4]);

		unsigned char rgb1[3];
		RGBA_To_RGB(data + i, rgb1);
		data[i]=rgb1[0];
		data[i+1]=rgb1[1];
		data[i+2]=rgb1[2];
	
	}

	delete[] alpha_f;
	delete[] alpha_g;
    //ClearToBlack();
    return true;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }
	 int* alpha_f=new int[height*width];
	int* alpha_g=new int[height*width];
	for(int i =0;i<height*width*4;i+=4)
	{
		alpha_f[i/4]=data[i+3];
		alpha_g[i/4]=pImage->data[i+3];
	}

	/*unsigned char* temp=new unsigned char[height*width*4];
	for(int i=0;i<height*width*4;i++)
	{
		temp[i]=255;
	}*/
	for(int i=0;i<height*width*4;i+=4)
	{
		data[i]=data[i]*alpha_g[i/4]/255;
		data[i+1]=data[i+1]*alpha_g[i/4]/255;
		data[i+2]=data[i+2]*alpha_g[i/4]/255;
		data[i+3]=alpha_f[i/4]*alpha_g[i/4]/255;

		unsigned char rgb1[3];
		RGBA_To_RGB(data + i, rgb1);
		data[i]=rgb1[0];
		data[i+1]=rgb1[1];
		data[i+2]=rgb1[2];

	
	}
	delete[] alpha_f;
	delete[] alpha_g;
    //ClearToBlack();
    return true;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
	// not working
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }
	double* alpha_f=new double[height*width];
	double* alpha_g=new double[height*width];
	for(int i =0;i<height*width*4;i+=4)
	{
		alpha_f[i/4]=data[i+3];
		alpha_g[i/4]=pImage->data[i+3];
	}

	/*unsigned char* temp=new unsigned char[height*width*4];
	for(int i=0;i<height*width*4;i++)
	{
		temp[i]=255;
	}*/
	for(int i=0;i<height*width*4;i+=4)
	{
		data[i]=floor(data[i]*(255-alpha_g[i/4])/255);
		data[i+1]=floor(data[i+1]*(255-alpha_g[i/4])/255);
		data[i+2]=floor(data[i+2]*(255-alpha_g[i/4])/255);
		data[i+3]=floor(alpha_f[i/4]*(255-alpha_g[i/4])/255);

		unsigned char rgb1[3];
		RGBA_To_RGB(data + i, rgb1);
		data[i]=rgb1[0];
		data[i+1]=rgb1[1];
		data[i+2]=rgb1[2];

	
	}
	delete[] alpha_f;
	delete[] alpha_g;
    //ClearToBlack();
    return true;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
	//not working
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }
	double* alpha_f=new double[height*width];
	double* alpha_g=new double[height*width];
	for(int i =0;i<height*width*4;i+=4)
	{
		alpha_f[i/4]=data[i+3];
		alpha_g[i/4]=pImage->data[i+3];
	}

	//unsigned char* temp=new unsigned char[height*width*4];

	for(int i=0;i<height*width*4;i+=4)
	{
		data[i]=floor((alpha_g[i/4]/255)*data[i]+(255-alpha_f[i/4])/255*(pImage->data[i]));
		data[i+1]=floor((alpha_g[i/4]/255)*data[i+1]+(255-alpha_f[i/4])/255*(pImage->data[i+1]));
		data[i+2]=floor((alpha_g[i/4]/255)*data[i+2]+(255-alpha_f[i/4])/255*(pImage->data[i+2]));
		data[i+3]=floor((alpha_g[i/4]/255)*alpha_f[i/4]+(255-alpha_f[i/4])/255*alpha_g[i/4]);

		unsigned char rgb1[3];
		RGBA_To_RGB(data + i, rgb1);
		data[i]=rgb1[0];
		data[i+1]=rgb1[1];
		data[i+2]=rgb1[2];
	
	}

	delete[] alpha_f;
	delete[] alpha_g;
    //ClearToBlack();
    return true;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
	//not working
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

	double* alpha_f=new double[height*width];
	double* alpha_g=new double[height*width];
	for(int i =0;i<height*width*4;i+=4)
	{
		alpha_f[i/4]=data[i+3];
		alpha_g[i/4]=pImage->data[i+3];
	}

	//unsigned char* temp=new unsigned char[height*width*4];

	for(int i=0;i<height*width*4;i+=4)
	{
		data[i]=floor((255-alpha_g[i/4])/255*data[i]+(255-alpha_f[i/4])/255*(pImage->data[i]));
		data[i+1]=floor((255-alpha_g[i/4])/255*data[i+1]+(255-alpha_f[i/4])/255*(pImage->data[i+1]));
		data[i+2]=floor((255-alpha_g[i/4])/255*data[i+2]+(255-alpha_f[i/4])/255*(pImage->data[i+2]));
		data[i+3]=floor((255-alpha_g[i/4])/255*alpha_f[i/4]+(255-alpha_f[i/4])/255*alpha_g[i/4]);

		unsigned char rgb1[3];
		RGBA_To_RGB(data + i, rgb1);
		data[i]=rgb1[0];
		data[i+1]=rgb1[1];
		data[i+2]=rgb1[2];
	
	}

	delete[] alpha_f;
	delete[] alpha_g;
    //ClearToBlack();
    return true;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
	int b_f[5][5]={{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1}};
	for(int j=0;j<height-5;j++)
	{
	for(int i=(2*width*4+8)+(j*width*4);i<(3*width*4-8)+(j*width*4);i+=4)
	{
		data[i]=
		(data[i-2*width*4-8]*b_f[0][0]+data[i-2*width*4-4]*b_f[0][1]+data[i-2*width*4]*b_f[0][2]+
		data[i-2*width*4+4]*b_f[0][3]+data[i-2*width*4+8]*b_f[0][4]+
		data[i-width*4-8]*b_f[1][0]+data[i-width*4-4]*b_f[1][1]+data[i-width*4]*b_f[1][2]+
		data[i-width*4+4]*b_f[1][3]+data[i-width*4+8]*b_f[1][4]+
		data[i-8]*b_f[2][0]+data[i-4]*b_f[2][1]+data[i]*b_f[2][2]+
		data[i+4]*b_f[2][3]+data[i+8]*b_f[2][4]+
		data[i+width*4-8]*b_f[3][0]+data[i+width*4-4]*b_f[3][1]+data[i+width*4]*b_f[3][2]+
		data[i+width*4+4]*b_f[3][3]+data[i+width*4+8]*b_f[3][4]+
		data[i+2*width*4-8]*b_f[4][0]+data[i+2*width*4-4]*b_f[4][1]+data[i+2*width*4]*b_f[4][2]+
		data[i+2*width*4+4]*b_f[4][3]+data[i+2*width*4+8]*b_f[4][4])/25;

		data[i+1]=
		(data[i+1-2*width*4-8]*b_f[0][0]+data[i+1-2*width*4-4]*b_f[0][1]+data[i+1-2*width*4]*b_f[0][2]+
		data[i+1-2*width*4+4]*b_f[0][3]+data[i+1-2*width*4+8]*b_f[0][4]+
		data[i+1-width*4-8]*b_f[1][0]+data[i+1-width*4-4]*b_f[1][1]+data[i+1-width*4]*b_f[1][2]+
		data[i+1-width*4+4]*b_f[1][3]+data[i+1-width*4+8]*b_f[1][4]+
		data[i+1-8]*b_f[2][0]+data[i+1-4]*b_f[2][1]+data[i+1]*b_f[2][2]+
		data[i+1+4]*b_f[2][3]+data[i+1+8]*b_f[2][4]+
		data[i+1+width*4-8]*b_f[3][0]+data[i+1+width*4-4]*b_f[3][1]+data[i+1+width*4]*b_f[3][2]+
		data[i+1+width*4+4]*b_f[3][3]+data[i+1+width*4+8]*b_f[3][4]+
		data[i+1+2*width*4-8]*b_f[4][0]+data[i+1+2*width*4-4]*b_f[4][1]+data[i+1+2*width*4]*b_f[4][2]+
		data[i+1+2*width*4+4]*b_f[4][3]+data[i+1+2*width*4+8]*b_f[4][4])/25;


		data[i+2]=
		(data[i+2-2*width*4-8]*b_f[0][0]+data[i+2-2*width*4-4]*b_f[0][1]+data[i+2-2*width*4]*b_f[0][2]+
		data[i+2-2*width*4+4]*b_f[0][3]+data[i+2-2*width*4+8]*b_f[0][4]+
		data[i+2-width*4-8]*b_f[1][0]+data[i+2-width*4-4]*b_f[1][1]+data[i+2-width*4]*b_f[1][2]+
		data[i+2-width*4+4]*b_f[1][3]+data[i+2-width*4+8]*b_f[1][4]+
		data[i+2-8]*b_f[2][0]+data[i+2-4]*b_f[2][1]+data[i+2]*b_f[2][2]+
		data[i+2+4]*b_f[2][3]+data[i+2+8]*b_f[2][4]+
		data[i+2+width*4-8]*b_f[3][0]+data[i+2+width*4-4]*b_f[3][1]+data[i+2+width*4]*b_f[3][2]+
		data[i+2+width*4+4]*b_f[3][3]+data[i+2+width*4+8]*b_f[3][4]+
		data[i+2+2*width*4-8]*b_f[4][0]+data[i+2+2*width*4-4]*b_f[4][1]+data[i+2+2*width*4]*b_f[4][2]+
		data[i+2+2*width*4+4]*b_f[4][3]+data[i+2+2*width*4+8]*b_f[4][4])/25;
	}
	}

    //ClearToBlack();
    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
	int b_f[5][5]={{1,2,3,2,1},{2,4,6,4,2},{3,6,9,6,3},{2,4,6,4,2},{1,2,3,2,1}};
	for(int j=0;j<height-5;j++)
	{
	for(int i=(2*width*4+8)+(j*width*4);i<(3*width*4-8)+(j*width*4);i+=4)
	{
		data[i]=
		(data[i-2*width*4-8]*b_f[0][0]+data[i-2*width*4-4]*b_f[0][1]+data[i-2*width*4]*b_f[0][2]+
		data[i-2*width*4+4]*b_f[0][3]+data[i-2*width*4+8]*b_f[0][4]+
		data[i-width*4-8]*b_f[1][0]+data[i-width*4-4]*b_f[1][1]+data[i-width*4]*b_f[1][2]+
		data[i-width*4+4]*b_f[1][3]+data[i-width*4+8]*b_f[1][4]+
		data[i-8]*b_f[2][0]+data[i-4]*b_f[2][1]+data[i]*b_f[2][2]+
		data[i+4]*b_f[2][3]+data[i+8]*b_f[2][4]+
		data[i+width*4-8]*b_f[3][0]+data[i+width*4-4]*b_f[3][1]+data[i+width*4]*b_f[3][2]+
		data[i+width*4+4]*b_f[3][3]+data[i+width*4+8]*b_f[3][4]+
		data[i+2*width*4-8]*b_f[4][0]+data[i+2*width*4-4]*b_f[4][1]+data[i+2*width*4]*b_f[4][2]+
		data[i+2*width*4+4]*b_f[4][3]+data[i+2*width*4+8]*b_f[4][4])/81;

		data[i+1]=
		(data[i+1-2*width*4-8]*b_f[0][0]+data[i+1-2*width*4-4]*b_f[0][1]+data[i+1-2*width*4]*b_f[0][2]+
		data[i+1-2*width*4+4]*b_f[0][3]+data[i+1-2*width*4+8]*b_f[0][4]+
		data[i+1-width*4-8]*b_f[1][0]+data[i+1-width*4-4]*b_f[1][1]+data[i+1-width*4]*b_f[1][2]+
		data[i+1-width*4+4]*b_f[1][3]+data[i+1-width*4+8]*b_f[1][4]+
		data[i+1-8]*b_f[2][0]+data[i+1-4]*b_f[2][1]+data[i+1]*b_f[2][2]+
		data[i+1+4]*b_f[2][3]+data[i+1+8]*b_f[2][4]+
		data[i+1+width*4-8]*b_f[3][0]+data[i+1+width*4-4]*b_f[3][1]+data[i+1+width*4]*b_f[3][2]+
		data[i+1+width*4+4]*b_f[3][3]+data[i+1+width*4+8]*b_f[3][4]+
		data[i+1+2*width*4-8]*b_f[4][0]+data[i+1+2*width*4-4]*b_f[4][1]+data[i+1+2*width*4]*b_f[4][2]+
		data[i+1+2*width*4+4]*b_f[4][3]+data[i+1+2*width*4+8]*b_f[4][4])/81;


		data[i+2]=
		(data[i+2-2*width*4-8]*b_f[0][0]+data[i+2-2*width*4-4]*b_f[0][1]+data[i+2-2*width*4]*b_f[0][2]+
		data[i+2-2*width*4+4]*b_f[0][3]+data[i+2-2*width*4+8]*b_f[0][4]+
		data[i+2-width*4-8]*b_f[1][0]+data[i+2-width*4-4]*b_f[1][1]+data[i+2-width*4]*b_f[1][2]+
		data[i+2-width*4+4]*b_f[1][3]+data[i+2-width*4+8]*b_f[1][4]+
		data[i+2-8]*b_f[2][0]+data[i+2-4]*b_f[2][1]+data[i+2]*b_f[2][2]+
		data[i+2+4]*b_f[2][3]+data[i+2+8]*b_f[2][4]+
		data[i+2+width*4-8]*b_f[3][0]+data[i+2+width*4-4]*b_f[3][1]+data[i+2+width*4]*b_f[3][2]+
		data[i+2+width*4+4]*b_f[3][3]+data[i+2+width*4+8]*b_f[3][4]+
		data[i+2+2*width*4-8]*b_f[4][0]+data[i+2+2*width*4-4]*b_f[4][1]+data[i+2+2*width*4]*b_f[4][2]+
		data[i+2+2*width*4+4]*b_f[4][3]+data[i+2+2*width*4+8]*b_f[4][4])/81;
	}
	}

    //ClearToBlack();
    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{

	int b_f[5][5]={{1,4,6,4,1},{4,16,24,16,4},{6,24,36,24,6},{4,16,24,16,4},{1,4,6,4,1}};
	for(int j=0;j<height-5;j++)
	{
	for(int i=(2*width*4+8)+(j*width*4);i<(3*width*4-8)+(j*width*4);i+=4)
	{
		data[i]=
		(data[i-2*width*4-8]*b_f[0][0]+data[i-2*width*4-4]*b_f[0][1]+data[i-2*width*4]*b_f[0][2]+
		data[i-2*width*4+4]*b_f[0][3]+data[i-2*width*4+8]*b_f[0][4]+
		data[i-width*4-8]*b_f[1][0]+data[i-width*4-4]*b_f[1][1]+data[i-width*4]*b_f[1][2]+
		data[i-width*4+4]*b_f[1][3]+data[i-width*4+8]*b_f[1][4]+
		data[i-8]*b_f[2][0]+data[i-4]*b_f[2][1]+data[i]*b_f[2][2]+
		data[i+4]*b_f[2][3]+data[i+8]*b_f[2][4]+
		data[i+width*4-8]*b_f[3][0]+data[i+width*4-4]*b_f[3][1]+data[i+width*4]*b_f[3][2]+
		data[i+width*4+4]*b_f[3][3]+data[i+width*4+8]*b_f[3][4]+
		data[i+2*width*4-8]*b_f[4][0]+data[i+2*width*4-4]*b_f[4][1]+data[i+2*width*4]*b_f[4][2]+
		data[i+2*width*4+4]*b_f[4][3]+data[i+2*width*4+8]*b_f[4][4])/256;

		data[i+1]=
		(data[i+1-2*width*4-8]*b_f[0][0]+data[i+1-2*width*4-4]*b_f[0][1]+data[i+1-2*width*4]*b_f[0][2]+
		data[i+1-2*width*4+4]*b_f[0][3]+data[i+1-2*width*4+8]*b_f[0][4]+
		data[i+1-width*4-8]*b_f[1][0]+data[i+1-width*4-4]*b_f[1][1]+data[i+1-width*4]*b_f[1][2]+
		data[i+1-width*4+4]*b_f[1][3]+data[i+1-width*4+8]*b_f[1][4]+
		data[i+1-8]*b_f[2][0]+data[i+1-4]*b_f[2][1]+data[i+1]*b_f[2][2]+
		data[i+1+4]*b_f[2][3]+data[i+1+8]*b_f[2][4]+
		data[i+1+width*4-8]*b_f[3][0]+data[i+1+width*4-4]*b_f[3][1]+data[i+1+width*4]*b_f[3][2]+
		data[i+1+width*4+4]*b_f[3][3]+data[i+1+width*4+8]*b_f[3][4]+
		data[i+1+2*width*4-8]*b_f[4][0]+data[i+1+2*width*4-4]*b_f[4][1]+data[i+1+2*width*4]*b_f[4][2]+
		data[i+1+2*width*4+4]*b_f[4][3]+data[i+1+2*width*4+8]*b_f[4][4])/256;


		data[i+2]=
		(data[i+2-2*width*4-8]*b_f[0][0]+data[i+2-2*width*4-4]*b_f[0][1]+data[i+2-2*width*4]*b_f[0][2]+
		data[i+2-2*width*4+4]*b_f[0][3]+data[i+2-2*width*4+8]*b_f[0][4]+
		data[i+2-width*4-8]*b_f[1][0]+data[i+2-width*4-4]*b_f[1][1]+data[i+2-width*4]*b_f[1][2]+
		data[i+2-width*4+4]*b_f[1][3]+data[i+2-width*4+8]*b_f[1][4]+
		data[i+2-8]*b_f[2][0]+data[i+2-4]*b_f[2][1]+data[i+2]*b_f[2][2]+
		data[i+2+4]*b_f[2][3]+data[i+2+8]*b_f[2][4]+
		data[i+2+width*4-8]*b_f[3][0]+data[i+2+width*4-4]*b_f[3][1]+data[i+2+width*4]*b_f[3][2]+
		data[i+2+width*4+4]*b_f[3][3]+data[i+2+width*4+8]*b_f[3][4]+
		data[i+2+2*width*4-8]*b_f[4][0]+data[i+2+2*width*4-4]*b_f[4][1]+data[i+2+2*width*4]*b_f[4][2]+
		data[i+2+2*width*4+4]*b_f[4][3]+data[i+2+2*width*4+8]*b_f[4][4])/256;
	}
	}
    //ClearToBlack();
    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    ClearToBlack();
   return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
	int x=0;
	int y=0;
	int z=0;
	unsigned char* temp_data=new unsigned char[height*width*4];


	int b_f[5][5]={{-1,-4,-6,-4,-1},{-4,-16,-24,-16,-4},{-6,-24,220,-24,-6},{-4,-16,-24,-16,-4},{-1,-4,-6,-4,-1}};
	for(int j=0;j<height-5;j++)
	{
	for(int i=(2*width*4+8)+(j*width*4);i<(3*width*4-8)+(j*width*4);i+=4)
	{
		x=
		floor((double)((data[i-2*width*4-8]*b_f[0][0]+data[i-2*width*4-4]*b_f[0][1]+data[i-2*width*4]*b_f[0][2]+
		data[i-2*width*4+4]*b_f[0][3]+data[i-2*width*4+8]*b_f[0][4]+
		data[i-width*4-8]*b_f[1][0]+data[i-width*4-4]*b_f[1][1]+data[i-width*4]*b_f[1][2]+
		data[i-width*4+4]*b_f[1][3]+data[i-width*4+8]*b_f[1][4]+
		data[i-8]*b_f[2][0]+data[i-4]*b_f[2][1]+data[i]*b_f[2][2]+
		data[i+4]*b_f[2][3]+data[i+8]*b_f[2][4]+
		data[i+width*4-8]*b_f[3][0]+data[i+width*4-4]*b_f[3][1]+data[i+width*4]*b_f[3][2]+
		data[i+width*4+4]*b_f[3][3]+data[i+width*4+8]*b_f[3][4]+
		data[i+2*width*4-8]*b_f[4][0]+data[i+2*width*4-4]*b_f[4][1]+data[i+2*width*4]*b_f[4][2]+
		data[i+2*width*4+4]*b_f[4][3]+data[i+2*width*4+8]*b_f[4][4])/256));

		y=
		floor((double)((data[i+1-2*width*4-8]*b_f[0][0]+data[i+1-2*width*4-4]*b_f[0][1]+data[i+1-2*width*4]*b_f[0][2]+
		data[i+1-2*width*4+4]*b_f[0][3]+data[i+1-2*width*4+8]*b_f[0][4]+
		data[i+1-width*4-8]*b_f[1][0]+data[i+1-width*4-4]*b_f[1][1]+data[i+1-width*4]*b_f[1][2]+
		data[i+1-width*4+4]*b_f[1][3]+data[i+1-width*4+8]*b_f[1][4]+
		data[i+1-8]*b_f[2][0]+data[i+1-4]*b_f[2][1]+data[i+1]*b_f[2][2]+
		data[i+1+4]*b_f[2][3]+data[i+1+8]*b_f[2][4]+
		data[i+1+width*4-8]*b_f[3][0]+data[i+1+width*4-4]*b_f[3][1]+data[i+1+width*4]*b_f[3][2]+
		data[i+1+width*4+4]*b_f[3][3]+data[i+1+width*4+8]*b_f[3][4]+
		data[i+1+2*width*4-8]*b_f[4][0]+data[i+1+2*width*4-4]*b_f[4][1]+data[i+1+2*width*4]*b_f[4][2]+
		data[i+1+2*width*4+4]*b_f[4][3]+data[i+1+2*width*4+8]*b_f[4][4])/256));


		z=
		floor((double)((data[i+2-2*width*4-8]*b_f[0][0]+data[i+2-2*width*4-4]*b_f[0][1]+data[i+2-2*width*4]*b_f[0][2]+
		data[i+2-2*width*4+4]*b_f[0][3]+data[i+2-2*width*4+8]*b_f[0][4]+
		data[i+2-width*4-8]*b_f[1][0]+data[i+2-width*4-4]*b_f[1][1]+data[i+2-width*4]*b_f[1][2]+
		data[i+2-width*4+4]*b_f[1][3]+data[i+2-width*4+8]*b_f[1][4]+
		data[i+2-8]*b_f[2][0]+data[i+2-4]*b_f[2][1]+data[i+2]*b_f[2][2]+
		data[i+2+4]*b_f[2][3]+data[i+2+8]*b_f[2][4]+
		data[i+2+width*4-8]*b_f[3][0]+data[i+2+width*4-4]*b_f[3][1]+data[i+2+width*4]*b_f[3][2]+
		data[i+2+width*4+4]*b_f[3][3]+data[i+2+width*4+8]*b_f[3][4]+
		data[i+2+2*width*4-8]*b_f[4][0]+data[i+2+2*width*4-4]*b_f[4][1]+data[i+2+2*width*4]*b_f[4][2]+
		data[i+2+2*width*4+4]*b_f[4][3]+data[i+2+2*width*4+8]*b_f[4][4])/256));


		if(x<=0)
			temp_data[i]=0;
		else
			temp_data[i]=x;

		if(y<=0)
			temp_data[i+1]=0;
		else
			temp_data[i+1]=y;

		if(z<=0)
			temp_data[i+2]=0;
		else
			temp_data[i+2]=z;

		temp_data[i+3]=255;
	
	}
	}

	for(int i=0;i<height*width*4;i++)
	{
		data[i]=temp_data[i];
	}

	delete[] temp_data;
	//cout<<x<<" "<<y<<" "<<z<<endl;
    //ClearToBlack();
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
	int x=0;
	int y=0;
	int z=0;
	unsigned char* temp_data=new unsigned char[height*width*4];
	unsigned char* temp_data2=new unsigned char[height*width*4];

	int b_f[5][5]={{-1,-4,-6,-4,-1},{-4,-16,-24,-16,-4},{-6,-24,220,-24,-6},{-4,-16,-24,-16,-4},{-1,-4,-6,-4,-1}};
	for(int j=0;j<height-5;j++)
	{
	for(int i=(2*width*4+8)+(j*width*4);i<(3*width*4-8)+(j*width*4);i+=4)
	{
		x=
		floor((double)((data[i-2*width*4-8]*b_f[0][0]+data[i-2*width*4-4]*b_f[0][1]+data[i-2*width*4]*b_f[0][2]+
		data[i-2*width*4+4]*b_f[0][3]+data[i-2*width*4+8]*b_f[0][4]+
		data[i-width*4-8]*b_f[1][0]+data[i-width*4-4]*b_f[1][1]+data[i-width*4]*b_f[1][2]+
		data[i-width*4+4]*b_f[1][3]+data[i-width*4+8]*b_f[1][4]+
		data[i-8]*b_f[2][0]+data[i-4]*b_f[2][1]+data[i]*b_f[2][2]+
		data[i+4]*b_f[2][3]+data[i+8]*b_f[2][4]+
		data[i+width*4-8]*b_f[3][0]+data[i+width*4-4]*b_f[3][1]+data[i+width*4]*b_f[3][2]+
		data[i+width*4+4]*b_f[3][3]+data[i+width*4+8]*b_f[3][4]+
		data[i+2*width*4-8]*b_f[4][0]+data[i+2*width*4-4]*b_f[4][1]+data[i+2*width*4]*b_f[4][2]+
		data[i+2*width*4+4]*b_f[4][3]+data[i+2*width*4+8]*b_f[4][4])/256));

		y=
		floor((double)((data[i+1-2*width*4-8]*b_f[0][0]+data[i+1-2*width*4-4]*b_f[0][1]+data[i+1-2*width*4]*b_f[0][2]+
		data[i+1-2*width*4+4]*b_f[0][3]+data[i+1-2*width*4+8]*b_f[0][4]+
		data[i+1-width*4-8]*b_f[1][0]+data[i+1-width*4-4]*b_f[1][1]+data[i+1-width*4]*b_f[1][2]+
		data[i+1-width*4+4]*b_f[1][3]+data[i+1-width*4+8]*b_f[1][4]+
		data[i+1-8]*b_f[2][0]+data[i+1-4]*b_f[2][1]+data[i+1]*b_f[2][2]+
		data[i+1+4]*b_f[2][3]+data[i+1+8]*b_f[2][4]+
		data[i+1+width*4-8]*b_f[3][0]+data[i+1+width*4-4]*b_f[3][1]+data[i+1+width*4]*b_f[3][2]+
		data[i+1+width*4+4]*b_f[3][3]+data[i+1+width*4+8]*b_f[3][4]+
		data[i+1+2*width*4-8]*b_f[4][0]+data[i+1+2*width*4-4]*b_f[4][1]+data[i+1+2*width*4]*b_f[4][2]+
		data[i+1+2*width*4+4]*b_f[4][3]+data[i+1+2*width*4+8]*b_f[4][4])/256));


		z=
		floor((double)((data[i+2-2*width*4-8]*b_f[0][0]+data[i+2-2*width*4-4]*b_f[0][1]+data[i+2-2*width*4]*b_f[0][2]+
		data[i+2-2*width*4+4]*b_f[0][3]+data[i+2-2*width*4+8]*b_f[0][4]+
		data[i+2-width*4-8]*b_f[1][0]+data[i+2-width*4-4]*b_f[1][1]+data[i+2-width*4]*b_f[1][2]+
		data[i+2-width*4+4]*b_f[1][3]+data[i+2-width*4+8]*b_f[1][4]+
		data[i+2-8]*b_f[2][0]+data[i+2-4]*b_f[2][1]+data[i+2]*b_f[2][2]+
		data[i+2+4]*b_f[2][3]+data[i+2+8]*b_f[2][4]+
		data[i+2+width*4-8]*b_f[3][0]+data[i+2+width*4-4]*b_f[3][1]+data[i+2+width*4]*b_f[3][2]+
		data[i+2+width*4+4]*b_f[3][3]+data[i+2+width*4+8]*b_f[3][4]+
		data[i+2+2*width*4-8]*b_f[4][0]+data[i+2+2*width*4-4]*b_f[4][1]+data[i+2+2*width*4]*b_f[4][2]+
		data[i+2+2*width*4+4]*b_f[4][3]+data[i+2+2*width*4+8]*b_f[4][4])/256));


		if(x<=0)
			temp_data[i]=0;
		else
			temp_data[i]=x;

		if(y<=0)
			temp_data[i+1]=0;
		else
			temp_data[i+1]=y;

		if(z<=0)
			temp_data[i+2]=0;
		else
			temp_data[i+2]=z;

		temp_data[i+3]=255;
	
	}
	}

	for(int i=0;i<height*width*4;i++)
	{
		int sum=data[i]+temp_data[i];
		if(sum>255)
			temp_data2[i]=255;
		else
			temp_data2[i]=sum;
	}
	for(int i=0;i<height*width*4;i++)
	{
		data[i]=temp_data2[i];
	}

	delete[] temp_data;
	delete[] temp_data2;
    //ClearToBlack();
    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

