//////////////////////////////////////////////////////////////////////////////
//                UWHDR: High Dynamic Range Imaging
//                *********************************
//
// Software by:  Chaman Singh Verma
//               Department of Computer Science  
//               University of Wisconsin, Madison.
//               USA
// Date:         29th September, 2009
//
// Information:  This software implements the radiance map of the Debevec's
//               method as described in the Siggraph 97 Paper in C++.
//               
//               Recovering High Dynamic Range Radiance Maps from Photographs
//               Paul E. Debevec and Jitendra Malik.
// External Software Dependencies:
// (1)  ImageMagick++ 
// (2)  RGBE 
// (3)  CLAPACK
//
//////////////////////////////////////////////////////////////////////////////`
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <limits>
#include <sys/time.h>

#include <fstream>

#include <iostream>
#include <Magick++.h>

#include <vector>
#include <list>

#include "rgbe.h"

extern "C" {
#include "blaswrap.h"
#include "f2c.h"
#include "clapack.h"
}
    
using namespace std;
using namespace Magick;

typedef unsigned short ushort;

class StopWatch
{
  public:
   void start() { gettimeofday( &start_time, 0); }
   void stop()  { gettimeofday( &end_time,   0); }
   int  elapsed()
   {
      return (end_time.tv_sec  - start_time.tv_sec) + 
             (end_time.tv_usec - start_time.tv_usec)*1.0E-06;
   }
  private:
   struct timeval  start_time, end_time;
};
 
///////////////////////////////////////////////////////////////////////////////

ushort eight_bit_pixel( const Quantum &q) 
{
    double dval = 255.0*((double) q/65535.0);
    if( dval < 0.0) dval = 0.0;
    if( dval > 255) dval = 255;
    return round(dval);
}

///////////////////////////////////////////////////////////////////////////////

int load_images_and_exposure_data( const string directory, 
                                   vector<Image> &images, vector<double> &exposeTime)
{
    string fname;
   
    images.clear();
    exposeTime.clear();
   
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(directory.c_str())) == NULL) {
        cout << "Error Cann't open directory " << directory << endl;
        return 1;
    }

    vector<string> files;
    while ((dirp = readdir(dp)) != NULL) {
        fname = dirp->d_name;
        if( fname.find(".jpg") != string::npos) 
            files.push_back(fname);
    }
    closedir(dp);

    sort( files.begin(), files.end() );

    double ev;
    Magick::Image image;

    string exif, strval;
    size_t pos;
    for( int i = 0; i < files.size(); i++) {
            image.read(directory + "/" + files[i] );
            exif = image.attribute("exif:ExposureTime");
            pos  = exif.find("/");
            if( pos != string::npos) {
                strval = exif.substr(0,pos);
                int numer = atoi( strval.c_str() );
                strval = exif.substr(pos+1);
                int denom = atoi( strval.c_str() );
                ev   = ( double)numer/(double)denom;
                int nrows = image.rows();
                int ncols = image.columns();
                image.resize( Geometry(0.5*nrows, 0.5*ncols) );            
                images.push_back(image);
                exposeTime.push_back(ev);
                cout << files[i] << " Exposure time " << ev << endl;
            } else 
                cout << "Warning: Exposure is not rational number " << files[i] << "  " << exif << endl;
    }
    cout << "Info: number of images loaded " << images.size() << endl;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void set_pixel_weights( vector<double> &pixelWeight, int zMin = 0, int zMax = 255)
{
   int nSize = zMax - zMin + 1;
   pixelWeight.resize( nSize );

   int zMid = 0.5*( zMin + zMax );

   for( int i = 0; i < nSize; i++) {
        int z = zMin + i;
        if( z <= zMid ) pixelWeight[i] =  z - zMin;
        if( z  > zMid ) pixelWeight[i] =  zMax - z;
   }
}

///////////////////////////////////////////////////////////////////////////////

int random_index( int minValue, int maxValue)
{
   return rand()%maxValue;
}

///////////////////////////////////////////////////////////////////////////////
void select_random_pos( const vector<ushort> &pixelValues, vector<int> &pixelPos)
{
   size_t numPixels = pixelValues.size();

   vector<int> pixelIndex(numPixels);
   int numSamePixels;

   pixelPos.clear();

   size_t numTotal;
   for( int ipixel = 5; ipixel < 250; ipixel += 1) {
        numSamePixels = 0;
        for( int j = 0; j < numPixels; j++) {
             if( pixelValues[j] == ipixel )
                 pixelIndex[numSamePixels++] = j;
        }

        if( numSamePixels ) {
            int pos = random_index( 0, numSamePixels);
            assert( pos >= 0 && pos < numSamePixels );
            pixelPos.push_back( pixelIndex[pos] );
        }
   }

   sort( pixelPos.begin(), pixelPos.end() );
   vector<int>::iterator ibegin, iter, iter_end;
   iter_end = unique( pixelPos.begin(), pixelPos.end() );

   vector<int> unique_pos;
   unique_pos.reserve( pixelPos.size() );
   for( iter = pixelPos.begin(); iter != iter_end; ++iter)
        unique_pos.push_back( *iter );

   pixelPos = unique_pos;
}

///////////////////////////////////////////////////////////////////////////////

int getIndex( int i, int j, int numRows, int numCols)
{
   return j*numRows + i;
}

///////////////////////////////////////////////////////////////////////////////

void generate_matrix( const vector<vector<ushort> > &imageFrames, 
                      const vector<int> &pixelPos, 
                      const vector<double> &exposeTime, double lambda, 
                      vector<double> &A, int &numRows, int &numCols, 
                      vector<double> &b)
{
    vector<double> pixelWeight(256);
    set_pixel_weights( pixelWeight );

    int numFrames = imageFrames.size();
    int numSamplePixels = pixelPos.size();

    numRows = numSamplePixels*numFrames + 255;
    numCols = 256 + numSamplePixels;

    assert( numRows > numCols );

    A.resize(numRows*numCols);
    for( int i = 0; i < numRows*numCols; i++) 
         A[i] = 0.0;

    b.resize(numRows);
    for( int i = 0; i < numRows; i++) 
         b[i] = 0.0;

    int offset, rowindex = 0;
    for( int ipos = 0; ipos < numSamplePixels; ipos++) {
         int currpos = pixelPos[ipos];
         for( int iframe = 0; iframe < numFrames; iframe++) {
              int zij =  imageFrames[iframe][currpos];
              int wij =  pixelWeight[zij];

              offset = getIndex(rowindex, zij, numRows, numCols);
              A[offset] = (double)wij;

              offset = getIndex(rowindex, 256 + ipos, numRows, numCols);
              A[offset] = -(double)wij;

              b[rowindex] = wij*log( exposeTime[iframe] );
              rowindex++;
         }
    }
    
    rowindex = numFrames*numSamplePixels;
    for( int ipixel = 1; ipixel < 255; ipixel++) {
         offset      = getIndex( rowindex, ipixel-1, numRows, numCols);
         A[offset]   = lambda*pixelWeight[ipixel];

         offset      = getIndex( rowindex, ipixel, numRows, numCols);
         A[offset]   = -2.0*lambda*pixelWeight[ipixel];

         offset      = getIndex( rowindex, ipixel+1, numRows, numCols);
         A[offset]   = lambda*pixelWeight[ipixel];

         rowindex++;
    }

    offset      = getIndex( numRows-1, 127, numRows, numCols);
    A[offset]   = 1.0;
}

////////////////////////////////////////////////////////////////////////////////////

void least_square_solution(double *A, double nrows, double ncols, double *b, double *x)
{
   integer  m = nrows;
   integer  n = ncols;
   integer  nrhs = 1;

   integer lda = m;
   integer ldb = max(m,n);

   vector<double> singular( max(m,n) );

   double rcond = 0.01;  // Remaining eigenvalues will be treated zero.
   integer rank, info;

   integer lwork = 3*min(m,n) + max( 2*min(m,n), max(m,n));
   vector<double> work(lwork);

   dgelss_(&m, &n, &nrhs, A, &lda, b, &ldb, &singular[0],
           &rcond, &rank, &work[0], &lwork, &info);

   for( int i = 0; i < m; i++) 
         x[i] = b[i];
}

///////////////////////////////////////////////////////////////////////////////
void radiance_map( const vector< vector<ushort> > &pixelFrames, 
                   const vector<double> &exposeTime, 
                   const vector<int> &pixelSamplePos, double lambda, 
                   vector<double> &lnE, vector<double> &g, 
                   vector<float> &hdr_color)
{
    vector<double> A, b, x;
    int numRows, numCols, nSize;

    generate_matrix( pixelFrames, pixelSamplePos, exposeTime, lambda, 
                     A, numRows, numCols, b);

    x.resize( numRows );
    least_square_solution(&A[0], numRows, numCols, &b[0], &x[0]);

    nSize = 256;
    g.resize(nSize);
    for( int i = 0; i < nSize; i++) g[i] = x[i];

    // Make corrections at the end points;
    double gmin = g[0];
    for( int i = 0; i < 50; i++) gmin = gmin < g[i] ? gmin : g[i];
    for( int i = 0; i < 50; i++) {
         if( g[i] == gmin) break;
         g[i] = gmin;
    }

    double gmax = g[255];
    for( int i = 254; i > 200; i--) gmax = gmax > g[i] ? gmax : g[i];
    for( int i = 255; i > 200; i--) {
         if( g[i] == gmax) break;
         g[i] = gmax;
    }

    vector<double> pixelWeight(256);
    set_pixel_weights( pixelWeight );

    int numFrames = pixelFrames.size();
    vector<double> E(256);
    lnE.resize(256);
    for( int ip = 1; ip < 255; ip++) {
         double sum_numer = 0.0;
         double sum_denom = 0.0;
         int  zij = ip;
         for( int j = 0; j < numFrames; j++) {
              double Ei = g[zij] - log( exposeTime[j] );
              sum_numer += pixelWeight[zij]*Ei;
              sum_denom += pixelWeight[zij];
         }
         assert( sum_denom > 0.0);
         lnE[ip] = sum_numer/sum_denom;
         E[ip] = exp(lnE[ip]);
    }

    E[0]   = E[1];
    lnE[0] = lnE[1];

    E[255]   = E[254];
    lnE[255] = lnE[254];
    
    int numPixels = pixelFrames[0].size(); // Same for all images.
    int midframe  = pixelFrames.size()/2;

    hdr_color.resize( numPixels );
    for( int ip = 0; ip < numPixels; ip++) {
         int zij  = pixelFrames[midframe][ip];
         hdr_color[ip] = E[zij];
    }
}

///////////////////////////////////////////////////////////////////////////////

void  extract_color_from_images( const vector<Image> &images, int whichColor, 
                                 vector< vector<ushort> > &pixelFrames)
{
   ushort c;
   int numRows, numColums, numFrames, index = 0;

   numFrames = images.size();
   pixelFrames.resize( numFrames );
   
   for( int k = 0; k < images.size(); k++) {
        int numRows = images[k].rows();
        int numCols = images[k].columns();
        assert( images[k].colorSpace() == RGBColorspace );
        assert( images[k].classType()  == DirectClass   );
        pixelFrames[k].resize( numRows*numCols );
        
        index = 0;
        for( int i = 0; i < numRows; i++) {
        for( int j = 0; j < numCols; j++) {
                Color pixelColor = images[k].pixelColor(j,i);
                switch( whichColor )
                {
                    case 0:
                         c = eight_bit_pixel( pixelColor.redQuantum() );
                         break;
                    case 1:
                         c = eight_bit_pixel( pixelColor.greenQuantum());
                         break;
                    case 2:
                         c = eight_bit_pixel( pixelColor.blueQuantum() );
                         break;
                 }
                 pixelFrames[k][index++] = c;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void write_response_curve(const vector<double> &lnE, const vector<double> &g, int whichColor)
{
   ostringstream oss1;
   oss1 << "lnE" << whichColor << ".dat";
   ofstream ofile1( oss1.str().c_str(), ios::out);
   for( int ip = 0; ip < 256; ip++) 
        ofile1 << ip << " " << lnE[ip] << endl;

/*
   ostringstream oss2;
   oss2 << "g" << whichColor << ".dat";
   ofstream ofile2( oss2.str().c_str(), ios::out);
   for( int ip = 0; ip < 256; ip++) 
        ofile2 << ip << " " << g[ip]  << endl;
*/

}

///////////////////////////////////////////////////////////////////////////////

void makehdr( const vector<Image> &images, const vector<double> &exposeTime, 
              double lambda, const string filename)
{
   int midplane;
   vector< vector<ushort> > pixelFrames;
   vector<int> pixelSamplePos;

   // Red Channel 
   cout << " Processing Red Channel ... " << endl;
   vector<float>  red_map;
   vector<double> lnRE, rg;
   extract_color_from_images( images, 0, pixelFrames);
   midplane = pixelFrames.size()/2;
   select_random_pos( pixelFrames[midplane], pixelSamplePos);
   cout << "# Sample Pixels : " << pixelSamplePos.size() << endl;
   radiance_map( pixelFrames, exposeTime, pixelSamplePos, lambda, lnRE, rg, red_map);
   write_response_curve( lnRE, rg, 0);

   // Green Channel 
   cout << " Processing Green Channel ... " << endl;
   vector<float>  green_map;
   vector<double> lnGE, gg;
   extract_color_from_images( images, 1, pixelFrames);
   midplane = pixelFrames.size()/2;
   select_random_pos( pixelFrames[midplane], pixelSamplePos);
   cout << "# Sample Pixels : " << pixelSamplePos.size() << endl;
   radiance_map( pixelFrames, exposeTime, pixelSamplePos, lambda, lnGE, gg, green_map);
   write_response_curve( lnGE, gg, 1);

   // Blue Channel 
   cout << " Processing Blue Channel ... " << endl;
   vector<float>  blue_map;
   vector<double> lnBE, bg;
   extract_color_from_images( images, 2, pixelFrames);
   midplane = pixelFrames.size()/2;
   select_random_pos( pixelFrames[midplane], pixelSamplePos);
   cout << "# Sample Pixels : " << pixelSamplePos.size() << endl;
   radiance_map( pixelFrames, exposeTime, pixelSamplePos, lambda, lnBE, bg, blue_map);
   write_response_curve( lnBE, bg, 2 );

   int numRows = images[0].rows();
   int numCols = images[0].columns();
   int numPixels = numRows*numCols;

   vector<float> hdr_data(3*numPixels );
   for( int i = 0; i < numPixels; i++) {
        hdr_data[3*i+0] = red_map[i];
        hdr_data[3*i+1] = green_map[i];
        hdr_data[3*i+2] = blue_map[i];
   }

   FILE *fpout = fopen( "scene.hdr", "w");
   assert(fpout != NULL );
   RGBE_WriteHeader(fpout, numCols, numRows, NULL );
   RGBE_WritePixels(fpout, &hdr_data[0], numPixels);

}

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) 
{
   if( argc != 3 ) {
       cout << "Usage: executable directory lambda " << endl;
       exit(0);
   }
   
   vector<Image> images;
   vector<double> exposeTime;

   string directory = string( argv[1] );
   double lambda    = atof( argv[2] );

   load_images_and_exposure_data( directory, images, exposeTime);

   string output = "scene.hdr";
   makehdr( images, exposeTime, lambda, output);
}

///////////////////////////////////////////////////////////////////////////////
