#include <iostream>
#include <cassert>
#include <cstring>
#include <vector>
#include <ctime>
#include <fstream>
#include <math.h>

#include <Eigen/Dense>
#include "ImageRGB.hpp"
#include "ioJPG.hpp"
#include "exif.h"

using namespace kn;




int getWeight(int z, int zMin = 0 , int zMax = 255){
  if( z <= (zMin + zMax)/2.0){
    return z - zMin;
  }
  else{
    return zMax - z;
  }
}



//////////////////////////////////////////////////////////////////////////////////////////////
// open the file "filename" and copy the file content into a string (required for exif reader)
std::string fileToString(const std::string& filename){
  std::ifstream file(filename.c_str());//, std::ios::binary);
  if (!file) return "";
  std::string str(std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>()));
  return str;
}


//////////////////////////////////////////////////////////////////////////////////////////////
void drawPixels(ImageRGB8u &image, const std::vector<Eigen::Vector2i> &pixels){
  for(uint i=0; i<pixels.size(); ++i){
    image(pixels[i][0],pixels[i][1])[0] = 255;
    image(pixels[i][0],pixels[i][1])[1] = 0;
    image(pixels[i][0],pixels[i][1])[2] = 0;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////
void exifParsingError(const int parseSuccess){
  switch(parseSuccess){
    case 1982 : // No JPEG markers found in buffer, possibly invalid JPEG file
      std::cout << "exif parsing error : PARSE_EXIF_ERROR_NO_JPEG" << std::endl;
    break;
    case 1983 : // No EXIF header found in JPEG file.
      std::cout << "exif parsing error : PARSE_EXIF_ERROR_NO_EXIF" << std::endl;
    break;
    case 1984 : // Byte alignment specified in EXIF file was unknown (not Motorola or Intel).
      std::cout << "exif parsing error : PARSE_EXIF_ERROR_UNKNOWN_BYTEALIGN" << std::endl;
    break;
    case 1985 : // EXIF header was found, but data was corrupted.
      std::cout << "exif parsing error : PARSE_EXIF_ERROR_CORRUPT" << std::endl;
    break;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////
void loadImages(const int argc, char **argv, std::vector<ImageRGB8u> &images, std::vector<double> &exposure){

  for(int i=1; i<argc; ++i){
    //load the image
    ImageRGB8u imgageTmp;
    std::cout << "loading '" << argv[i] << "' ...";
    loadJPG(imgageTmp,argv[i]);
    images.push_back(imgageTmp);
    std::cout << " done" << std::endl;

    // load the exposure duration from the exif
    EXIFInfo exifReader;
    int parseSuccess = exifReader.parseFrom(fileToString(argv[i]));
    if(parseSuccess != PARSE_EXIF_SUCCESS){
      exifParsingError(parseSuccess);
      exit(0);
    }
    std::cout << " wxh : " << exifReader.ImageWidth << " x " << exifReader.ImageHeight << std::endl;
    std::cout << " exposure : " << exifReader.ExposureTime << " s" << std::endl;
    std::cout << " flash : " << ((exifReader.Flash==0)?"no":"yes") << std::endl;
    std::cout << " camera : " << exifReader.Model << std::endl;
    std::cout << " ISO : " << exifReader.ISOSpeedRatings << std::endl;
    std::cout << " apperture : " << exifReader.FNumber << std::endl;
    std::cout << std::endl;

    // update exposure
    exposure.push_back((double)exifReader.ExposureTime);
    //exposure.push_back((i+1) * 5.0);
  }
}


Eigen::VectorXd responseRecovery (const std::vector<Eigen::MatrixXi> &images,
const std::vector<double> &exposure,
const std::vector<Eigen::Vector2i> &pixels,
const std::vector<double> &zWeight,
const int valueMin,
const int valueMax,
const double lambda){

  int image_nb = images.size();
  int sample_nb = pixels.size();
  std::cout<<"image_nb : "<<image_nb<<std::endl;
  std::cout<<"sample_nb : "<<sample_nb<<std::endl;
  std::cout<<"Creation vecteur b"<<std::endl;

  /*********************************/
  /* Création Vecteur b */
  /*********************************/
  Eigen::VectorXd b = Eigen::VectorXd::Zero(sample_nb * image_nb + 254 + 1); //+1??

  //partie supérieure de b
  for(int current_image = 0; current_image < image_nb; ++current_image){
    for(int current_px = 0; current_px < sample_nb; ++current_px){
      std::cout << "(" << current_image << "e image "<< current_px<<"e pixel) : (" << pixels[current_px](0) << ", " << pixels[current_px](1) << ")"<< std::endl;
      int currentZ = images[current_image](pixels[current_px](0), pixels[current_px](1));
      b(current_image * sample_nb + current_px) = zWeight[currentZ] * log(exposure[current_image]); //* calcul_poids(pixel(current_image,current_px))
    }

  }
  std::cout << "OK bb" << std::endl;
  b(sample_nb * image_nb + 254) = 1;

  std::cout<<"C'EST GENIAL"<<std::endl;
  std::cout<<"matrice A"<<std::endl;
  /*********************************/
  /* Création Matrice A */
  /*********************************/
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(254.0 + sample_nb * image_nb + 1, 256.0 + sample_nb);

  std::cout<<"part sup gauche"<<std::endl;
  //partie supérieure gauche
  for (int current_image = 0; current_image < image_nb; ++current_image){
    for (int current_px = 0; current_px < sample_nb; ++current_px){
      int z = images[current_image](pixels[current_px](0),pixels[current_px](1));

      A(current_px + current_image*sample_nb, z) = zWeight[z];
    }
  }

  std::cout<<"part sup droite"<<std::endl;
  //partie supérieure droite
  for (int current_image = 0; current_image < image_nb; ++current_image){
    for (int current_px = 0; current_px < sample_nb; ++current_px){
      int currentZ = images[current_image](pixels[current_px](0), pixels[current_px](1));
      A(current_px + current_image*sample_nb, 256 + current_px) = - zWeight[currentZ];
    }
  }

  std::cout<<"partie inf gauche"<<std::endl;
  //partie inférieure gauche
  for (int i = 0; i < 254; ++i){
    A(sample_nb*image_nb + i, i) = lambda * zWeight[i+1];
    A(sample_nb*image_nb + i, i+1) = -2 * lambda * zWeight[i+1];
    A(sample_nb*image_nb + i, i+2) = lambda * zWeight[i+1];
  }

  A(sample_nb*image_nb+254, 127) = 1;


  /*********************************/
  /* Vecteur x */
  /*********************************/
  std::cout<<"vecteur X creation"<<std::endl;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(256+sample_nb);

  x = (A.transpose() * A).inverse() * A.transpose() * b;

  std::cout << std::endl << std::endl << std::endl;
  return x;

}

void setWeights(std::vector<double> &zWeight, int zMin = 0, int zMax = 255)
{
   int nSize = zMax - zMin + 1;
   zWeight.resize( nSize );

   int zMid = (zMin + zMax)/2;

   for( int i = 0; i < nSize; i++) {
        int z = zMin + i;
        if( z <= zMid ) zWeight[i] =  z - zMin;
        if( z  > zMid ) zWeight[i] =  zMax - z; 
        if(zWeight[i] <= 0) zWeight[i] = 1;
   }
}


Eigen::MatrixXd calculIrradiance(const int width,
    const int height,
    const std::vector<Eigen::MatrixXi> &imagesMatrix,
    const std::vector<double> &zWeight,
    const Eigen::VectorXd &x,
    const std::vector<double> &exposure,
    double &emin,
    double &emax
){

    //--------Reconstruction-------------------------------
  //--Irradiance-----------------------------------------

  Eigen::MatrixXd imageHDR = Eigen::MatrixXd::Zero(height, width);

  double logEtotal = 0;
  double totalWeight = 0;
  double E;
  double lnE = 0;
  int weightIJ = 1;


  std::cout<<"  "<<std::endl;
  for(int i = 0; i < height; ++i){
    for(int j = 0; j < width; ++j){
      logEtotal = 0;
      totalWeight = 0;
      for (unsigned int current_img = 0; current_img < imagesMatrix.size(); ++current_img){
        int currentZ = imagesMatrix[current_img](i,j);
        if(currentZ >= 0 && currentZ <= 255){
          weightIJ = zWeight[currentZ];
        }
        else{
          weightIJ = 1;
        }
        double lnDeltaTj = log(exposure[current_img]);
        double Ei = (x(currentZ) - lnDeltaTj);

        logEtotal += weightIJ * Ei;
        totalWeight += weightIJ;
      }

      assert(totalWeight > 0.0);
      lnE = (logEtotal/totalWeight);
      if(lnE > emax)  emax = lnE;
      if(lnE < emin)  emin = lnE;

      imageHDR(i,j) = lnE;

    }
  }

  return imageHDR;
}


//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  // check arguments
  if(argc < 3){
    std::cerr << "usage : " << argv[0] << " image_1.jpg ... image_n.jpg" << std::endl;
    std::cerr << "or : " << argv[0] << " dirname/*.jpg" << std::endl;
    exit(0);
  }

//------------------------------------------
  std::vector<ImageRGB8u> images;
  std::vector<double> exposure;
  loadImages(argc, argv, images, exposure);

  std::cout<<"createMatrix"<<std::endl;

  //vecteur de matrice images (eigen)
  std::vector<Eigen::MatrixXi> imagesMatrixRed(images.size());
  std::vector<Eigen::MatrixXi> imagesMatrixGreen(images.size());
  std::vector<Eigen::MatrixXi> imagesMatrixBlue(images.size());

  for (unsigned int current_img = 0; current_img < images.size(); ++current_img){
    imagesMatrixRed[current_img] = Eigen::MatrixXi::Zero(images[current_img].height(),images[current_img].width());
    imagesMatrixGreen[current_img] = Eigen::MatrixXi::Zero(images[current_img].height(),images[current_img].width());
    imagesMatrixBlue[current_img] = Eigen::MatrixXi::Zero(images[current_img].height(),images[current_img].width());
    for (int i = 0; i < imagesMatrixRed[current_img].rows(); ++i){
      for(int j = 0; j< imagesMatrixRed[current_img].cols(); ++j){
        imagesMatrixRed[current_img](i,j)=images[current_img](j,i)[0];
        imagesMatrixGreen[current_img](i,j)=images[current_img](j,i)[1];
        imagesMatrixBlue[current_img](i,j)=images[current_img](j,i)[2];
      }
    }
  }

  std::cout<<imagesMatrixRed[0].cols()<<std::endl<<std::endl;
  std::cout<<imagesMatrixRed[0].rows()<<std::endl<<std::endl;

  float width = imagesMatrixRed[0].cols();
  float height = imagesMatrixRed[0].rows();

  int taille = 10;

  //vecteur poids
  std::vector<double> zWeight(256);
  setWeights(zWeight);

  //création du vecteur pixels--
  //...
  std::vector<Eigen::Vector2i> pixels;
  for (int i = 0; i < taille; ++i)//imagesMatrixRed[0].cols()
  {
    for (int j = 0; j < taille; ++j)//imagesMatrixRed[0].rows()
    {
      pixels.push_back(Eigen::Vector2i(height * j/taille, width * i/taille));
    }
  }

  std::cout<<"vecteur pixel créé"<<std::endl;


  //RESPONSE RECOVERY
  std::cout<<"Calcul de responseRecovery RGB...";
  Eigen::VectorXd x_red = responseRecovery(imagesMatrixRed, exposure, pixels, zWeight, 0, 255 , 1000);
  Eigen::VectorXd x_green = responseRecovery(imagesMatrixGreen, exposure, pixels, zWeight, 0, 255 , 1000);
  Eigen::VectorXd x_blue = responseRecovery(imagesMatrixBlue, exposure, pixels, zWeight, 0, 255 , 1000);
  std::cout<<"... OK"<<std::endl;





  //IRRADIANCE
  double emax = 0, emin = 255;
  std::cout << "calcul de l'image d'irradiance ..." << std::endl;
    std::cout << "imageHDR_red";
  Eigen::MatrixXd imageHDR_red = calculIrradiance(width, height, imagesMatrixRed, zWeight, x_red, exposure, emin, emax);
    std::cout << "...OK." << std::endl;
    std::cout << "imageHDR_green ...";
  Eigen::MatrixXd imageHDR_green = calculIrradiance(width, height, imagesMatrixGreen, zWeight, x_green, exposure, emin, emax);
    std::cout << "...OK." << std::endl;
    std::cout << "imageHDR_blue ...";
  Eigen::MatrixXd imageHDR_blue = calculIrradiance(width, height, imagesMatrixBlue, zWeight, x_blue, exposure, emin, emax);
    std::cout << "...OK." << std::endl;
  std::cout << "...calcul de l'image d'irradiance OK" << std::endl;



  //--Tone-mapping---------------------------------
  //Recherche équation linéaire

  std::cout << "...TONE MAPPING..." << std::endl;

  double a,b;
  int zmin=0,zmax=255;

  a = (zmax - zmin) / (emax - emin);
  b = ( zmax / (emax + emin) ) * emin;
  
  ImageRGB8u result = images[0];
  for(uint x=0; x<images[0].width(); ++x){
    for(uint y=0; y<images[0].height(); ++y){
      
      //result(x,y)[0] = a*imageHDR_red(y,x) +b;
      result(x,y)[0] = ((imageHDR_red(y,x)-emin) * (zmax-zmin)) / (emax - emin);
      result(x,y)[1] = ((imageHDR_green(y,x)-emin) * (zmax-zmin)) / (emax - emin);
      result(x,y)[2] = ((imageHDR_blue(y,x)-emin) * (zmax-zmin)) / (emax - emin);

      result(x,y)[0] = imageHDR_red(y,x)

      float v = result(x,y)[0];
      //if(x== 90 && y < 100) std::cout << v << std::endl;
    }
  }

  //std::cout << "a= " << a << " ------ b = " << b << std::endl;
  std::cout << "emax = " << emax << " -----  emin = " << emin << std::endl;

  std::cout<<"END_"<<std::endl;

  // save the final image
  //saveJPG(result,"output/Tango-Charly_3.00.jpg");
  saveJPG(result,"output/escalier_v1.47.jpg");

  return 0;
}