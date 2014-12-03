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
/*TEST 010 */

////////////////////////////////////////////////////////////////////////////////////////
////////////////////RESPONSE RECOVERY//////////////////////////////////////////////////

int getWeight(int z, int zMin, int zMax){
  if( (double)z <= (zMin + zMax)/2.0){
    return z - zMin;
  }
  else{
    return zMax - z;
  }
}


Eigen::VectorXd responseRecovery (const std::vector <Eigen::MatrixXi> &images ,
const std :: vector < double > & exposure ,
const std :: vector < Eigen::Vector2i > & pixels ,
const int valueMin ,
const int valueMax ,
const double lambda ){

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
      b(current_image * sample_nb + current_px) = getWeight(images[current_image](pixels[current_px](0), pixels[current_px](1)), 0, 255) * log(exposure[current_image]); //* calcul_poids(pixel(current_image,current_px))
    }
  }
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
      int z = images[current_image](pixels[current_px](0),pixels[current_px](1));//TO DO
      A(current_px + current_image*sample_nb, z) = getWeight(images[current_image](pixels[current_px](0), pixels[current_px](1)), 0, 255); //ajouter poids
    }
  }

  std::cout<<"part sup droite"<<std::endl;
  //partie supérieure droite
  for (int current_image = 0; current_image < image_nb; ++current_image){
    for (int current_px = 0; current_px < sample_nb; ++current_px){
      A(current_px + current_image*sample_nb, 256 + current_px) = -getWeight(images[current_image](pixels[current_px](0), pixels[current_px](1)), 0, 255); //ajouter poids
    }
  }

  std::cout<<"partie inf gauche"<<std::endl;
  //partie inférieure gauche
  for (int i = 0; i < 254; ++i){
    A(sample_nb*image_nb + i, i) = lambda * getWeight(i+1, 0, 255);
    A(sample_nb*image_nb + i, i+1) = -2 * lambda * getWeight(i+1, 0, 255);
    A(sample_nb*image_nb + i, i+2) = lambda * getWeight(i+1, 0, 255);
  }

  A(sample_nb*image_nb+254, 128) = 1;


  /*********************************/
  /* Vecteur x */
  /*********************************/
  std::cout<<"vecteur X creation"<<std::endl;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(256+sample_nb);

  x = (A.transpose() * A).inverse() * A.transpose() * b;

  std::cout << std::endl << std::endl << std::endl;
  //std::cout<<"b: " << b<<std::endl<<std::endl<<std::endl<<std::endl;
  std::cout<<"x"<<x<<std::endl;

  return x;

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

  for (int current_img = 0; current_img < images.size(); ++current_img){
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

  //création du vecteur pixels--TO DO!!
  //...
  std::vector<Eigen::Vector2i> pixels;
  for (int i = 0; i < 10; ++i)//imagesMatrixRed[0].cols()
  {
    for (int j = 0; j < 10; ++j)//imagesMatrixRed[0].rows()
    {
      pixels.push_back(Eigen::Vector2i(i,j));
    }
  }

  std::cout<<"OK"<<std::endl;

  std::cout<<"create A"<<std::endl;
  responseRecovery(imagesMatrixRed, exposure, pixels, 0, 255 , 1000);
  std::cout<<"OK"<<std::endl;





  // draw something
  /*for(uint x=images[0].width()/4; x<images[0].width()*3/4.0; ++x)
  for(uint y=images[0].height()/4; y<images[0].height()*3/4.0; ++y){
  images[0](x,y)[0] = 255; // R
  images[0](x,y)[1] = 0; // G
  images[0](x,y)[2] = 0; // B
  }*/

  // save the final image
  //saveJPG(images[0],"output/test465656.jpg");

  return 0;
}