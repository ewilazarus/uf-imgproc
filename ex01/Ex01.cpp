/**
 * Image Processing and Computer Graphics
 *
 * Exercise 1: Noise, basic operators and filters
 */

#include <cstdlib>   /// contains: EXIT_SUCCESS
#include <iostream>  /// contains: std::cout etc.
#include "CMatrix.h"
#include "CTensor.h"
#include "Noise.h"
#include "Filter.h"

#define IMG_COUNT 50

void noise(void)
{
  /// Define array of images
  CMatrix<float> images[IMG_COUNT + 1];

  /// Read image from a PGM file
  images[0].readFromPGM("lena.pgm");

  /// Instantiate noise producer
  NoiseProducer nProducer(images[0]);

  /// Define avgImage
  CMatrix<float> avgImage(images[0].xSize(), images[0].ySize(), 0.0);


  /// Populate images with noise using variance equal to index times 10
  for (int i = 1; i <= IMG_COUNT; i++) {
  	  images[i] = nProducer.applyGaussianNoise(0, i * 10);
  	  avgImage += images[i];
  }
  avgImage *= (1.0/IMG_COUNT);

  float PSNR;

  // Calculate PSNR with variance of 10
  PSNR = calculatePSNR(images[0], images[1]);
  printf("PSNR (variance = 10): %f\n", PSNR);
  images[1].writeToPGM("lenaNoisy10.pgm");

  // Calculate PSNR with variance of 20
  PSNR = calculatePSNR(images[0], images[2]);
  printf("PSNR (variance = 20): %f\n", PSNR);
  images[2].writeToPGM("lenaNoisy20.pgm");

  // Calculate PSNR of averaged image
  PSNR = calculatePSNR(images[0], avgImage);
  printf("PSNR of averaged image: %f\n", PSNR);
  avgImage.writeToPGM("lenaNoisyAvg.pgm");
}

void difference(void)
{
	CTensor<float> sidenbladh;
	CTensor<float> sidenbladhBG;

	sidenbladh.readFromPPM("Sidenbladh.ppm");
	sidenbladhBG.readFromPPM("SidenbladhBG.ppm");

	//sidenbladhBG *= -1;
	//sidenbladh += sidenbladhBG;
	//sidenblad
	//sidenbladh.writeToPPM("SidenbladhDiff.ppm");
	
	int xSize = sidenbladh.xSize();
	int ySize = sidenbladh.ySize();
	int zSize = sidenbladh.zSize();
	
	CTensor<float> sidenbladhDiff(xSize, ySize, zSize);

	for (int x = 0; x < xSize; x++) {
		for (int y = 0; y < ySize; y++) {
			for (int z = 0; z < zSize; z++) {
				//sidenbladhDiff(x, y, z) = sidenbladh(x, y, z) - sidenbladhBG(x, y, z);
				sidenbladhDiff(x, y, z) = fabs(sidenbladh(x, y, z) - sidenbladhBG(x, y, z));
			}
		}
	}
	sidenbladhDiff.writeToPPM("SidenbladhDiff.ppm");
}

void filter(void)
{
	CMatrix<float> chinaToilet;
	CMatrix<float> *chinaToiletFilterG;
	CMatrix<float> *chinaToiletFilterB;
	CMatrix<float> *chinaToiletFilterR;

	chinaToilet.readFromPGM("chinaToilet.pgm");

	chinaToiletFilterG = gaussianFilter(chinaToilet, 5.0);
	chinaToiletFilterG->writeToPGM("chinaToiletFilterG.pgm");

	chinaToiletFilterB = boxFilter(chinaToilet, 5.0);
	chinaToiletFilterB = boxFilter(*chinaToiletFilterB, 5.0);
	chinaToiletFilterB = boxFilter(*chinaToiletFilterB, 5.0);
	chinaToiletFilterB->writeToPGM("chinaToiletFilterB.pgm");

	chinaToiletFilterR = recursiveFilter(chinaToilet, 5.0);
	//chinaToiletFilterR = recursiveFilter(chinaToiletFilterR, 1.0);
	//chinaToiletFilterR = recursiveFilter(chinaToiletFilterR, 1.0);
	chinaToiletFilterR->writeToPGM("chinaToiletFilterR.pgm");


	CTensor<float> fallingMangoes;
	CTensor<float> *fallingMangoesFilter;
	CMatrix<float> *tempFallingMangoes;

	fallingMangoes.readFromPPM("fallingMangoes.ppm");
	fallingMangoesFilter = new CTensor<float>(fallingMangoes.xSize(), fallingMangoes.ySize(), 3);

	for (int i = 0; i < 3; ++i) {
		tempFallingMangoes = gaussianFilter(fallingMangoes.getMatrix(i), 5.0);
		fallingMangoesFilter->putMatrix(*tempFallingMangoes, i);
	}
	fallingMangoesFilter->writeToPPM("fallingMangoesFilter.ppm");
}

int main(int argc, char** args)
{
  /// Tell the compiler not to throw warnings for unused variables
  /// Remove these lines if you want to use command line arguments.
  (void)argc;
  (void)args;

  noise();
  difference();
  filter();

  return EXIT_SUCCESS;
}

