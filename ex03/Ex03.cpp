/**
 * Image Processing and Computer Graphics
 *
 * Exercise 3: Local descriptors
 */

#include <cstdlib>   /// contains: EXIT_SUCCESS
#include <iostream>  /// contains: std::cout etc.
#include "CMatrix.h"
#include "CTensor.h"
#include "CFilter.h"

#define LAYER_COUNT 3


void computeBhaskara(float *roots, float a, float b, float c)
{
	// r = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}
	float delta = b*b - 4*a*c;
	roots[0] = (-b + sqrt(delta))/(2*a);
	roots[1] = (-b - sqrt(delta))/(2*a);
}

float computeSecondLargestEigenVal(CMatrix<float> J)
{
	// compute weights a, b and c where: ax^2 + bx + c = 0
	float a, b, c;
	a = 1.0;
	b = abs(J(0, 0) + J(1, 1));
	c = -1*J(0, 1)*J(1, 0) + J(0, 0)*J(1, 1);

	// compute roots
	float eigenVals[2];
	computeBhaskara(eigenVals, a, b, c);

	// return the "second largest" (smaller?) eigenvalue
	return (eigenVals[0] < eigenVals[1]) ? eigenVals[1] : eigenVals[0];
}

void applyFilterToLayers(CFilter<float> f, CTensor<float> image)
{
	// apply the filter f to every image layer
	for (int l = 0; l < LAYER_COUNT; l++) {
		CMatrix<float> layer = image.getMatrix(l);
		NFilter::filter(layer, f, f);
		image.putMatrix(layer, l);
	}
}

bool isPeak(CMatrix<float> M, int x, int y)
{
	float u = -1, d = -1, l = -1, r = -1;
	float px = M(x, y);

	if (x - 1 < M.xSize()) u = M(x - 1, y);
	if (x + 1 >= 0) d = M(x + 1, y);
	if (y - 1 >= 0) l = M(x, y - 1);
	if (y + 1 < M.ySize()) r = M(x, y + 1);

	return px > u && px > d && px > l && px > r;
}

void cornerDetector(char* imagePath, int k, int s, int p, int d1, int d2)
{
	int kernelSize = k;

	CTensor<float> inImage;
	inImage.readFromPPM(imagePath);

	CTensor<float> outImage(inImage);

	// smooths inImage
	CSmooth<float> smooth(s, p);  // CSmooth(float aSigma, float aPrecision)
	applyFilterToLayers(smooth, inImage);

	// compute derivatives
	CDerivative<float> derivative1(d1);
	CMatrix<float> sum(inImage.xSize(), inImage.ySize(), 0.0);
	CMatrix<float> Ix(inImage.xSize(), inImage.ySize());
	CMatrix<float> Iy(inImage.xSize(), inImage.ySize());
	for (int l = 0; l < 3; l++) {
		sum += inImage.getMatrix(l); 
	}
	NFilter::filter(sum, Ix, derivative1, 1);
	NFilter::filter(sum, Iy, 1, derivative1);

	// initialize matrix to hold the "cornerness" of each pixel
	CMatrix<float> cornerness(inImage.xSize(), inImage.ySize(), 0.0);

	#pragma omp parallel for
	for (int x = kernelSize; x < inImage.xSize() - kernelSize; x++) {
		for (int y = kernelSize; y < inImage.ySize() - kernelSize; y++) {

			// initialize J with zeros and convolve it with a box-window
			CMatrix<float> J(2, 2, 0);
			for (int i = x - kernelSize; i < x + kernelSize; i++) {
				for (int j = y - kernelSize; j < y + kernelSize; j++) {
					for (int k = 0; k < LAYER_COUNT; k++) {
						J(0, 0) += Ix(i, j)*Ix(i, j);
						J(0, 1) += Ix(i, j)*Iy(i, j);
						J(1, 0) += Ix(i, j)*Iy(i, j);
						J(1, 1) += Iy(i, j)*Iy(i, j);
					}
				}
			}
			cornerness(x, y) = computeSecondLargestEigenVal(J);
		}
	}
	CDerivative<float> derivative2(d2);
	CMatrix<float> cDerivatives(cornerness);
	NFilter::filter(cornerness, derivative2, derivative2);

	#pragma omp parallel for
	for (int x = 0; x < inImage.xSize(); x++) {
		for (int y = 0; y < inImage.ySize(); y++) {
			if (!cornerness(x, y) && isPeak(cornerness, x, y)) {
				outImage.drawRect(x - 1, y - 1, x + 1, y + 1);
			}
		}
	}
	
	std::string outImagePath = "Cornered" + std::string(imagePath);
	outImage.writeToPPM(outImagePath.c_str());
}

float computeGradientDirectrion(float partialX, float partialY)
{
	return atan(partialY / partialX);
}


void SIFT()
{
	CTensor<float> inImage;
	inImage.readFromPPM("tennis500.ppm");

	CTensor<float> outImage(inImage);

	// smooths inImage
	CSmooth<float> smooth(3, 3);  // CSmooth(float aSigma, float aPrecision)
	applyFilterToLayers(smooth, inImage);

	// compute derivatives
	CDerivative<float> derivative(2);
	CMatrix<float> sum(inImage.xSize(), inImage.ySize(), 0.0);
	CMatrix<float> Ix(inImage.xSize(), inImage.ySize());
	CMatrix<float> Iy(inImage.xSize(), inImage.ySize());
	for (int l = 0; l < 3; l++) {
		sum += inImage.getMatrix(l); 
	}

	NFilter::filter(sum, Ix, derivative, 1);
	NFilter::filter(sum, Iy, 1, derivative);
	
	printf("%f\n", computeGradientDirectrion(23.0, 13.0));
}

int main(int argc, char** args) 
{  
	int values[7];
	for (int i = 2; i < argc; i++) {
		values[i] = atoi(args[i]);
	}

	cornerDetector(args[1],
			values[2],
			values[3],
			values[4],
			values[5],
			values[6]);

	/* SIFT(); */

	return EXIT_SUCCESS;
}

