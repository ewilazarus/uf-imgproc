#include <math.h>
#include "CMatrix.h"
//#include "ImageDisplay.h"

#define X(a) _fixBoundaries(a, xSize)
#define Y(a) _fixBoundaries(a, ySize)


/**************************** HELPER FUNCTIONS *******************************/

int _fixBoundaries(int index, int axisSize)
{
	int rv = index;
	if (index < 0) {
		rv = -index;	
	} else if (index > axisSize - 1) {
		rv = 2*axisSize - index - 2;
	}
	return rv;
}

void _computeKernelIxs(int *kIxs, int const size) {
	for (int i = 0, s = -(size/2); i < size; i++, s++) {
		if (s == 0) s++;
		kIxs[i] = s;	
	}
}


/**************************** GAUSSIAN FILTER ********************************/

void _computeGKernel(float const sigma, int const size, int * const bounds, float *kernel) 
{
	float sum = 0.0;
	for (int i = 0, x = bounds[i]; i < size; i++) {
		float s = exp(-((x*x) / (2*sigma*sigma))) / (sqrt(2*M_PI)*sigma);
		kernel[i] = s;
		sum += s;
	}
	for (int i = 0; i < size; i++) {  // Normalization
		kernel[i] /= sum;
	}
}

CMatrix<float> *gaussianFilter(CMatrix<float> image, float sigma)
{
	CMatrix<float> *filteredImage = new CMatrix<float>(image);

	int kSize = 6*sigma;

	int *kIxs = (int *)malloc(kSize*sizeof(int));
	_computeKernelIxs(kIxs, kSize);

	float *kernel = (float *)malloc(kSize*sizeof(float));
	_computeGKernel(sigma, kSize, kIxs, kernel);

	int xSize = filteredImage->xSize();
	int ySize = filteredImage->ySize();

	for(int x = 0; x < xSize; x++) {
		for (int y = 0; y < ySize; y++) {
			float pxout = 0.0;
			for (int k = 0; k < kSize; k++) {
				pxout += (*filteredImage)(X(x + kIxs[k]), y)*kernel[k];
			}
			(*filteredImage)(x, y) = pxout;	
		}
	}		
	for (int y = 0; y < ySize; y++) {
		for(int x = 0; x < xSize; x++) {
			float pxout = 0.0;
			for (int k = 0; k < kSize; k++) {
				pxout += (*filteredImage)(x, Y(y + kIxs[k]))*kernel[k];
			}
			(*filteredImage)(x, y) = pxout;	
		}
	}		
	free(kIxs);
	free(kernel);
	return filteredImage;
}


/******************************* BOX FILTER **********************************/

CMatrix<float> *boxFilter(CMatrix<float> image, float sigma)
{
	CMatrix<float> *filteredImage = new CMatrix<float>(image);

	int kSize = 2*sigma;

	int *kIxs = (int *)malloc(kSize*sizeof(int));
	_computeKernelIxs(kIxs, kSize);

	float kernel = 1.0/kSize;

	int xSize = filteredImage->xSize();
	int ySize = filteredImage->ySize();
	for(int x = 0; x < xSize; x++) {
		float pxout = 0.0;
		for (int k = 0; k < kSize; k++) {
			pxout += (*filteredImage)(X(x + kIxs[k]), 0)*kernel;
		}
		(*filteredImage)(x, 0) = pxout;
		for (int y = 1; y < ySize; y++) {
			pxout = (*filteredImage)(x, y - 1);
			pxout -= (*filteredImage)(x, Y(y - sigma))*kernel;
			pxout += (*filteredImage)(x, Y(y + sigma))*kernel;
			(*filteredImage)(x, y) = pxout;
		}
	}		
	for (int y = 0; y < ySize; y++) {
		float pxout = 0.0;
		for (int k = 0; k < kSize; k++) {
			pxout += (*filteredImage)(0, Y(y + kIxs[k]))*kernel;
		}
		(*filteredImage)(0, y) = pxout;
		for(int x = 1; x < xSize; x++) {
			pxout = (*filteredImage)(x - 1, y);
			pxout -= (*filteredImage)(X(x - sigma), y)*kernel;
			pxout += (*filteredImage)(X(x + sigma), y)*kernel;
			(*filteredImage)(x, y) = pxout;
		}
	}		
	free(kIxs);
	return filteredImage;
}


/**************************** RECURSIVE FILTER *******************************/

int xSize, ySize;
float alpha, eAlpha, e2Alpha, N;
CMatrix<float> *helperF, *helperG, *originalImg, *filteredImg;

void _setup(float sigma) {
	alpha = 5.0 / sqrt(4*M_PI*sigma*sigma);
	eAlpha = exp(-alpha);
	e2Alpha = exp(-2*alpha);
	N = pow(1 - eAlpha, 2) / (1 + 2*alpha*eAlpha - e2Alpha);
}

void _teardown()
{
	delete helperF;
	delete helperG;
	delete originalImg;
}

float _rfX(int x, int y)
{
	if (y < 0) return 0.0;
	if ((*helperF)(x, y) >= 0) return (*helperF)(x, y);
	//printf("Calculating F (%d, %d)\n", x, y);

	float i = (*originalImg)(x, y);
	float iM1 = (y - 1) < 0 ? 0.0 : (*originalImg)(x, y - 1);
	float fIM1 = (y - 1) < 0 ? 0.0 : (*helperF)(x, y - 1) = _rfX(x, y - 1);
	float fIM2 = (y - 2) < 0 ? 0.0 : (*helperF)(x, y - 2) = _rfX(x, y - 2);

	return N * (i + eAlpha*(alpha - 1)*iM1) + 2*eAlpha*fIM1 - e2Alpha*fIM2;
}

float _rgX(int x, int y)
{
	if (y > ySize - 1) return 0.0;
	if ((*helperG)(x, y) >= 0) return (*helperG)(x, y);
	//printf("Calculating G (%d, %d)\n", x, y);

	float iP1 = (y + 1) > (ySize - 1) ? 0.0 : (*originalImg)(x, y + 1);
	float iP2 = (y + 2) > (ySize - 1) ? 0.0 : (*originalImg)(x, y + 2);
	float gIP1 = (y + 1) > (ySize - 1) ? 0.0 : (*helperG)(x, y + 1) = _rgX(x, y + 1);
	float gIP2 = (y + 2) > (ySize - 1) ? 0.0 : (*helperG)(x, y + 2) = _rgX(x, y + 2);

	return N * (eAlpha*(alpha + 1)*iP1 - e2Alpha*iP2) + 2*eAlpha*gIP1 - e2Alpha*gIP2;
}

float _rfY(int x, int y)
{
	if (x < 0) return 0.0;
	if ((*helperF)(x, y) >= 0) return (*helperF)(x, y);
	//printf("Calculating F (%d, %d)\n", x, y);

	float i = (*originalImg)(x, y);
	float iM1 = (x - 1) < 0 ? 0.0 : (*originalImg)(x - 1, y);
	float fIM1 = (x - 1) < 0 ? 0.0 : (*helperF)(x - 1, y) = _rfY(x - 1, y);
	float fIM2 = (x - 2) < 0 ? 0.0 : (*helperF)(x - 2, y) = _rfY(x - 2, y);

	return N * (i + eAlpha*(alpha - 1)*iM1) + 2*eAlpha*fIM1 - e2Alpha*fIM2;
}

float _rgY(int x, int y)
{
	if (x > xSize - 1) return 0.0;
	if ((*helperG)(x, y) >= 0) return (*helperG)(x, y);
	//printf("Calculating G (%d, %d)\n", x, y);

	float iP1 = (x + 1) > (xSize - 1) ? 0.0 : (*originalImg)(x + 1, y);
	float iP2 = (x + 2) > (xSize - 1) ? 0 : (*originalImg)(x + 2, y);
	float gIP1 = (x + 1) > (xSize - 1) ? 0.0 : (*helperG)(x + 1, y) = _rgY(x + 1, y);
	float gIP2 = (x + 2) > (xSize - 1) ? 0.0 : (*helperG)(x + 2, y) = _rgY(x + 2, y);

	return N * (eAlpha*(alpha + 1)*iP1 - e2Alpha*iP2) + 2*eAlpha*gIP1 - e2Alpha*gIP2;
}

CMatrix<float> *recursiveFilter(CMatrix<float> image, float sigma)
{
	filteredImg = new CMatrix<float>(image);

	_setup(sigma);

	xSize = image.xSize();
	ySize = image.ySize();

	originalImg = new CMatrix<float>(image);
	helperF = new CMatrix<float>(xSize, ySize, -1);
	helperG = new CMatrix<float>(xSize, ySize, -1);

	for (int y = 0; y < ySize; y++) {
		for (int x = 0; x < xSize; x++) {
			(*filteredImg)(x, y) = _rfX(x, y) + _rgX(x, y);
		}
	}

	*helperF = -1;
	*helperG = -1;
	delete originalImg;
	originalImg = new CMatrix<float>(*filteredImg);

	for (int x = 0; x < xSize; x++) {
		for (int y = 0; y < ySize; y++) {
			(*filteredImg)(x, y) = _rfY(x, y) + _rgY(x, y);
		}
	}

	_teardown();

	return filteredImg;
}
