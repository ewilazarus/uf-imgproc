#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "NMath.h"
#include "CMatrix.h"


class NoiseProducer {
	private:
		CMatrix<float> templateImage;
		float U;
		float V;
		bool switcher;

		float _clip(float px)
		{
			if (px < 0) {
				px = 0;
			} else if (px > 255) {
				px = 255;
			}

			return px;
		}

		float _generateGaussianNoise(float mean, float variance)
		{
			// Compute Box-Muller
			float M = sqrt((-2) * log(U)) * sin(2 * M_PI * V);
			float N = sqrt((-2) * log(U)) * cos(2 * M_PI * V);

			float noise = (switcher = !switcher) ? M : N;
			return noise * variance + mean;
		}

		float _applyGNoise(float mean, float variance, float px) {
			return _clip(px + _generateGaussianNoise(mean, variance));
		}

	public:
		NoiseProducer(CMatrix<float> image)
		{
			srand((unsigned)time(NULL));
			templateImage = image;
			switcher = 0;
		}

		CMatrix<float> applyGaussianNoise(float mean, float variance)
		{
			CMatrix<float> noisedImage(templateImage);

			// Iterate over pixels
			for (int y = 0; y < templateImage.ySize(); ++y) {
				for (int x = 0; x < templateImage.xSize(); ++x) {
					// Create two independent random variables
					U = NMath::random();
					V = NMath::random();
					
					float px = templateImage(x, y);
					noisedImage(x, y) = _applyGNoise(mean, variance, px);
				}
			}
			return noisedImage;
		}
};


float calculatePSNR(CMatrix<float> originalImage, CMatrix<float> noisedImage)
{
	int N = originalImage.size();

	float n = N * pow((noisedImage.max() - noisedImage.min()), 2);
	float d = 0;

	CMatrix<float> diffImage = noisedImage - originalImage;
	for (int y = 0; y < diffImage.ySize(); ++y) {
		for (int x = 0; x < diffImage.xSize(); ++x) {
			d += pow(diffImage(x, y), 2);
		}
	}

	return 10*log10(n / d);
}
