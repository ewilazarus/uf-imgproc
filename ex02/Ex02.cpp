/**
 * Image Processing and Computer Graphics
 *
 * Exercise 2: Motion estimation
 */

#include <cstdlib>   /// contains: EXIT_SUCCESS
#include <iostream>  /// contains: std::cout etc.
#include <algorithm>
#include <Eigen/Sparse>
//#include <Eigen/IterativeLinearSolvers>
//#include <Eigen/SparseLU>
//#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <vector>
#include "CMatrix.h"
#include "CTensor.h"
#include "CFilter.h"
#include "flowToImage.h"

#define U 0
#define V 1
#define ITER 1000

CMatrix<float> _invert(CMatrix<float> M)
{
	CMatrix<float> I(M.xSize(), M.ySize());
	I(0, 0) = M(1, 1);
	I(0, 1) = -M(0, 1);
	I(1, 0) = -M(1, 0);
	I(1, 1) = M(0, 0);
	
	float det = M(0, 0)*M(1, 1) - M(0, 1)*M(1, 0);
	I *= (1 / det);

	return I;
}

void LucasKanade() 
{
	// Define kernel size (?)
	int range = 7;

	// Read images from file
	CMatrix<float> street1;
	street1.readFromPGM("cropped-street_000009.pgm");

	CMatrix<float> street2;
	street2.readFromPGM("cropped-street_000010.pgm"); 
	// Creates "smoother"
	CSmooth<float> smoother(3, 3);

	// Smooths input images
	NFilter::filter(street1, smoother, smoother);
	NFilter::filter(street2, smoother, smoother);

	// Specify x, y, t sizes
	int xSize = street1.xSize();
	int ySize = street1.ySize();
	int tSize = 2;

	// Creates CTensor to hold the two images loaded before
	CTensor<float> street(xSize, ySize, tSize);
	street.putMatrix(street1, 0);
	street.putMatrix(street2, 1);

	// Create derivative filter
	CDerivative<float> derivative(5);	

	// Create CTensors to hold the x, y, t derivatives
	CTensor<float> Ix(xSize, ySize, tSize);
	CTensor<float> Iy(xSize, ySize, tSize);
	CTensor<float> It(xSize, ySize, tSize);

	// Compute x, y, t derivatives
	NFilter::filter(street, Ix, derivative, 1, 1);
	NFilter::filter(street, Iy, 1, derivative, 1);
	NFilter::filter(street, It, 1, 1, derivative);

	// Creates a CTensor to hold the resulting flow (z = 2 => "u and v")
	CTensor<float> flow(xSize, ySize, 2);

	// Computes vectors U and V for every pixel (x, y)
	for (int x = range; x < xSize - range; x++) {
		for (int y = range; y < ySize - range; y++) {

			CMatrix<float> A(2, 2, 0);
			CMatrix<float> B(2, 1, 0);
			
			for (int i = x - range; i < x + range; i++) {
				for (int j = y - range; j < y + range; j++) {
					A(0, 0) += Ix(i, j, 0)*Ix(i, j, 0);
					A(0, 1) += Ix(i, j, 0)*Iy(i, j, 0);
					A(1, 0) += Ix(i, j, 0)*Iy(i, j, 0);
					A(1, 1) += Iy(i, j, 0)*Iy(i, j, 0);

					B(0, 0) += Ix(i, j, 0)*It(i, j, 1);
					B(1, 0) += Iy(i, j, 0)*It(i, j, 1);
				}
			}

			B *= -1;
			CMatrix<float> Ainv = _invert(A);

			flow(x, y, U) = Ainv(0, 0)*B(0, 0) + Ainv(0, 1)*B(1, 0);
			flow(x, y, V) = Ainv(1, 0)*B(0, 0) + Ainv(1, 1)*B(1, 0);
		}
	}

	// Creates a RGB empty image
	CTensor<float> flowResult(xSize, ySize, 3);

	// Depicts flow
	flowToImage(flow, flowResult);
	
	// Outputs result
	flowResult.writeToPPM("croppedStreetMotionLK.ppm");
}

typedef Eigen::SparseMatrix<float> MyMatrix;
typedef Eigen::VectorXf MyVector;
typedef Eigen::Triplet<float> T;

int get1DCoord(int x, int y, int size)
{
	return x*size + y;
}

void HornSchunck()
{
	// Read images from file
	CMatrix<float> street1;
	street1.readFromPGM("cropped-street_000009.pgm");

	CMatrix<float> street2;
	street2.readFromPGM("cropped-street_000010.pgm");

	// Creates "smoother"
	CSmooth<float> smoother(3, 3);

	// Smooths input images
	NFilter::filter(street1, smoother, smoother);
	NFilter::filter(street2, smoother, smoother);

	// Specify x, y, t sizes
	int size = street1.size();
	int xSize = street1.xSize();
	int ySize = street1.ySize();
	int tSize = 2;

	std::cout << "Matrix: " << xSize << "X" << ySize << "\n";

	// Creates CTensor to hold the two images loaded before
	CTensor<float> street(xSize, ySize, tSize);
	street.putMatrix(street1, 0);
	street.putMatrix(street2, 1);

	// Create derivative filter
	CDerivative<float> derivative(5);	

	// Create CTensors to hold the x, y, t derivatives
	CTensor<float> Ix(xSize, ySize, tSize);
	CTensor<float> Iy(xSize, ySize, tSize);
	CTensor<float> It(xSize, ySize, tSize);

	// Compute x, y, t derivatives
	NFilter::filter(street, Ix, derivative, 1, 1);
	NFilter::filter(street, Iy, 1, derivative, 1);
	NFilter::filter(street, It, 1, 1, derivative);

	// Initializes w with 0, this vector is going to hold the answer for u and v
	MyVector w(2*size);
	w.setOnes();
	float alpha = 1;

	for (int iter = 0; iter < ITER; iter++) {
		std::vector<T> coordinateS;
		std::vector<T> coordinateD;
		MyVector b(2*size);

		// Iterate over the image twice to fill the coordinates needed in
		// the matrix A and the vector b of our system. (Ax=b)
		for (int k = 0; k < 2; k++) {
			for (int x = 0; x < xSize; x++) {
				for (int y = 0; y < ySize; y++) {

					int offset = k*size;
					int i = get1DCoord(x, y, ySize) + offset;
					int j;

					coordinateS.push_back(T(i, i, (-4)*w[i]));
					/* printf("\n'S' Processing (x: %d, y: %d) = (i:%d, j:%d): %f ...\n", x, y, i, i, (-4)*w[i]); */

					if (x + 1 < xSize) {
						j = get1DCoord(x + 1, y, ySize) + offset;
						coordinateS.push_back(T(i, j, w[j]));
						/* printf("'S' Processing (x + 1: %d, y: %d) = (i:%d, j:%d): %f ...\n", x+1, y, i, j, w[j]); */
					}
					if (x - 1 >= 0) {
						j = get1DCoord(x - 1, y, ySize) + offset;
						coordinateS.push_back(T(i, j, w[j]));
						/* printf("'S' Processing (x - 1: %d, y: %d) = (i:%d, j:%d): %f ...\n", x-1, y, i, j, w[j]); */
					}
					if (y + 1 < ySize) {
						j = get1DCoord(x, y + 1, ySize) + offset;
						coordinateS.push_back(T(i, j, w[j]));
						/* printf("'S' Processing (x: %d, y + 1: %d) = (i:%d, j:%d): %f ...\n", x, y+1, i, j, w[j]); */
					}
					if (y - 1 >= 0) {
						j = get1DCoord(x, y - 1, ySize) + offset;
						coordinateS.push_back(T(i, j, w[j]));
						/* printf("'S' Processing (x: %d, y - 1: %d) = (i: %d, j:%d): %f ...\n", x, y-1, i, j, w[j]); */

					}

					if (k == 0) {
						/* printf("\n'D' Processing (x: %d, y: %d) = (i: %d, j: %d): %f\n", x, y, i, i, -(Ix(x, y)*Ix(x, y)*w[i])/alpha); */
						coordinateD.push_back(T(i, i, -(Ix(x, y)*Ix(x, y)*w[i])/alpha));
						/* printf("'D' Processing (x: %d, y: %d) = (i: %d, j: %d): %f\n", x, y, i, i + size, -(Ix(x, y)*Iy(x, y)*w[i + size])/alpha); */
						coordinateD.push_back(T(i, i + size, -(Ix(x, y)*Iy(x, y)*w[i + size])/alpha));
						/* printf("\n'b' Processing (x: %d, y: %d): %f\n", x, y, Ix(x, y)*It(x, y)/alpha); */
						b[i] = Ix(x, y)*It(x, y);
					} else {
						/* printf("\n'D' Processing (x: %d, y: %d) = (i: %d, j: %d): %f\n", x, y, i, i - size, -(Ix(x, y)*Iy(x, y)*w[i - size])/alpha); */
						coordinateD.push_back(T(i, i - size, -(Ix(x, y)*Iy(x, y)*w[i - size])/alpha));
						/* printf("'D' Processing (x: %d, y: %d) = (i: %d, j: %d): %f\n", x, y, i, i, -(Iy(x, y)*Iy(x, y)*w[i])/alpha); */
						coordinateD.push_back(T(i, i, -(Iy(x, y)*Iy(x, y)*w[i])/alpha));
						/* printf("\n'b' Processing (x: %d, y: %d): %f\n", x, y, Iy(x, y)*It(x, y)/alpha); */
						b[i] = Iy(x, y)*It(x, y)/alpha;
					}
				}
			}
			/* int tmp; */
			/* std::cin >> tmp; */
		}

		// Adds matrices for smoothing (As) and data (Ad)
		MyMatrix As(2*size, 2*size);
		As.setFromTriplets(coordinateS.begin(), coordinateS.end());

		MyMatrix Ad(2*size, 2*size);
		Ad.setFromTriplets(coordinateD.begin(), coordinateD.end());

		MyMatrix A(2*size, 2*size);
		A = As + Ad;

		// Solve
		Eigen::SimplicialLDLT<MyMatrix> solver;
		solver.compute(A);
		if (solver.info()!=1) {
			printf("Failed to compute.\n");
		}
		w = solver.solve(b);
		if (solver.info()!=1) {
			printf("Failed to solve.\n");
		}
	}

	// Creates a CTensor to hold the resulting flow (z = 2 => "u and v")
	CTensor<float> flow(xSize, ySize, 2);
	
	// Populate the flow with the computed values
	for (int z = 0; z < 2; z++) {
		int offset = z*size;
		for (int x = 0; x < xSize; x++) {
			for (int y = 0; y < ySize; y++) {
				float value = w[get1DCoord(x, y, ySize) + offset];
				flow(x, y, z) = value;
			}
		}
	}

	// Creates a RGB empty image
	CTensor<float> flowResult(xSize, ySize, 3);

	// Depicts flow
	flowToImage(flow, flowResult);
	
	// Outputs result
	flowResult.writeToPPM("croppedStreetMotionHS.ppm");
}


int main(int argc, char** args) 
{  
  /// Tell the compiler not to throw warnings for unused variables
  /// Remove these lines if you want to use command line arguments.
  (void)argc;
  (void)args;

  //LucasKanade();
  HornSchunck();

  return EXIT_SUCCESS;
}

