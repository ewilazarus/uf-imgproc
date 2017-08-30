Exercises for the class of Image Processing & CG @ Uni-Freiburg
===============================================================

Questions
~~~~~~~~~

Exercise 1
----------

1. What am I supposed to see with the average PSNR of the image sequence? I'm currently seeing a black image
   output image values for one px every iteration

2. When filtering I get two "lanes" in the bottom of the images, one gray for the Gaussian and one black for the Box. What's going on?

3. Posted at the class' forum

Exercise 2
----------

1. Lucas-Kanade:

- So far, so good? I took the derivatives Ix, Iy and It (using CTensor, CDerivative, NFilter::filter), from now on, whenever I need them, is it going to be I(x, y, 0) or I(x, y, 1) [or I(x, y) - bilinear interpolation]?

- What is the aSize parameter of the CDerivative constructor?

- What is the aPrecision parameter of the CSmooth constructor? It doesn't seem to affect my end result

- Will I get better results if I use the Gaussian convolution instead of the sum to calculate vectors (u, v)?

2. Horn-Schunck:

- What is h in (A^h)*w=b? Is it the grid finess?
