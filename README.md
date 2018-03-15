# winograd minimal filtering for Convolutional Neural Network
I'm going to use the Winograd’s minimal ﬁltering algorithms to introduce a new class of fast algorithms for convolutional neural networks using C and OpenBLAS. It is a first implementation of the [Fast Algorithms for Convolutional Neural Networks - Paper](https://arxiv.org/pdf/1509.09308.pdf)
## Theoretical Background
The algorithm that we are going to implement relize the operation in the image below.
![Convolutional NN layer - kernel and tile](https://i.stack.imgur.com/9Iu89.gif)

We have extended the algorithm by using more dimension: <br/>
***Image***
* W x H x Channel -> It can be seen as a cube with a width, an height and a number of channels
* In our implementation we consider W = H
<br/>***Kernel***
* R x R x Channel x K -> It can be seen as a "list of cube". We have K cubes with a dimension of R x R x Channel
<br/>***Tile***
* A tile is a "Portion of the image" with a dimension D X D x Channel.
* The only requirement that we have is that D should be divisor of W (and H)
* We can consider (W / D) * (W / D) different tiles in a single image.
<br/>***Output***
* The output of a tile is a matrix of a dimension of M x M
* The program will give us as an output a cube with a dimension of M * (W / D) x M * (W / D) x K
* By using a tensorflow notation a single output tile is calculated as:
<br/>output[k,h,w] = sum_{c,dh,dw} input[c, h + dh, w + dw] * filter[k,c,h,k]
<br/>***but instead of using this formula we have implemented it with the winograd algorithm.***
### Parameters
* M - dimension of an output tile
* R - dimension of the kernel
* Channel - number of channel of the image and of the kernel
* K - number of kernels
* W - width of the image
* H - height of the image
* kernel.txt - file of the input kernel
* input.txt - file of the input image
## Getting Started
These instructions allow you to run the program on your computer.
### Prerequisites
* Install openBlas on your computer
* Install python with np on your computer
### Installing
We need to compile the program, you should give to the compiler the path of the openBLAS library.
I use this line: 
```
cc -static -o test main.c -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran
```
### Running an example
To run an example of the program you should follow these steps:
* Create the matrix parameters
```
python calcMc.py M R
```
* Create the input image
```
python createInput.py W H CHANNEL input.txt
```
* Create the input kernel
```
python createKernel.py R CHANNEL K kernel.txt
```
* Call the C executable program with the correct parameters
```
./test M R CHANNEL K W H kernel.txt input.txt
```
# Contributing
This project has been developed by [Me](https://github.com/istoony) and [Andrea Facchini](https://github.com/AndreF010203)
