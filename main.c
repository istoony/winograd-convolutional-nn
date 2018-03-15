/*
F(4,3)
cc -static -o test_cblas_open main43.c -I /opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas -lpthread -lgfortran
	->For each cube of kernel
		-> for each tile (a tile is a cube d x d x channel)
			-> for each channel
				-> Apply the winograd algorithm
			-> sum the result of the previous channel to the result of the next channell
		-> write the tile in a single output layer
	->Add the output layer in che cube layer and change the kernel cube
*/

#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#define MAX_ELEMENTS 100
#define FILENAMELEN 10
#define DDIMENSION(m,r) (m+r-1)
#define OUTPUTDIMENSION(m,w,d) (m * (w / d))

typedef float * Matrix;
typedef Matrix * Cube;
typedef Cube * Hypercube;

void productABAt(float result[], float A[], float B[], int rowsA, int colsB, int colsA);
void twoDcorrelation(float C[], float g[], float D[], int m, int r, int d);
Matrix elaborateTile(Cube tile, Cube g, int m, int r, int d, int channel);
Cube fillTheTile(Cube image, int starth, int startw, int w, int d, int channel);
void saveOutputTile(Matrix output, Matrix tileOutput, int starth, int startw, int w, int d, int m);
Matrix elaborateKernel(Cube image, Cube g, int m, int r, int channel, int w, int h, int d);

/*START FUNCTION TO CREATE VARIABLES*/
	Cube generateCube(int w, int h, int ch);
	Matrix generateMatrix(int w, int h);
	Hypercube readKernels(int w, int h, int ch, int k, char fileName[]);
	Cube readInput(int w, int h, int ch, char fileName[]);

	void printCube(Cube c, int w, int channel);
/*END FUNCTION TO CREARE VARIABLES*/

/*START FUNCTIONS TO READ THE FILES*/
	void readFile(float result[], char fileName[]);
	void generateNameAT(char* at, int m, int r);
	void generateNameG(char* g, int m, int r);
	void generateNameBT(char* g, int m, int r);
/* END FUNCRIONS TO READ THE FILES*/


static Matrix A;	//Parameters that are used to make the lowest operation.
static Matrix B;	//
static Matrix G;	//

void main(int argc, char *argv[])
{
	
	int m = atoi(argv[1]);			//dimension of the output tile
	int r = atoi(argv[2]);			//dimension of the kernel rxr
	int channel = atoi(argv[3]);	//dimension of the channel
	int k = atoi(argv[4]);			//number of kernel
	int w = atoi(argv[5]);			//width of the image -> multiple of d
	int h = atoi(argv[6]);			//height of the image -> multiple of d

	char * inputFilename = argv[7];
	char * kernelFilename = argv[8];
		
	int d = DDIMENSION(m,r);

	char fnameAT[FILENAMELEN];
	char fnameG[FILENAMELEN];
	char fnameBT[FILENAMELEN];

	generateNameAT(fnameAT, m,r);
	generateNameG(fnameG, m,r);
	generateNameBT(fnameBT, m,r);
	
	A = generateMatrix(m, d);
	B = generateMatrix(d, d);
	G = generateMatrix(r, d);

	readFile(A, fnameAT);
	readFile(B, fnameBT);
	readFile(G, fnameG);
	
	Cube image = readInput(w, h, channel, inputFilename);
	Hypercube kernel = readKernels(r, r, channel, k, kernelFilename);

	printf("IMAGE\n");
	printCube(image, w, channel);

	printf("\nKERNELS\n");
	for(int i = 0; i < k; i++)
	{
		printf("Kernel number %d\n", i);
		printCube(kernel[i], r, channel);
	}
	printf("\n------------------\n");
	Cube output = (Cube)malloc(k * sizeof(Matrix));	//k is the number of filters
	
	//printf("START THE COMPUTATION - SLEEP FOR 3 SECONDS\n");
 	//sleep(3);
	/*
	 *	START COMPUTATION.
	 */
	clock_t start = clock();
	for(int i = 0; i < k; i ++)
		output[i] = elaborateKernel(image, kernel[i], m, r, channel, w, h, d);
	
	clock_t end = clock();
	/*
	 *	END COMPUTATION.
	 */
	
	printf("OUTPUT:\n");
	printCube(output, OUTPUTDIMENSION(m,w,d), k);
	
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	
	printf("\nTIME IN SECONDS TO DO F(%d,%d) = %f\n",m,r,seconds);
	
}

void printCube(Cube c, int w, int channel)
{
	for(int i = 0; i < channel; i++)
	{
		printf("channel %d:\n",i);
		for(int j = 0; j < w * w; j ++)
		{
			printf("%.2f ", c[i][j]);
			if((j+1) % w == 0) printf("\n");
		}
	}
}


Matrix elaborateKernel(Cube image, Cube g, int m, int r, int channel, int w, int h, int d)
{
	Matrix output = generateMatrix(OUTPUTDIMENSION(m,w,d), OUTPUTDIMENSION(m,w,d));
	
	for(int i = 0; i < h; i = i + d)
	{
		for (int j = 0; j < w; j = j + d)
		{
			Cube tile = fillTheTile(image, i, j, w, d, channel);
			
			Matrix tileOutput = elaborateTile(tile, g, m, r, d, channel);
			saveOutputTile(output, tileOutput, i, j, w, d, m);

			free(tileOutput);
			free(tile);
		}
	}
	return output;
}


/*
	Write the output tile on the output image.
*/

void saveOutputTile(Matrix output, Matrix tileOutput, int starth, int startw, int w, int d, int m)
{
		int outputWidth = OUTPUTDIMENSION(m,w,d);
		int offsetw = (startw / d) * m;
		int offseth = (starth / d) * m;
		for(int y = 0; y < m; y ++)
			for(int x = 0; x < m; x++)
				output[(x + offsetw) + (outputWidth * (y + offseth))] = tileOutput[x + (m * y)];

}

/*
	This function takes an image, it divides the image based on starth and startw and it returns a tile.
		- image  -> the all image
		- starth -> height coordinate of the image where to start to make the tile
		- startw -> width coordinate of the image where to start to make the tile
		- w -> total width of the image
		- d -> dimension of the tile -> d x d
		- channel -> number of channels of the image
	RESULT
		- A cube d x d x channell that represent a single tile.
*/
Cube fillTheTile(Cube image, int starth, int startw, int w, int d, int channel)
{
	Cube tile = generateCube(d, d, channel);
	for(int c = 0; c < channel; c++)
		for(int y = 0; y < d; y ++)
			for(int x = 0; x < d; x++)
				tile[c][x + (d * y)] = image[c][(startw + x) + (w * (y + starth))];
	return tile;
}

/*
	This function is used to make the combination with a tile and a kernel.
		- tile has the dimension d x d and has 'channel' channels
		- g is the kernel, it has the dimension of r x r and has 'channel' channels
		- m -> dimension of the output
		- r -> dimension of the kernel
		- d -> dimension of the tile
		- channell -> number of channels
	RESULT
		- it returns a matrix calculated by an element wise sum of the result of the function twoDcorrelation 
			calculated for each channel.
*/
Matrix elaborateTile(Cube tile, Cube g, int m, int r, int d, int channel)
{
	Matrix temp = generateMatrix(d,d);
	Matrix tileOutputBig = generateMatrix(d,d);
	Matrix tileOutput = generateMatrix(m,m);

	for(int i = 0; i < m * m; i++)
			tileOutput[i] = 0;

	for(int i = 0; i < channel; i++)
	{
	 	twoDcorrelation(temp, g[i], tile[i], m, r, d);
		for(int j = 0 ; j < d * d; j++)
			tileOutputBig[j] = tileOutputBig[j] + temp[j];
		
	}
	productABAt(tileOutput, A, tileOutputBig, m, d, d);
	free(temp);
	free(tileOutputBig);

	return tileOutput;
}

Cube generateCube(int w, int h, int ch)
{
	Cube cube = (Cube)malloc(ch * sizeof(Matrix));
	for(int i = 0; i < ch; i++)
		cube[i] = generateMatrix(w,h);

	return cube;
}

Matrix generateMatrix(int w, int h)
{
	return (Matrix)malloc(w * h * sizeof(float));
}

/*
	Calculate a single channel for a single kernel for a single tile
	Input Parameter:
		- result -> matrix of the result
		- g -> the considered kernel
		- D -> the considered tile
		- m -> Dimension of the output
		- r -> dimension of the kernel
		- d -> dimension of the tile (m + r - 1)
*/
void twoDcorrelation(float C[], float g[], float D[], int m, int r, int d)
{
	
	
	float GgG[MAX_ELEMENTS];
	float BDB[MAX_ELEMENTS];
	int i=0;

	productABAt(GgG, G, g, d, r, r);
	productABAt(BDB, B, D, d, d, d);

 	for(i=0;i<d*d;i++)
		C[i] = GgG[i] * BDB[i];
  	//productABAt(C, A, GgG, m, d, d);
}


/*
 *   Make the operation R = ABA**T
 *   Parameters:
 *		- Matrix Result
 *		- Matrix A and B
 *		- Number of rows of A -> rowsA
 *		- Number of cols of B -> colsB 
 *		- Number of cols of A that is equal to the numer of rows of B -> colsA
 *  Result
 *		- Matrix result -> that is an rowsA x rowsA matrix
 */
void productABAt(float result[], float A[], float B[], int rowsA, int colsB, int colsA)
{
	float C[MAX_ELEMENTS];

	cblas_sgemm(CblasRowMajor,	//Modo in cui è salvata la matrice, legge i numeri riga per riga 
		CblasNoTrans,	//Non fare trasposta 
		CblasNoTrans, 	//Non fare trasposta
		rowsA, 		//numero righe matrice A
		colsB, 		//numero colonne matrice B
		colsA, 		//numero colonne A & numero righe B
		1, 		//moltiplicatore prima matrice
		A, 		//Matrice A
		colsA, 		//numero colonne matrice A
		B, 		//Matrice B
		colsB, 		//numero colonne matrice B
		0, 		//Moltiplicatore matrice C
		C,		//matrice C
		colsB);		//numero righe matrice C

	//C has rowsA rows and colsB cols
	//the result has rowsA rows and rowsA cols 

	cblas_sgemm(CblasRowMajor,
		CblasNoTrans,
		CblasTrans,
		rowsA,	//numero righe matrice C
		rowsA,
		colsB,
		1,
		C,
		colsB,
		A,
		colsA,
		0,
		result,
		rowsA);
}

Cube readInput(int w, int h, int ch, char fileName[])
{

    FILE *f;
    f = fopen(fileName, "r");
    if(f==NULL)
	{
    	printf("I can not read the file %s\n", fileName);
		exit(0);
	}
			
    int fi = 0;
    int linesize = w*h;
    Cube matrix = (Cube)malloc(ch * sizeof(Matrix));

    for(int i = 0; i < ch; i++)
	{
        matrix[i] = (Matrix)malloc(linesize * sizeof(float));
        while(!feof(f)&&fi<linesize)
		{
            fscanf(f, "%f", &matrix[i][fi]);
            fi++;
        }
        fi=0;
    }
    fclose(f);
    return matrix;

}

Hypercube readKernels(int w, int h, int ch, int k, char fileName[])
{

    FILE *f;
    f = fopen(fileName, "r");
    if(f==NULL)
    {
		printf("I can not read the file %s\n", fileName);
		exit(0);
	}
    int fi = 0;
    int linesize = w*h;
    int kernelsize = linesize*ch;
    Hypercube matrix = (Hypercube)malloc(k * sizeof(Cube));

    for(int i = 0; i < k; i++)
	{
        matrix[i] = (Cube)malloc(kernelsize * sizeof(Matrix));
        for(int j=0; j< ch; j++)
		{
            matrix[i][j] = (Matrix)malloc(linesize * sizeof(float));
            while(!feof(f) && fi<linesize)
			{
                   fscanf(f, "%f", &matrix[i][j][fi]);
                    fi++;
            }
            fi=0;
        }
    }
    fclose(f);
    return matrix;
}



void readFile(float result[], char fileName[])
{
	FILE *f; 
	f = fopen(fileName, "r");
	if(f == NULL)
	{
		printf("I can not read the file %s\n", fileName);
		exit(0);
	}
	int i = 0;
	while(!feof(f))
	{
		fscanf(f, "%f", &result[i]);	
		i++;
	}
	fclose(f);
}

/* FILENAMES FUNCTIONS */
void generateNameAT(char* at, int m, int r){
    at[0] = 'A';
    at[1] = 'T';
    at[2] = m + '0';
    at[3] = 'x';
    at[4] = r + '0';
    at[5] = '.';
    at[6] = 't';
    at[7] = 'x';
    at[8] = 't';
    at[9] = '\0';
}

void generateNameG(char* g, int m, int r){
    g[0] = 'G';
    g[1] = m + '0';
    g[2] = 'x';
    g[3] = r + '0';
    g[4] = '.';
    g[5] = 't';
    g[6] = 'x';
    g[7] = 't';
    g[8] = '\0';
}

void generateNameBT(char* g, int m, int r){
    g[0] = 'B';
    g[1] = 'T';
    g[2] = m + '0';
    g[3] = 'x';
    g[4] = r + '0';
    g[5] = '.';
    g[6] = 't';
    g[7] = 'x';
    g[8] = 't';
    g[9] = '\0';
}
