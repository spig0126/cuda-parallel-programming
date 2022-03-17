/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  Ch03 image blur program
 *
 *        Version:  1.0
 *        Created:  02/08/2021 10:41:21 PM
 *       Revision:  none
 *       Compiler:  nvcc
 *
 *         Author:  Yeonhee Jung
 *
 * =====================================================================================
 */

#include <iostream>
#include "mkPpm.h"
#include "mkCuda.h"
#include "mkClockMeasure.h"

using namespace std;

const int MAX_ITER = 10;

void cpuCode(unsigned char *outArray, const unsigned char *inArray, const int w, const int h, int blurSize){
	float rVal, gVal, bVal;
	int curIndex, curRow, curCol, indexOffset;

	int totalSize = w * h * 3;
	int wSize = w * 3;
	int rblurSize = blurSize * wSize;
	int cblurSize = blurSize * 3;

	for (int r = 0; r < totalSize; r += wSize)
	{
		for (int c = 0; c < wSize; c += 3){
			rVal = gVal = bVal = 0;
			int pixels = 0;

			for (int rOffset = -rblurSize; rOffset <= rblurSize; rOffset+=wSize){
				for (int cOffset = -cblurSize; cOffset <= cblurSize; cOffset+=3){
					curRow = r + rOffset;
					curCol = c + cOffset;

					if (-1 < curRow && curRow < totalSize && -1 < curCol && curCol < wSize){
						indexOffset = curRow + curCol;
						rVal += inArray[indexOffset];
						gVal += inArray[indexOffset + 1];
						bVal += inArray[indexOffset + 2];
						pixels++;
					}
				}
			}

			curIndex = r + c;
			outArray[curIndex] = rVal / pixels;
			outArray[curIndex + 1] = gVal / pixels;
			outArray[curIndex + 2] = bVal / pixels;
		}
	}
}

__global__
void gpuCode(unsigned char *outArray, const unsigned char *inArray, const int w, const int h, const int blurSize){
	int c = threadIdx.x + blockIdx.x * blockDim.x;
	int r = threadIdx.y + blockIdx.y * blockDim.y;

	if(c < w && r < h){
		int curRow, curCol, curIndex, indexOffset, pixels = 0;
		float rVal = 0, gVal = 0, bVal = 0;

		for (int rOffset = -blurSize; rOffset <= blurSize; rOffset++){
			for (int cOffset = -blurSize; cOffset <= blurSize; cOffset++){
				curRow = r + rOffset;
				curCol = c + cOffset;
				if(-1 < curRow && curRow < h && -1 < curCol && curCol < w){
					indexOffset = (curRow * w + curCol) * 3;
					rVal += inArray[indexOffset];
					gVal += inArray[indexOffset + 1];
					bVal += inArray[indexOffset + 2];
					pixels++;	
				}
			}
		}

		curIndex = (r * w + c) * 3;
		outArray[curIndex] = rVal / pixels;
		outArray[curIndex + 1] = gVal / pixels;
		outArray[curIndex + 2] = bVal / pixels;
	}
}

int main(){
	int w, h, blurSize = 10;
	unsigned char *h_imageArray;
	unsigned char *h_outImageArray;
	unsigned char *d_imageArray;
	unsigned char *d_outImageArray;
	unsigned char *h_outImageArray2;

	ppmLoad("./data/ewha_picture.ppm", &h_imageArray, &w, &h);

	size_t arraySize = sizeof(unsigned char) * h * w * 3;

	h_outImageArray = (unsigned char*)malloc(arraySize);
	h_outImageArray2 = (unsigned char*)malloc(arraySize);

	cudaError_t err = cudaMalloc((void **) &d_imageArray, arraySize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_outImageArray, arraySize);
	checkCudaError(err);

	err = cudaMemcpy(d_imageArray, h_imageArray, arraySize, cudaMemcpyHostToDevice);
	checkCudaError(err);

	const int tSize = 16;
	dim3 blockSize(tSize, tSize, 1);
	dim3 gridSize(ceil((float)w/tSize), ceil((float)h/tSize), 1);

	mkClockMeasure *ckCpu = new mkClockMeasure("CPU CODE");
	mkClockMeasure *ckGpu = new mkClockMeasure("GPU CODE");

	ckCpu->clockReset();
	ckGpu->clockReset();


	for(int i = 0; i < MAX_ITER; i++){
		
		ckCpu->clockResume();
		cpuCode(h_outImageArray, h_imageArray, w, h, blurSize);
		ckCpu->clockPause();

		ckGpu->clockResume();
		gpuCode<<<gridSize, blockSize>>>(d_outImageArray, d_imageArray, w, h, blurSize);
		err=cudaDeviceSynchronize();
		ckGpu->clockPause();
		checkCudaError(err);

	}
	ckCpu->clockPrint();
	ckGpu->clockPrint();

	err = cudaMemcpy(h_outImageArray2, d_outImageArray, arraySize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
			
	ppmSave("ewha_picture_cpu.ppm", h_outImageArray, w, h);
	ppmSave("ewha_picture_gpu.ppm", h_outImageArray2, w, h);
	return 0;
}


