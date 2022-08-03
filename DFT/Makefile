CC=nvcc 

a.out:  mkClockMeasure.h  mkCuda.h DFT.cu
	$(CC) -o dft.out DFT.cu

clean:
	rm -rf dft.out
