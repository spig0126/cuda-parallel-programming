CC=nvcc 

a.out:  mkClockMeasure.h  mkCuda.h FFT.cu
	$(CC) -o fft.out FFT.cu

clean:
	rm -rf fft.out
