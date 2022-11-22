CC=nvcc 

a.out:  mkClockMeasure.h  mkCuda.h FFT.cu
	# $(CC) -o fft.out FFT.cu --ptxas-options=-v --use_fast_math -lcufft

clean:
	rm -rf temp.out
