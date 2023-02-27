CC=nvcc 

main.out: mkClockMeasure.h mkCuda.h main.cu
	$(CC) -o main.out main.cu

clean:
	rm -rf ./main.out