CC=nvcc 

# test.out: mkClockMeasure.h  mkCuda.h test.cu
# 	$(CC) -o test.out test.cu

# clean:
# 	rm -rf ./test.out

# 	# rm -rf RTTNW.out *.ppm

main.out: mkClockMeasure.h mkCuda.h bfp.cuh bfpStruct.cuh main.cu
	$(CC) -o main.out main.cu

clean:
	rm -rf ./main.out