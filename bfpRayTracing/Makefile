CC=nvcc 

RTTNW.out:  mkClockMeasure.h  mkCuda.h main.cu 
	$(CC) -o RTTNW.out main.cu

# test.out: mkClockMeasure.h  mkCuda.h test.cu
# 	$(CC) -o test.out test.cu

clean:
	rm -rf test.out

	# rm -rf RTTNW.out *.ppm