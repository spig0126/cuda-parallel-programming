CC=nvcc 

RTTNW.out:  mkClockMeasure.h  mkCuda.h main.cu 
	$(CC) -o ../makefiles/RTTNW_400_3_5.out main.cu

# test.out: mkClockMeasure.h  mkCuda.h test.cu
# 	$(CC) -o test.out test.cu

clean:
	rm -rf test.out

	# rm -rf RTTNW.out *.ppm