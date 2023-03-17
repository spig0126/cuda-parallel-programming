CC=nvcc   

# num.out: bfp.cuh bfpStruct.cuh mkCuda.h error_analysis_program.cu 
# 	$(CC) -o errors.out error_analysis_program.cu

main.out: bfp.cuh bfpStruct.cuh mkCuda.h main.cu 
	$(CC) -o div.out main.cu

clean:
	rm -rf mult.out
