#include <stdio.h>
#include <string.h>
#include "timer.h"

#define image_height 102400
#define image_width 255

#define filter_height 5
#define filter_width 5

#define border_height ((filter_height/2)*2)
#define border_width ((filter_width/2)*2)
#define input_height (image_height + border_height)
#define input_width (image_width + border_width)

#define block_size_x 32
#define block_size_y 16

#define SEED 1234

using namespace std;

void convolutionSeq(float *output, float *input, float *filter) {
	//for each pixel in the output image

	timer sequentialTime = timer("Sequential");

	sequentialTime.start();

	for (int y=0; y < image_height; y++) {
		for (int x=0; x < image_width; x++) {

			//for each filter weight
			for (int i=0; i < filter_height; i++) {
				for (int j=0; j < filter_width; j++) {
					output[y*image_width+x] += input[(y+i)*input_width+x+j] * filter[i*filter_width+j];
				}
			}

		}
	}
	sequentialTime.stop();
	cout << "convolution (sequential): \t\t" << sequentialTime << endl;

}


__global__ void convolution_kernel_naive(float *output, float *input, float *filter) {
	unsigned x = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned y = blockDim.y * blockIdx.y + threadIdx.y;
	float sum = 0;
	unsigned int image_height_block = image_height / ceilf(image_height  / (1024.) /* optimal size */);

	if(y >= image_height_block || x > image_width)
		return;
	//for each filter weight
	for (int i=0; i < filter_height; i++) {
		for (int j=0; j < filter_width; j++) {
			sum += input[(y+i)*input_width+x+j] * filter[i*filter_width+j];
		}
	}

	output[y*image_width+x] = sum;

}

void convolutionCUDA(float *output, float *input, float *filter) {
	float *d_input, *d_input2; float *d_output; float *d_filter;
	cudaStream_t streams[2];
	cudaStream_t *stream1 = streams, *stream2 = streams + 1;
	cudaError_t err;
	unsigned int image_height_block = image_height / ceilf(image_height  / (1024.) /* optimal size */);
	unsigned int input_height_block = (image_height_block + border_height);
	timer kernelTime = timer("kernelTime");
	timer memoryTime = timer("memoryTime");
	cudaStreamCreate(stream1);
	cudaStreamCreate(stream2);
	// memory allocation
	err = cudaMalloc((void **)&d_input, input_height_block*input_width*sizeof(float));
	if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_input: %s\n", cudaGetErrorString( err )); }
	err = cudaMalloc((void **)&d_input2, input_height_block*input_width*sizeof(float));
	if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_input: %s\n", cudaGetErrorString( err )); }
	err = cudaMalloc((void **)&d_output, image_height_block*image_width*sizeof(float));
	if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_output: %s\n", cudaGetErrorString( err )); }
	err = cudaMalloc((void **)&d_filter, filter_height*filter_width*sizeof(float));
	if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMalloc d_filter: %s\n", cudaGetErrorString( err )); }

	memoryTime.start();
	err = cudaMemcpyAsync(d_input, input, input_height_block*input_width*sizeof(float), cudaMemcpyHostToDevice, *stream2);
	if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy host to device input: %s\n", cudaGetErrorString( err ));  }
	err = cudaMemcpyAsync(d_filter, filter, filter_height*filter_width*sizeof(float), cudaMemcpyHostToDevice, *stream2);
	if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy host to device filter: %s\n", cudaGetErrorString( err ));  }

	err = cudaMemsetAsync(d_output, 0, image_height_block*image_width*sizeof(float), *stream2);
	if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemset output: %s\n", cudaGetErrorString( err ));  }
	memoryTime.stop();

	{float* tmp = d_input; d_input = d_input2; d_input2 = tmp;}


	kernelTime.start();
	for(int i = 0; i < image_height; ) {
		{float* tmp = d_input; d_input = d_input2; d_input2 = tmp;}
		dim3 threads(block_size_x, block_size_y);
		dim3 grid(int(ceilf(image_width/(float)threads.x)), int(ceilf(image_height_block/(float)threads.y)) );
		cudaStreamSynchronize(*stream2);
		convolution_kernel_naive<<<grid, threads, 0, *stream1>>>(d_output, d_input, d_filter);
		memoryTime.start();
		i += image_height_block;
		if(i >= image_height) {
			kernelTime.stop();
			memoryTime.stop();
			goto cont;
		}
		err = cudaMemcpyAsync(d_input2, input + (i/image_height_block) * image_height_block * input_width, input_height_block*input_width*sizeof(float), cudaMemcpyHostToDevice, *stream2);
		if (err != cudaSuccess) { fprintf(stderr, "Error in cudaMemcpy host to device input: %s\n", cudaGetErrorString( err ));}

		memoryTime.stop();

cont:
		cudaStreamSynchronize(*stream1);
		err = cudaMemcpyAsync(output + (i - image_height_block) * image_width, d_output, image_height_block*image_width*sizeof(float), cudaMemcpyDeviceToHost, *stream2);
	}

	cudaDeviceSynchronize();
	err = cudaFree(d_input);
	if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_input: %s\n", cudaGetErrorString( err )); }
	err = cudaFree(d_output);
	if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_output: %s\n", cudaGetErrorString( err )); }
	err = cudaFree(d_filter);
	if (err != cudaSuccess) { fprintf(stderr, "Error in freeing d_filter: %s\n", cudaGetErrorString( err )); }

	cout << "convolution (kernel): \t\t" << kernelTime << endl;
	cout << "convolution (memory): \t\t" << memoryTime << endl;

}

int compare_arrays(float *a1, float *a2, int n) {
	int errors = 0;
	int print = 0;

	for (int i=0; i<n; i++) {

		if (isnan(a1[i]) || isnan(a2[i])) {
			errors++;
			if (print < 10) {
				print++;
				fprintf(stderr, "Error NaN detected at i=%d,\t a1= %10.7e \t a2= \t %10.7e\n",i,a1[i],a2[i]);
			}
		}

		float diff = (a1[i]-a2[i])/a1[i];
		if (diff > 1e-6f) {
			errors++;
			if (print < 10) {
				print++;
				fprintf(stderr, "Error detected at i=%d, \t a1= \t %10.7e \t a2= \t %10.7e \t rel_error=\t %10.7e\n",i,a1[i],a2[i],diff);
			}
		}

	}

	return errors;
}


int main() {
	int i;
	int errors=0;

	//allocate arrays and fill them
	float *input = (float *) malloc(input_height * input_width * sizeof(float));
	float *output1 = (float *) calloc(image_height * image_width, sizeof(float));
	float *output2 = (float *) calloc(image_height * image_width, sizeof(float));
	float *filter = (float *) malloc(filter_height * filter_width * sizeof(float));

	for (i=0; i< input_height * input_width; i++) {
		input[i] = (float) (i % SEED);
	}

//THis is specific for a W==H smoothening filteri, where W and H are odd.
	for (i=0; i<filter_height * filter_width; i++) {
		filter[i] = 1.0;
	}

	for (i=filter_width+1; i<(filter_height - 1) * filter_width; i++) {
		if (i % filter_width > 0 && i % filter_width < filter_width-1) filter[i]+=1.0;
	}

	filter[filter_width*filter_height/2]=3.0;
//end initialization

	//measure the CPU function
	convolutionSeq(output1, input, filter);
	//measure the GPU function
	convolutionCUDA(output2, input, filter);


	//check the result
	errors += compare_arrays(output1, output2, image_height*image_width);
	if (errors > 0) {
		printf("TEST FAILED! %d errors!\n", errors);
	} else {
		printf("TEST PASSED!\n");
	}


	free(filter);
	free(input);
	free(output1);
	free(output2);

	return 0;
}
