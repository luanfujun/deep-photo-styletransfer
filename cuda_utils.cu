extern "C" {
	#include "lua.h"
	#include "lualib.h"
	#include "lauxlib.h"
}

#include "luaT.h"
#include "THC.h"

#include <stdio.h>
#include <assert.h>
#include <math_constants.h>
#include <math_functions.h>
#include <stdint.h>
#include <unistd.h>

#define TB 256
#define EPS 1e-4

THCState* getCutorchState(lua_State* L)
{
	lua_getglobal(L, "cutorch");
	lua_getfield(L, -1, "getState");
	lua_call(L, 0, 1);
	THCState *state = (THCState*) lua_touserdata(L, -1);
	lua_pop(L, 2);
	return state;
}

void checkCudaError(lua_State *L) {
	cudaError_t status = cudaPeekAtLastError();
	if (status != cudaSuccess) {
		luaL_error(L, cudaGetErrorString(status));
	}
}

THCudaTensor *new_tensor_like(THCState *state, THCudaTensor *x)
{
	THCudaTensor *y = THCudaTensor_new(state);
	THCudaTensor_resizeAs(state, y, x);
	return y;
}

__global__ void matting_laplacian_kernel(
	float *input, float *grad, int h, int w, 
	int *CSR_rowIdx, int *CSR_colIdx, float *CSR_val,
	int N
)
{	
	int size = h * w;
	int _id = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (_id < size) {
		int x = _id % w, y = _id / w;
		int id = x * h + y;

		/// Because matting laplacian L is systematic, sum row is sufficient
		// 1.1 Binary search
		int start = 0;
		int end = N-1;
		int mid = (start + end)/2;
		int index = -1;
		while (start <= end) {
			int rowIdx = (CSR_rowIdx[mid]) - 1;
			
			if (rowIdx == id) {
				index = mid;    break;
			}
			if (rowIdx > id) {
				end = mid - 1;
				mid = (start + end)/2;
			} else {
				start = mid + 1;
				mid = (start + end)/2;
			}
		}
		if (index != -1) {
			// 1.2 Complete range
			int index_s = index, index_e = index;
			while ( index_s >= 0 && ((CSR_rowIdx[index_s] - 1) == id) )
				index_s--;
			while ( index_e <  N && ((CSR_rowIdx[index_e] - 1) == id) )
				index_e++;
			// 1.3 Sum this row
			for (int i = index_s + 1; i < index_e; i++) {
				//int rowIdx = CSR_rowIdx[i] - 1; 
				int _colIdx = (CSR_colIdx[i]) - 1;
				float val  = CSR_val[i];

				int _x = _colIdx / h, _y = _colIdx % h;
				int colIdx = _y *w + _x;

				grad[_id] 			+= 2*val * input[colIdx];
				grad[_id + size] 	+= 2*val * input[colIdx + size];
				grad[_id + 2*size]  += 2*val * input[colIdx + 2*size];
			}

		}	
	}

	return ;
}
 
//cuda_utils.matting_laplacian(input, h, w, CSR_rowIdx, CSR_colIdx, CSR_val, CSC_rowIdx, CSC_colIdx, CSC_val, N)
int matting_laplacian(lua_State *L) {
	THCState *state = getCutorchState(L);
	THCudaTensor *input 			= (THCudaTensor*)luaT_checkudata(L, 1, "torch.CudaTensor");
	int h 							= luaL_checknumber(L, 2);
	int w 							= luaL_checknumber(L, 3);
	THCudaIntTensor *CSR_rowIdx		= (THCudaIntTensor*)luaT_checkudata(L, 4, "torch.CudaIntTensor");
	THCudaIntTensor *CSR_colIdx		= (THCudaIntTensor*)luaT_checkudata(L, 5, "torch.CudaIntTensor");
	THCudaTensor *CSR_val			= (THCudaTensor*)luaT_checkudata(L, 6, "torch.CudaTensor");
	int N 							= luaL_checknumber(L, 7);
	THCudaTensor *grad = new_tensor_like(state, input);
	THCudaTensor_zero(state, grad);

	matting_laplacian_kernel<<<(h*w-1)/TB+1, TB>>>(
		THCudaTensor_data(state, input),
		THCudaTensor_data(state, grad),
		h, w,
		THCudaIntTensor_data(state, CSR_rowIdx),
		THCudaIntTensor_data(state, CSR_colIdx),
		THCudaTensor_data(state, CSR_val),
		N
	);
	checkCudaError(L);
	luaT_pushudata(L, grad, "torch.CudaTensor");
	return 1;
}

__device__ bool InverseMat4x4(double m_in[4][4], double inv_out[4][4]) {
	double m[16], inv[16];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			m[i * 4 + j] = m_in[i][j];
		}
	}

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    double det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (abs(det) < 1e-9) {
        return false;
    }


    det = 1.0 / det; 

    for (int i = 0; i < 4; i++) {
    	for (int j = 0; j < 4; j++) {
    		inv_out[i][j] = inv[i * 4 + j] * det;
    	}
    }

    return true;
}

__global__ void best_local_affine_kernel(
	float *output, float *input, float *affine_model,
	int h, int w, float epsilon, int kernel_radius
)
{
	int size = h * w;
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (id < size) {
		int x = id % w, y = id / w;
		
		double Mt_M[4][4] = {}; // 4x4
		double invMt_M[4][4] = {};
		double Mt_S[3][4] = {}; // RGB -> 1x4
		double A[3][4] = {}; 
		for (int i = 0; i < 4; i++) 
			for (int j = 0; j < 4; j++) {
				Mt_M[i][j] = 0, invMt_M[i][j] = 0;
				if (i != 3) {
					Mt_S[i][j] = 0, A[i][j] = 0;
					if (i == j)
			    		Mt_M[i][j] = 1e-3;
			    }
			}

		for (int dy = -kernel_radius; dy <= kernel_radius; dy++) {
			for (int dx = -kernel_radius; dx <= kernel_radius; dx++) {
				
				int xx = x + dx, yy = y + dy;
				int id2 = yy * w + xx;

				if (0 <= xx && xx < w && 0 <= yy && yy < h) {

					Mt_M[0][0] += input[id2 + 2*size] * input[id2 + 2*size];
					Mt_M[0][1] += input[id2 + 2*size] * input[id2 + size];
					Mt_M[0][2] += input[id2 + 2*size] * input[id2];
					Mt_M[0][3] += input[id2 + 2*size];

					Mt_M[1][0] += input[id2 + size] * input[id2 + 2*size];
					Mt_M[1][1] += input[id2 + size] * input[id2 + size];
					Mt_M[1][2] += input[id2 + size] * input[id2];
					Mt_M[1][3] += input[id2 + size];

					Mt_M[2][0] += input[id2] * input[id2 + 2*size];
					Mt_M[2][1] += input[id2] * input[id2 + size];
					Mt_M[2][2] += input[id2] * input[id2];
					Mt_M[2][3] += input[id2];
					 
					Mt_M[3][0] += input[id2 + 2*size];
					Mt_M[3][1] += input[id2 + size];
					Mt_M[3][2] += input[id2];
					Mt_M[3][3] += 1;

					Mt_S[0][0] += input[id2 + 2*size] * output[id2 + 2*size];
					Mt_S[0][1] += input[id2 + size] * output[id2 + 2*size];
					Mt_S[0][2] += input[id2] * output[id2 + 2*size];
					Mt_S[0][3] += output[id2 + 2*size];

					Mt_S[1][0] += input[id2 + 2*size] * output[id2 + size];
					Mt_S[1][1] += input[id2 + size] * output[id2 + size];
					Mt_S[1][2] += input[id2] * output[id2 + size];
					Mt_S[1][3] += output[id2 + size];

					Mt_S[2][0] += input[id2 + 2*size] * output[id2];
					Mt_S[2][1] += input[id2 + size] * output[id2];
					Mt_S[2][2] += input[id2] * output[id2];
					Mt_S[2][3] += output[id2];
				}
			}
		}

		bool success = InverseMat4x4(Mt_M, invMt_M);

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++) {
					A[i][j] += invMt_M[j][k] * Mt_S[i][k];
				}
			}
		}

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
				int affine_id = i * 4 + j;
				affine_model[12 * id + affine_id] = A[i][j];
			}
		}



	}
	return ;
}


__global__ void bilateral_smooth_kernel(
	float *affine_model, float *filtered_affine_model, float *guide, 
	int h, int w, int kernel_radius, float sigma1, float sigma2
)	
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int size = h * w;
	if (id < size) {
		int x = id % w;
		int y = id / w;

		double sum_affine[12] = {};
		double sum_weight = 0;
		for (int dx = -kernel_radius; dx <= kernel_radius; dx++) {
			for (int dy = -kernel_radius; dy <= kernel_radius; dy++) {
				int yy = y + dy, xx = x + dx;
				int id2 = yy * w + xx;
				if (0 <= xx && xx < w && 0 <= yy && yy < h) {
					float color_diff1 = guide[yy*w + xx] - guide[y*w + x];
					float color_diff2 = guide[yy*w + xx + size] - guide[y*w + x + size];
					float color_diff3 = guide[yy*w + xx + 2*size] - guide[y*w + x + 2*size];
					float color_diff_sqr = 
						(color_diff1*color_diff1 + color_diff2*color_diff2 + color_diff3*color_diff3) / 3;

					float v1 = exp(-(dx * dx + dy * dy) / (2 * sigma1 * sigma1));
					float v2 = exp(-(color_diff_sqr) / (2 * sigma2 * sigma2));
					float weight = v1 * v2;

					for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 4; j++) {
							int affine_id = i * 4 + j;
							sum_affine[affine_id] += weight * affine_model[id2*12 + affine_id];
						}
					}
					sum_weight += weight;
				}
			}
		}

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
				int affine_id = i * 4 + j;
				filtered_affine_model[id*12 + affine_id] = sum_affine[affine_id] / sum_weight;
			}
		}
	}
	return ;
}


__global__ void reconstruction_best_kernel(
	float *input, float *filtered_affine_model, float *filtered_best_output,
	int h, int w
) 
{	
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int size = h * w;
	if (id < size) {
		double out1 = 
			input[id + 2*size] * filtered_affine_model[id*12 + 0] + // A[0][0] + 
			input[id + size]   * filtered_affine_model[id*12 + 1] + // A[0][1] + 
			input[id]          * filtered_affine_model[id*12 + 2] + // A[0][2] +
								 filtered_affine_model[id*12 + 3]; //A[0][3];
		double out2 = 
			input[id + 2*size] * filtered_affine_model[id*12 + 4] + //A[1][0] + 
			input[id + size]   * filtered_affine_model[id*12 + 5] + //A[1][1] + 
			input[id]          * filtered_affine_model[id*12 + 6] + //A[1][2] +
								 filtered_affine_model[id*12 + 7]; //A[1][3];
		double out3 = 
			input[id + 2*size] * filtered_affine_model[id*12 + 8] + //A[2][0] + 
			input[id + size]   * filtered_affine_model[id*12 + 9] + //A[2][1] + 
			input[id]          * filtered_affine_model[id*12 + 10] + //A[2][2] +
								 filtered_affine_model[id*12 + 11]; // A[2][3];

		filtered_best_output[id] = out1;
		filtered_best_output[id + size] = out2;
		filtered_best_output[id + 2*size] = out3;
	}
	return ;
}

// local best01 = cuda_utils.smooth_local_affine(output01, input01, epsilon, patch, h, w, filter_radius, sigma1, sigma2)
int smooth_local_affine(lua_State *L) {
	THCState *state = getCutorchState(L);
	THCudaTensor *output  = (THCudaTensor*)luaT_checkudata(L, 1, "torch.CudaTensor");
	THCudaTensor *input   = (THCudaTensor*)luaT_checkudata(L, 2, "torch.CudaTensor");
	float         epsilon = luaL_checknumber(L, 3);
	int 		  patch   = luaL_checknumber(L, 4);
	int 		  h       = luaL_checknumber(L, 5);
	int 		  w       = luaL_checknumber(L, 6);
	int 		  f_r     = luaL_checknumber(L, 7);
	float 		  sigma1  = luaL_checknumber(L, 8);
	float 		  sigma2  = luaL_checknumber(L, 9);

	THCudaTensor *filtered_best_output = new_tensor_like(state, input);
	THCudaTensor_zero(state, filtered_best_output);

	THCudaTensor *affine_model = THCudaTensor_new(state);
	THCudaTensor_resize2d(state, affine_model, h*w, 12);
	THCudaTensor_zero(state, affine_model);

	THCudaTensor *filtered_affine_model = THCudaTensor_new(state);
	THCudaTensor_resize2d(state, filtered_affine_model, h*w, 12);
	THCudaTensor_zero(state, filtered_affine_model);

	int radius = (patch-1) / 2;

	best_local_affine_kernel<<<(h*w)/TB+1, TB>>>(
		THCudaTensor_data(state, output),
		THCudaTensor_data(state, input),
		THCudaTensor_data(state, affine_model),
		h, w, epsilon, radius
	);
	checkCudaError(L);

	bilateral_smooth_kernel<<<(h*w)/TB+1, TB>>>(
		THCudaTensor_data(state, affine_model),
		THCudaTensor_data(state, filtered_affine_model),
		THCudaTensor_data(state, input),
		h, w, f_r, sigma1, sigma2
	);
	checkCudaError(L);
	THCudaTensor_free(state, affine_model);

	reconstruction_best_kernel<<<(h*w)/TB+1, TB>>>(
		THCudaTensor_data(state, input),
		THCudaTensor_data(state, filtered_affine_model),
		THCudaTensor_data(state, filtered_best_output),
		h, w
	);
	checkCudaError(L);
	THCudaTensor_free(state, filtered_affine_model);
	

	luaT_pushudata(L, filtered_best_output, "torch.CudaTensor");

	return 1;
}

static const struct luaL_Reg funcs[] = {
	{"matting_laplacian", matting_laplacian},
	{"smooth_local_affine", smooth_local_affine},
	{NULL, NULL}
};

extern "C" int luaopen_libcuda_utils(lua_State *L) {
	luaL_openlib(L, "cuda_utils", funcs, 0);
	return 1;
}