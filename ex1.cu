#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "common/book.h"
using namespace std;


void load_file(vector<double>* vec, string fileName) {
	// string fileName = "data.txt";
	ifstream reader;
	reader.open(fileName, ios::in);
	if (reader.is_open()) {
		char linea[60];
		double valor;
		while (!reader.eof()) {
			reader >> linea;
			if (strlen(&linea[0]) != 0) {
				valor = stod(linea);
				vec->push_back(valor);
			}
		}

	}
	else {
		cout << "No se pudo abrir el fichero..." << endl;
	}
	reader.close();
}

void load_data_array(double* vec, int size,string fileName) {
	// string fileName = "data.txt";
	int count =0;
	// int size = sizeof(vec)/sizeof(*vec);
	ifstream reader;
	reader.open(fileName, ios::in);
	if (reader.is_open()) {
		char linea[60];
		double valor;
		while (!reader.eof()&&count<size) {
			reader >> linea;
			if (strlen(&linea[0]) != 0) {
				valor = stod(linea);
				vec[count] = valor;
				// cout<<"añadio: "<<vec[count]<<endl;
				count++;
				
				// vec->push_back(valor);
			}
		}

	}
	else {
		cout << "No se pudo abrir el fichero..." << endl;
	}
	reader.close();
}

const int N = 16; 
const int blocksize = 16; 
 
__global__ void hello(char *a, int *b) 
{
	a[threadIdx.x] += b[threadIdx.x];
}
////////////////////////////////////////////////////// C++ ///////////////////////////////////////////////////////////////
template <class RandomAccessIterator>
void print_vector(RandomAccessIterator inicio, RandomAccessIterator fin, string titulo) {
	RandomAccessIterator index;
	cout<<"=========="<<titulo<<"=========="<<endl;
	for (index = inicio; index != fin; index++) {
		cout << "indx: " << index - inicio << "->" << inicio[index - inicio]
			<< endl;
	}
}
void print_array(double* array, int size,string titulo) {
	cout<<"=========="<<titulo<<"=========="<<endl;
	// int size = (sizeof(array))/sizeof(*array);
	for (int  i = 0; i < size; i++) {
		cout << "indx: " << i << "->" << array[i]<< endl;
	}
}

void quickSort(vector<double>* array, int left, int right){
	// cout << "usando serial quicksort..." << endl;
	int i = left, j = right;
	double tmp;
	double pivot = (*array)[(left + right) / 2];
		/* PARTITION PART */
		while (i <= j) {
			while ((*array)[i] < pivot)
				i++;
			while ((*array)[j] > pivot)
				j--;
			if (i <= j) {
				tmp = (*array)[i];
				(*array)[i] = (*array)[j];
				(*array)[j] = tmp;
				i++;
				j--;
			}
		}
	if (left < j){ quickSort(array, left, j); }			
	if (i < right){ quickSort(array, i, right); }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////777
 
void cuda_info(){
	cudaDeviceProp  prop;

    int count;
    // HANDLE_ERROR( cudaGetDeviceCount( &count ) );
	cudaGetDeviceCount( &count );
    for (int i=0; i< count; i++) {
        // HANDLE_ERROR( cudaGetDeviceProperties( &prop, i ) );
		cudaGetDeviceProperties( &prop, i ) ;
        printf( "   --- General Information for device %d ---\n", i );
        printf( "Name:  %s\n", prop.name );
        printf( "Compute capability:  %d.%d\n", prop.major, prop.minor );
        printf( "Clock rate:  %d\n", prop.clockRate );
        printf( "Device copy overlap:  " );
        if (prop.deviceOverlap)
            printf( "Enabled\n" );
        else
            printf( "Disabled\n");
        printf( "Kernel execution timeout :  " );
        if (prop.kernelExecTimeoutEnabled)
            printf( "Enabled\n" );
        else
            printf( "Disabled\n" );

        printf( "   --- Memory Information for device %d ---\n", i );
        printf( "Total global mem:  %ld\n", prop.totalGlobalMem );
        printf( "Total constant Mem:  %ld\n", prop.totalConstMem );
        printf( "Max mem pitch:  %ld\n", prop.memPitch );
        printf( "Texture Alignment:  %ld\n", prop.textureAlignment );

        printf( "   --- MP Information for device %d ---\n", i );
        printf( "Multiprocessor count:  %d\n",
                    prop.multiProcessorCount );
        printf( "Shared mem per mp:  %ld\n", prop.sharedMemPerBlock );
        printf( "Registers per mp:  %d\n", prop.regsPerBlock );
        printf( "Threads in warp:  %d\n", prop.warpSize );
        printf( "Max threads per block:  %d\n",
                    prop.maxThreadsPerBlock );
        printf( "Max thread dimensions:  (%d, %d, %d)\n",
                    prop.maxThreadsDim[0], prop.maxThreadsDim[1],
                    prop.maxThreadsDim[2] );
        printf( "Max grid dimensions:  (%d, %d, %d)\n",
                    prop.maxGridSize[0], prop.maxGridSize[1],
                    prop.maxGridSize[2] );
        printf( "\n" );
    }
}
int main()
{
	// char a[N] = "Hello \0\0\0\0\0\0";
	// int b[N] = {15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
 
	// char *ad;
	// int *bd;
	// const int csize = N*sizeof(char);
	// const int isize = N*sizeof(int);
 
	// printf("%s", a);
 
	// cudaMalloc( (void**)&ad, csize ); 
	// cudaMalloc( (void**)&bd, isize ); 
	// cudaMemcpy( ad, a, csize, cudaMemcpyHostToDevice ); 
	// cudaMemcpy( bd, b, isize, cudaMemcpyHostToDevice ); 
	
	// dim3 dimBlock( blocksize, 1 );
	// dim3 dimGrid( 1, 1 );
	// hello<<<dimGrid, dimBlock>>>(ad, bd);
	// cudaMemcpy( a, ad, csize, cudaMemcpyDeviceToHost ); 
	// cudaFree( ad );
	// cudaFree( bd );
	
	// printf("%s\n", a);
	
	// cuda_info();
	// cudaDeviceProp  prop;
	// int nDevices=-1;
	// // cudaGetDeviceCount( &count );
	// cudaError_t err = cudaGetDeviceCount(&nDevices);
  	// if (err != cudaSuccess) printf("%s\n", cudaGetErrorString(err));
	// cout<<"no_devices: "<<nDevices<<endl;
	// cudaGetDeviceProperties( &prop, 1 ) ;
	// cout<<"nombre: "<< prop.name<<endl;
	
	// cuda_info();

	//host_variables

	int size = 1000;
	double* data = new double[size];
	data[50]=108;

	//load variable
	string file="1000Int.txt";
	load_data_array(data,size,file);
	// print_array(data,size,"data array raw");
	// size = *(&data + 1) - data;
	// load_file(&vec,file);
	// size =vec.size() ;
	// cout<<"size: "<<data[50]<<endl;
	// quickSort(&vec,0,vec.size()-1);
	// print_vector(vec.begin(),vec.end(),"ordenamiento secuencial");

	//device_variables
	double* d_data;

	// vector<double>* result;
	
	HANDLE_ERROR(cudaMalloc( (void**)&d_data, size * sizeof(double) ));
	HANDLE_ERROR(cudaMemcpy(d_data,data,size * sizeof(double),cudaMemcpyHostToDevice));

	double result[size];
	HANDLE_ERROR(cudaMemcpy(result,d_data,size * sizeof(double),cudaMemcpyDeviceToHost));
	print_array(&result[0],size,"data array device");

	cudaFree(d_data);
	
	return EXIT_SUCCESS;
}