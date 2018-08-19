// Sort2013.cpp: define el punto de entrada de la aplicaci�n de consola.




// opmp1.cpp: define el punto de entrada de la aplicaci�n de consola.

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <pthread.h>
#include <chrono>
//#include "mpi.h"

using namespace std;

void calcular_resultados(char* titulo, double t_serial, double t_inicio,
	double t_final) {
	int no_nucleos = omp_get_num_procs();

	double time_parallel = t_final - t_inicio;
	double speedUp = t_serial / time_parallel;
	double eficiencia = speedUp / no_nucleos;
	cout << endl;
	cout << "=============resultados: " << titulo << " ============" << endl;
	cout << "tiempo: " << time_parallel << " sec. Speedup: " << speedUp
		<< ". Eficiencia: " << eficiencia << endl;
}


template <class RandomAccessIterator>
void print_vector(RandomAccessIterator inicio, RandomAccessIterator fin, string titulo) {
	RandomAccessIterator index;
	cout<<"=========="<<titulo<<"=========="<<endl;
	for (index = inicio; index != fin; index++) {
		cout << "indx: " << index - inicio << "->" << inicio[index - inicio]
			<< endl;
	}
}

bool isSorted(vector<double>* vec) {
	int size = vec->size();
	// Array has one or no element
	if (size == 0 || size == 1) {
		cout << "order" << endl;
		return true;
	}
	for (int i = 1; i < size; i++) {
		// Unsorted pair found
		if ((*vec)[i - 1] >(*vec)[i]) {
			cout << "no order" << endl;
			return false;
		}
	}
	cout << "order" << endl;
	// No unsorted pair found
	return true;
}

//------------declarations
void swapChunks(vector<double>* array,vector<double>* tmp_array,int left, int right,
			 int global_size, int no_threads, vector<int>* pivots,vector<int>* inicios);
void chunkPartition(vector<double>* array,int left, int right, double pivot, int* total_menores, int my_rank,vector<int>*pivotes);
void calcular_inicios(vector<int>* pivots,int local_size,vector<int>* inicios);
//--------------------------------

void distributedQuickSort(vector<double>* array,int no_threads,int global_size,int left, int right, 
						double global_pivot, vector<int>* pivotes,int* total_menores,
						vector<double>* tmparray,vector<int>* inicios){
	int my_rank = omp_get_thread_num();
	int local_size = global_size / no_threads;
	int local_begin = my_rank * local_size;
	int local_end = local_begin + local_size -1;
	
	// #pragma omp critical
	// cout << "soy el thread: " << my_rank << ", inicio: " << local_begin << " fin: " << local_end << endl;
	
	chunkPartition(array,local_begin, local_end, global_pivot,total_menores,my_rank,pivotes);
	// (*pivotes)[my_rank]=indice_mitad_local;
	swapChunks( array, tmparray, left,  right,
			  global_size,  no_threads, pivotes,inicios);
}

void chunkPartition(vector<double>* array,int left, int right, double pivot, int* total_menores, int my_rank,vector<int>*pivotes){
	int i = left, j = right;
	double tmp;
	// double pivot = (*array)[(left + right) / 2];
		/* PARTITION PART */
		while (i <= j) {
			while ((*array)[i] <= pivot && i<= j)
				i++;
			while ((*array)[j] > pivot && i<= j)
				j--;
			if (i <= j) {
				tmp = (*array)[i];
				(*array)[i] = (*array)[j];
				(*array)[j] = tmp;
				i++;
				j--;
			}
		}
		*total_menores += (j>=left)?j-left+1:0;

	(*pivotes)[my_rank]= (j>=left)?j:-1;
}


void calcular_inicios(vector<int>* pivots,int local_size,vector<int>* inicios){
	
	cout<<"calculando inicios by thread: "<<omp_get_thread_num()<<endl;
	print_vector(pivots->begin(),pivots->end(),"pivotrs");

	int inicio_anterior=0;
	
	int mitad = inicios->size()/2;
	
	int no_menores_anterior, no_menores_local,inicio_local;
	
	// (*pivots)[0]=0;
	(*inicios)[0]=0;
	cout<<"mitad de inicios-----: "<< mitad <<endl;
	for(int i=1;i<mitad;i++){
		// inicio_local = 0;
		 if((*pivots)[i-1]>=0){
				inicio_local = i*local_size;
				inicio_anterior= inicio_local-local_size;
				//  no_menores_local = (*pivots)[i]-inicio_local+1;
				no_menores_anterior = (*pivots)[i-1]-inicio_anterior+1;
		} 
		else{no_menores_anterior=0;} 
		
		(*inicios)[i]=no_menores_anterior+(*inicios)[i-1];
	}
	int no_mayores_anterior, no_mayores_local,fin_local,fin_anterior;
	int rank_ultimo = mitad-1;
	if((*pivots)[rank_ultimo]>=0){ //existen menores en el ultimo thread?
		no_menores_anterior = (*pivots)[rank_ultimo]-rank_ultimo*local_size+1;
		 cout<<"no menores ultimo-----: "<< no_menores_anterior <<endl;
	}else no_menores_anterior=0;
	
	(*inicios)[mitad]=(*inicios)[mitad-1]+no_menores_anterior;

	int local_rank=0;
	for(int i=mitad+1;i<2*mitad;i++){
		if((*pivots)[i-1]>=0){
				local_rank= (i-mitad);
				inicio_local = local_rank*local_size;
				// inicio_anterior= inicio_local-local_size;
				fin_anterior = inicio_local-1;
				//  no_menores_local = (*pivots)[i]-inicio_local+1;
				no_mayores_anterior = fin_anterior-(*pivots)[(i-mitad)-1];
				
		} 
		else no_mayores_anterior=local_size;
		(*inicios)[i]=no_mayores_anterior+(*inicios)[i-1];
	}
	
	// print_vector((*inicios).begin(),(*inicios).end(),"inicios");
}
void swapChunks(vector<double>* array,vector<double>* tmp_array,int left, int right,
			 int global_size, int no_threads, vector<int>* pivots,vector<int>* inicios){
	
	int my_rank = omp_get_thread_num();
	int local_size = global_size / no_threads;
	int local_begin = my_rank * local_size;
	int local_end = local_begin + local_size -1;
	int no_menores_local = (*pivots)[my_rank]-local_begin+1;
	int no_mayores_local = ((*pivots)[my_rank]>=0)?local_end-(*pivots)[my_rank]:local_size;
	// #pragma omp master{
	// 	for(int i=0;i<no_menores_local;i++){
	// 		(*tmp_array)[i]=(*array)[i];
	// 	}
	// }
	// vector<int> inicios(no_threads);
	#pragma omp single
		calcular_inicios(pivots, local_size,inicios);
	 #pragma omp barrier
	// #pragma omp critical
	// print_vector(inicios.begin(),inicios.end(),"inicios");
	// #pragma omp critical
	// cout << "soy el thread: " << my_rank <<", inicio:"<<local_begin<<", fin: "<<local_end <<", tengo: " <<  no_menores_local
	// << " menores " <<"inicio en: "<<(*inicios)[my_rank]<< endl;

//copiar menores a temporal
	if(no_menores_local>0){
		// cout << "soy el thread: "<<endl;
		int inicio_en_tmp =(*inicios)[my_rank];
		for(int i = 0;i<no_menores_local;i++){
			(*tmp_array)[inicio_en_tmp+i]=(*array)[local_begin+i];
		}
	}

	#pragma omp barrier
	int mitad=inicios->size()/2;
		//copiar mayores a tmp
		int inicio_en_tmp =(*inicios)[my_rank + mitad];
		
		for(int i = 0;i<no_mayores_local;i++){
			(*tmp_array)[inicio_en_tmp+i]=(*array)[local_end-i];
		}
	
	
}

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

void quickSort(vector<double>* array, int left, int right){
	cout << "usando serial quicksort..." << endl;
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

void quickSort_parallel_internal(vector<double>* array, int left, int right, int cutoff);
void quickSort_parallel(vector<double>* array, int lenArray, int numThreads){

	int cutoff = 10000;

	#pragma omp parallel num_threads(numThreads)
	{	
		#pragma omp single nowait
		{
			quickSort_parallel_internal(array, 0, lenArray-1, cutoff);	
		}
	}	

}

void quickSort_parallel_internal(vector<double>* array, int left, int right, int cutoff) 
{
	int i = left, j = right;
	double tmp;
	double pivot = (*array)[(left + right) / 2];
	// int my_rank = omp_get_thread_num();
	// cout<<"i'am the thread: "<<my_rank<<", inicio: "<<left
	// <<", fin: "<<right<<endl;
	{
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
				// cout<<"intercambio: i="<<i<<", j="<<j<<endl;
				i++;
				j--;
				
			}
		}

	}


	if ( ((right-left)<cutoff) ){
		if (left < j){ 
			// quickSort_parallel_internal(array, left, j, cutoff); 
			quickSort(array, left, j); 
		}			
		if (i < right){ 
			// quickSort_parallel_internal(array, i, right, cutoff); 
			quickSort(array, i, right); 
		}

	}else{
		#pragma omp task 	
		{ quickSort_parallel_internal(array, left, j, cutoff); }
		#pragma omp task 	
		{ quickSort_parallel_internal(array, i, right, cutoff); }
		// #pragma omp task 	
		// { quickSort_parallel_internal(array, left, j, cutoff); }
		// #pragma omp task 	
		// { quickSort_parallel_internal(array, i, right, cutoff); }		
		// #pragma omp task 	
		// { quickSort_parallel_internal(array, left, j, cutoff); }
		// #pragma omp task 	
		// { quickSort_parallel_internal(array, i, right, cutoff); }
	}

}

int main(int argc, char* argv[]) {
	cout<<"corriendo en ubuntu.."<<endl;
 	vector<double> vec;
string file="dataInt.txt";
	load_file(&vec, file);
//Serial QuickSort
// 	 cout<<"cargado listo...tam: "<<vec.size()<<endl;
 double time_serial_inicio,time_serial_fin,time_parallel_inicio,time_parallel_fin;
// time_serial_inicio = omp_get_wtime();
// quickSort(&vec,0,vec.size()-1);
// time_serial_fin = omp_get_wtime();
// 	 cout<<"tiempo serial: "<<time_serial_fin-time_serial_inicio<<endl;
	
	//::::::Paralell quickSort load
	vec.clear();
	load_file(&vec,file);
	// vec.pop_back();

	cout<<"cargado listo...tam: "<<vec.size()<<endl;
	int no_threads = 8;
	using namespace std;
	using namespace std:: chrono;
	using clk=chrono::high_resolution_clock;


	//::::::Paralell quickSort
	time_parallel_inicio = omp_get_wtime();

	double pivot = vec[vec.size()/2];
	cout<<"global_pivot: "<<pivot<<endl;
 	vector<int> pivotes(no_threads);
	vector<int> inicios(no_threads*2);
	vector<double> vector_desplazado(vec.size());
 	int total_menores = 0;
#pragma omp parallel num_threads(no_threads)
	distributedQuickSort(&vec,no_threads,vec.size(),0,vec.size(),pivot,&pivotes,&total_menores,&vector_desplazado,&inicios);

	// auto t1 = clk::now();
	//  quickSort_parallel(&vec,vec.size(),no_threads);
	//  auto t2 = clk::now();
	//  auto diff = duration_cast<microseconds>(t2-t1);
	time_parallel_fin = omp_get_wtime();
	cout<<"tiempo paralelo: "<<time_parallel_fin-time_parallel_inicio<<endl;
// //=======resultado
//  	print_vector(&vec);

 	print_vector(vec.begin(),vec.end(),"array raw");
	print_vector(vector_desplazado.begin(),vector_desplazado.end(),"array tmp");
print_vector(inicios.begin(),inicios.end(),"inicios");

// print_vector(pivotes.begin(),pivotes.end(),"mis picotessssssssss");
pivotes.clear();
inicios.clear();
vec.clear();
vector_desplazado.clear();
	return 0;
}