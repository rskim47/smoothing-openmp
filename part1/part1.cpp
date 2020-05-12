#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <time.h>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define	INDEX(i,j,n) i+j*n // Array Indexing  

using namespace std;

// Initialize (Subprogram) 
// x - Array, n - # of elements 
void initialize(float* x, long n)
{
  int i,j;
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j) schedule(dynamic,100)	
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
		  x[INDEX(i,j,n)] = (float)abs( i%11 - j%5 ) / (float)( i%7 + j%3 + 1);
		}
	}
  #endif
}

// Smooth (Subprogram)
// x - Input Array, y - Output Array, n - # of elements, (a,b,c) - Smoothing Constants 
void smooth(float* x, float* y, long n, float a, float b, float c)
{	
  int i,j; 
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j) schedule(dynamic,5000)	
  for(i=1; i<n-1; i++){ // Calculating the Inner Loops 
	 	for(j=1; j<n-1; j++){
			y[INDEX(i,j,n)] = a*(x[INDEX((i-1),(j-1),n)] + x[INDEX((i-1),(j+1),n)] + x[INDEX((i+1),(j-1),n)] + x[INDEX((i+1),(j+1),n)]) + b*(x[INDEX((i-1),j,n)] + x[INDEX((i+1),j,n)] + x[INDEX(i,(j-1),n)] + x[INDEX(i,(j+1),n)]) + c*x[INDEX(i,j,n)];
		}
	}	
	#endif
}

// Count (Subprogram)
// array - Input array, n - # of elements, t - Threshold, count - # of elements below threshold 
void count(float* array, long n, float t,long long &count)
{
  int i,j;
  int tempCount = 0;
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j) schedule(dynamic,5000) reduction (+:tempCount)
  for(i=1; i<n-1; ++i){
		for(j=1; j<n-1; ++j){
			if(array[INDEX(i,j,n)] < t){ // Checking Threshold Values 
				tempCount = tempCount + 1;   
			}
		}
	}
	#endif
  count = tempCount;
}

double get_wall_time() {
  struct timeval time; 
  gettimeofday(&time,NULL);
  return (double)time.tv_sec + (double)time.tv_usec * .000001; 
}

int main(){
  printf("There are %d processors available\n", omp_get_num_procs());
  omp_set_num_threads(8);
	#ifdef _OPENMP
	#pragma omp parallel
	{
		printf("This is thread %d\n", omp_get_thread_num());
	}
	#endif 
	
	// Timing variable
	double start,end,time_initialize,time_smooth,time_count_x,time_count_y,time_alloc_x,time_alloc_y;
	
	// Initialize Comstants 
	float a = .05;
	float b = .1;
	float c = .4;

  // Initialize Threshold
  float t = .1; 

	// Array Elements in Row/Columns  
	long n = 98306; // 
	
	// Elements 
	long long count_x = 0;
	long long count_y = 0;

	// Declaring Array X & Y - Dynamic Allocation 
  start = get_wall_time();
	float* x = (float*) malloc(n*n*sizeof(long long)); 
  end = get_wall_time();
  time_alloc_x = end - start;
  start = get_wall_time();
	float* y = (float*) malloc(n*n*sizeof(long long));
  end = get_wall_time();
  time_alloc_y = end - start; 

	// Initialize 
	start = get_wall_time();
  #ifdef _OPENMP
	//#pragma omp parallel
  { 
	  initialize(x,n);
  }
	#endif
	end = get_wall_time();
	time_initialize = end - start;
	 
	// Smooth 
	start = get_wall_time();
	#ifdef _OPENMP
	//#pragma omp parallel 
	{
		smooth(x,y,n,a,b,c);
  }
	#endif
	end = get_wall_time();
 	time_smooth = end-start;

	// Count X below threshold 
	start = get_wall_time();
	#ifdef _OPENMP
	//#pragma omp parallel
	{
		count(x,n,t,count_x);
	}
	#endif 
	end = get_wall_time();
	time_count_x = end-start;

	// Count Y below threshold 
	start = get_wall_time();
	#ifdef _OPENMP
	//#pragma omp parallel
	{
		count(y,n,t,count_y);
	}
	#endif
	end = get_wall_time();
	time_count_y = end - start;
	
	// Result Calculation 
	long Memory = 4*n*n*pow(10,-9);
	float fraction_x = (float) count_x/(n*n);
	float fraction_y = (float) count_y/(n*n);		

	// Summary of Results 
  cout << "Summary" << endl;
  cout << "-------" << endl; 
	cout << "Number of elements in row/column	 :: " << n << endl;
	cout << "Number of inner elements in row/column   :: " << n-2 << endl;
	cout << "Total number of elements		 :: " << n*n<< endl;
	cout << "Total number of inner elements	         :: " << (n-2)*(n-2) << endl;
	cout << "Memory (GB) used per array		 :: " << Memory << endl;
	cout << "Threshold	                         :: " << t << endl;
	cout << "Smoothing constant (a,b,c)		 :: " << a << ", " << b << ", " << c << endl;
	cout << "Number   of elements below threshold (X) :: " << count_x << endl;
	cout << "Fraction of elements below threshold     :: " << fraction_x<< endl;
	cout << "Number   of elements below threshold (Y) :: " << count_y << endl;
	cout << "Fraction of elements below threshold     :: " << fraction_y << endl;
  cout << endl;
  cout << "Action           ::   time/s     Time resolution = " << (double)1.0/CLOCKS_PER_SEC << endl;  
	cout << "------" << endl; 
  cout << "CPU: Alloc-X     :: " << (double)time_alloc_x << endl; 
  cout << "CPU: Alloc-Y     :: " << (double)time_alloc_y << endl; 
  cout << "CPU: Init-X       :: " << (double)time_initialize << endl;
  cout << "CPU: Smooth       :: " << (double)time_smooth << endl; 
  cout << "CPU: Count-X     :: " << (double)time_count_x << endl;
  cout << "CPU: Count-Y     :: " << (double)time_count_y << endl; 
	
	// Deallocating arrays
	delete[] x;
	delete[] y;
}
