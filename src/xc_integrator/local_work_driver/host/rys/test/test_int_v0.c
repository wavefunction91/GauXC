#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>

#include "rys_integral.h"

uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

int main(int argc, char **argv) {

  if(argc < 4) {
    printf("Correct Usage: ./exe.x <shell file> <points file> <nruns>\n");
    return 1;
  }

  shells *shell_list = NULL;
#if 0
  point C[2];

  // the grid point
  C[0].x = -1.46097;
  C[0].y = 0.186863;
  C[0].z = -0.00156859;

  C[1].x = -1.46097;
  C[1].y = 0.186863;
  C[1].z = -0.00156859;
#endif
  
  // Read in shells
  FILE *fin = fopen(argv[1], "r");

  int n = 0;
  int nn = 0;
  fscanf(fin, "%d", &n);
  shell_list = (shells*) malloc(n * sizeof(shells));
  
  for(int i = 0; i < n; ++i) {
    double x, y, z;
    int m, L;
    
    fscanf(fin, "%lf,%lf,%lf,%d", &x, &y, &z, &L);
    fscanf(fin, "%d", &m);

    nn += ((L + 1) * (L + 2) / 2);
    
    shell_list[i].origin.x = x;
    shell_list[i].origin.y = y;
    shell_list[i].origin.z = z;
    shell_list[i].m = m;
    shell_list[i].L = L;

    shell_list[i].coeff = (coefficients*) malloc(m * sizeof(coefficients));
    
    for(int j = 0; j < m; ++j) {
      double a, c;
      fscanf(fin, "%lf,%lf", &a, &c);

      shell_list[i].coeff[j].alpha = a;
      shell_list[i].coeff[j].coeff = c;
    }
  }
  	 
  fclose(fin);


  // Read in points
  fin = fopen( argv[2], "r" );

  int npts = 0;
  fscanf(fin,"%d",&npts);
  point* C = (point*) malloc(npts * sizeof(point));
  for(int p = 0; p < npts; ++p) {
    double x,y,z;
    fscanf(fin,"%lf,%lf,%lf",&x,&y,&z);
    C[p].x = x;
    C[p].y = y;
    C[p].z = z;
  }

  fclose(fin);

  double *matrix = (double*) malloc(npts * nn * nn * sizeof(double));
  memset((void*) matrix, 0, npts * nn * nn * sizeof(double));

  int runs = atoi(argv[3]);

  long long t0, t1, sum = 0;
  
  for(int r = 0; r < runs; ++r) {
    t0 = rdtsc();
    compute_integral(n, shell_list, npts, C, matrix);
    t1 = rdtsc();

    sum += (t1 - t0);
  }

  printf("Exec: %lf\n", sum / ((double) (runs * 1.0)));

#ifdef DEBUG

  for(int p = 0; p < npts; ++p) {
    for(int j = 0; j < nn; ++j) {
      for(int i = 0; i < nn; ++i) {
	printf("%lf\t", *(matrix + nn * nn * p + nn * j + i));
      }
      printf("\n");
    }
    printf("\n\n");
  }

#endif

  for(int i = 0; i < n; ++i) {
    free(shell_list[i].coeff);
  }
  free(shell_list);
  free(C);

  free(matrix);
  
  return 0;
}
