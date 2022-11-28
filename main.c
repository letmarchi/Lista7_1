#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double **LeMatriz(char nome[20], int *m, int *n){
  double **a;
  int i, j;
  FILE *arq;
  
  arq=fopen(nome, "r");
  fscanf(arq, "%d", m);
  fscanf(arq, "%d", n);
  a = (double **) malloc(*m*sizeof(double *));
  for (i=0;i<*m;i++) a[i]= (double *) malloc((*n)*sizeof(double));
  for (i=0;i<*m;i++){
    for (j=0;j<*n;j++){
      fscanf(arq, "%lf", &a[i][j]);
    }
  }
  return a;
}

void **ImprimeMatriz(double **M, int *m, int *n){
  int i, j;
  
  for (i=0;i<*m;i++){
    for (j=0;j<*n;j++){
      printf("%lf ", M[i][j]);
    }
  puts("");
  }
  
}

void SeparaMatriz(double **M, int m, int n, double ***S, double **b)
{
		int i, j;

	double lvd;

	*S=malloc (m*sizeof(double *));
		for (i=0 ; i<m ; i++)
		{
			(*S)[i]=malloc(m*sizeof(double));
		}
		*b=malloc(m*sizeof(double));
		for(i=0 ; i< m ; i++)
		{
				(*b)[i] = M[i][n-1];
		}

		for(i=0 ; i<m ; i++)
		 {
				for( j=0; j<m ; j++)
				{
						(*S)[i][j]=M[i][j];
				}
		 }

		return;
}

void ImprimeVetor(double *v, int j){
  int i;
  for(i=0; i<j; i++) printf("v[%d] = %g\n", i, v[i]);
}

double NormaMatriz(double **v, int m, int n, int p){
  double max=0, aux;
  int i,j;

  if(p==0){
    for(i=0;i<m;i++){
      aux=0;
      for(j=0;j<n;j++){
        aux+=fabs(v[i][j]);
      }
      if(aux>max) max=aux;
    }
  }
  else if(p==1){
    for(i=0;i<n;i++){
      aux=0;
      for(j=0;j<m;j++){
        aux+=fabs(v[j][i]);
      }
      if(aux>max) max=aux;
    }
  }
  else{
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
          max+=pow(v[i][j],2);
        }
    }
    max=sqrt(max);
  }
  return max;
}

double NormaVetor(double *v, int m, int p){
  
  int i;
  double max=0;

  if(p==0){
    for(i=0;i<m;i++) {
      if(fabs(v[i])>max) max=fabs(v[i]);
    }
  }
  else{
    for(i=0;i<m;i++) {
      max+=pow(fabs(v[i]),p); 
     }
    max=pow(max,1.0/p); 
  }
  return max;
}

double Jacobi(double **Matriz, int m, int n, double *x0, int p){
  int i,j;
  double sum, *v, *r;

  v = (double *) malloc(m *sizeof(double));

  for(i=0;i<m;i++){
    sum=0;
    for(j=0;j<n-1;j++){
      if(i!=j) sum+=Matriz[i][j]*x0[j];
    }
    v[i]=(Matriz[i][n-1]-sum)/Matriz[i][i];
  }
  r = (double *) malloc(m *sizeof(double));
  for(i=0;i<m;i++){
    r[i]=x0[i]-v[i];
  }
  memcpy(x0, v, m*sizeof(double));
  return NormaVetor(r, m, p);
}

double Gauss(double **Matriz, int m, int n, double *x0, int p){
  int i,j;
  double sum, *v, *r;

  v = (double *) malloc(m *sizeof(double));
  memcpy(v, x0, m*sizeof(double));
  for(i=0;i<m;i++){
    sum=0;
    for(j=0;j<n-1;j++){
      if(i!=j) sum+=Matriz[i][j]*x0[j];
    }
    x0[i]=(Matriz[i][n-1]-sum)/Matriz[i][i];
  }
  r = (double *) malloc(m *sizeof(double));
  for(i=0;i<m;i++){
    r[i]=x0[i]-v[i];
  }
  return NormaVetor(r, m, p);
}

double Relaxacao(double **Matriz, int m, int n, double *x0, double omega, int p){
  int i,j;
  double sum, *v, *r;

  v = (double *) malloc(m *sizeof(double));
  memcpy(v, x0, m*sizeof(double));
  for(i=0;i<m;i++){
    sum=0;
    for(j=0;j<n-1;j++){
      if(i!=j) sum+=Matriz[i][j]*x0[j];
    }
    x0[i]=((1-omega)*x0[i]) + omega*(Matriz[i][n-1]-sum)/Matriz[i][i];
  }
  r = (double *) malloc(m *sizeof(double));
  for(i=0;i<m;i++){
    r[i]=x0[i]-v[i];
  }
  return NormaVetor(r, m, p);
}

int main() {
double **MA, **M, *VI, dx, tolerance=1e-7, w;
int m, n, l, i, it=0, p=0;
FILE *arq;

MA = LeMatriz("Matrix.dat",&m, &n);
printf("Matriz Aumentada\n");
ImprimeMatriz(MA, &m, &n);
SeparaMatriz(MA, m, n, &M, &VI);
printf("\nMatriz \n");
ImprimeMatriz(M, &m, &m);
ImprimeVetor(VI, m);
//for( i=0; i<m; i++) printf("%11.6g ", VI[i]);
//puts("");

//w=1;
//while(w<2){
//  do{
//  it++;
  //dx=Jacobi(M, m, n, v, p);
  //dx=Gauss(M, m, n, v, p);
//  dx = Relaxacao(M, m, n, VI, w, p);
  //printf("%d %8.4g ", it,dx);
  //for( i=0; i<m; i++) printf("%11.6g ", v[i]);
  //puts("");
//  } while (dx > tolerance);
//  printf("%d %8.4g %lf ", it,dx,w);
//  for( i=0; i<m; i++) printf("%11.6g ", VI[i]);
//  puts("");
//  it=0;
//  for( i=0; i<m; i++) VI[i]=0;
//  w+=0.1;
//}
return 0;
}


