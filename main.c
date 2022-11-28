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

double *MultiVetMat(double **M, double *x, int m, int n){
  int i, j, k;
  double *r;

  r=malloc(m*sizeof(double *));

  for(i=0; i<m; i++){
    for(j=0; j<m; j++){
      r[i]+=M[i][j]*x[j];
    }
  }
  return r;
}

double MultVetor(double *x, double *y, int t){
  int i;
  double l=0;
  for(i=0;i<t;i++){
    l+=x[i]*y[i];
  }
  return l;
}

double *Residuos(double **M, int m, int n, double *x){
  int i;
  double *r,*s;

  r=malloc(m*sizeof(double *));
  s=MultiVetMat(M, x, m, n);
  for(i=0;i<m;i++)
    r[i]=M[i][n-1]-s[i];

  return r;
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

double Gradiente(double **Matriz, int m, int n, double *x0, int p){
  int i,j;
  double sum, lambda, *r;

  r=(double *) malloc(m *sizeof(double));
  r=Residuos(Matriz, m, n, x0);
  lambda=MultVetor(r,r,m)/MultVetor(r, MultiVetMat(Matriz, r, m, n), m);
  for(i=0;i<m;i++)
    x0[i]+=lambda*r[i];
  return NormaVetor(r, m, p);
}

double Conjugado(double **Matriz, int m, int n, double *x0, double *d, int p){
  int i,j;
  double a, *Ad, rr, *r, B, *r1;

  r=r1=(double *) malloc(m *sizeof(double));
  r=r1=Residuos(Matriz, m, n, x0); 
  Ad=MultiVetMat(Matriz, d, m, n);
  rr=MultVetor(r,r,m);
  a=rr/MultVetor(d, Ad, m);
  
  for(i=0;i<m;i++){
    x0[i]+=a*d[i];
    r[i]-=a*Ad[i];
  }
  B=MultVetor(r,r,m)/rr;
  for(i=0;i<m;i++){
    d[i]=r[i]+B*d[i];
  }
  return NormaVetor(r, m, p);
}


int main() {
double **MA, **M, *VI, dx, *d, *r, tolerance=1e-8, w;
int m, n, l, i, it=0, p=0;
FILE *arq;

MA = LeMatriz("Matrix.dat",&m, &n);
printf("Matriz Aumentada\n");
ImprimeMatriz(MA, &m, &n);
SeparaMatriz(MA, m, n, &M, &VI);
printf("\nMatriz \n");
ImprimeMatriz(M, &m, &m);
printf("\nVetor \n");
ImprimeVetor(VI, m);

d=r=Residuos(M, m, n, VI);
{
  do{
  it++;
  //dx=Jacobi(M, m, n, v, p);
  //dx=Gauss(M, m, n, v, p);
  //dx = Relaxacao(M, m, n, VI, w, p);
  //dx=Gradiente(M, m, n, v, p);
  dx=Conjugado(M, m, n, VI, d, p);
  printf("%d %8.4g ", it,dx);
  for( i=0; i<m; i++) printf("%11.6g ", VI[i]);
  puts("");
  } while (dx > tolerance);
  
}
return 0;
}


