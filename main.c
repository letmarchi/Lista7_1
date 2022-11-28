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

double *LeVetor(char *nome, int *m){
  FILE *fp = fopen(nome, "r");
  int i, j;
  double *vetor;  

  fscanf(fp, "%d", m);
   
  
  vetor = (double *)malloc( *m*sizeof(double *));
  
  for(i=0;i<*m;i++){
          fscanf(fp, "%lf", &vetor[i]);
  }

  return vetor;
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


double Conjugado(double **M, int m, int n, double *x0, double *d, int p, double *l1){
  int i,j;
  double *r, a, b, *aux, *aux1, ort;
   
  r = (double *)malloc(m*sizeof(double));
  aux = (double *)malloc(m*sizeof(double));
 
  aux1 = (double *)malloc(m*sizeof(double));

  
  r = Residuos(M, m, n, x0);

	for(i=0;i<m;i++){
		for(j=0;j<m;j++){
			if(i==j){
				ort = r[i]*l1[j] + ort;
			} 

		}
	}

  d = Residuos(M, m, n, x0);

  aux1 = r;

  a = (MultVetor(r,r,m))/(MultVetor(r,MultiVetMat(M, r, m, n),m));

  for(i=0;i<m;i++) x0[i] = x0[i] + a*d[i];
  
  aux = MultiVetMat(M, d, m, n);

  l1 =r;
	if(ort == 0){
		printf("Os residuos sao ortogonais\n");
	}

  for(i=0;i<m;i++) r[i] = r[i] - a*aux[i];

  b = (MultVetor(r,r,m))/(MultVetor(aux1,aux1,m));

  for(i=0;i<m;i++) d[i] = d[i] + b*d[i];

  return(NormaVetor(r, m, p));
}

double MaximaDescida(double **M, int m, int n, double *x0, int p, double *l1){

  double *lt, *t, lambda, ort=0;
  
  int i,j;

  t = (double *)malloc(m*sizeof(double));
  
  t = Residuos(M, m, n,x0);

  for(i=0;i<m;i++){
		for(j=0;j<m;j++){
			if(i==j){
				ort = t[i]*l1[j] + ort;
			} 

		}
	}

  lambda=  (MultVetor(t,t,m))/(MultVetor(t,MultiVetMat(M, t, m, n),m));

  l1 =t;

  for(i=0;i<m;i++) x0[i]+=lambda*t[i];

	if(ort == 0){
		printf("Os residuos sao ortogonais\n");
	}

  return(NormaVetor(t, m, p));
}


int main() {
double **M, **C,*b ,*v,*d,*l1, dx, tolerance=1e-8;
int m, n, l, i, it=0, p=0;
FILE *arq;

M = LeMatriz("Matrix.dat",&m, &n);
v = LeVetor("vetor.dat", &l);

l1 = (double*)calloc(m,sizeof(double));

printf("\n Maxima descida\n");
do{
  it++;
  dx = MaximaDescida(M,m,n,v,p,l1);
  printf("It:%d %8.4g\n", it, dx);
  for( i=0; i<m; i++) printf("%11.6g ", v[i]);
  puts("");
}while (dx > tolerance);

free(l1);

l1 = (double*)calloc(m,sizeof(double));

v = LeVetor("vetor.dat", &l);
d = v;
it =0;
printf("\ngradiente conjugado\n");
do{
  it++;
  dx = Conjugado(M,m,n,v,d,p, l1);
  printf("It:%d %8.4g\n", it, dx);
  for( i=0; i<m; i++) printf("%11.6g ", v[i]);
  puts("");
}while (dx > tolerance);

free(l1);

v = LeVetor("vetor.dat", &l);
it =0;
printf("\nGauss\n");
do{
  it++;
  dx = Gauss(M,m,n,v,p);
  printf("It:%d %8.4g\n", it, dx);
  for( i=0; i<m; i++) printf("%11.6g ", v[i]);
  puts("");
}while (dx > tolerance);

free(l1);

printf("Letra c");

M = LeMatriz("Matrix1.dat",&m, &n);
v = LeVetor("vetor.dat", &l);

l1 = (double*)calloc(m,sizeof(double));

printf("\n Maxima descida\n");
do{
  it++;
  dx = MaximaDescida(M,m,n,v,p,l1);
  printf("It:%d %8.4g\n", it, dx);
  for( i=0; i<m; i++) printf("%11.6g ", v[i]);
  puts("");
}while (dx > tolerance);

free(l1);

l1 = (double*)calloc(m,sizeof(double));

v = LeVetor("vetor.dat", &l);
d = v;
it =0;
printf("\ngradiente conjugado\n");
do{
  it++;
  dx = Conjugado(M,m,n,v,d,p, l1);
  printf("It:%d %8.4g\n", it, dx);
  for( i=0; i<m; i++) printf("%11.6g ", v[i]);
  puts("");
}while (dx > tolerance);

free(l1);

v = LeVetor("vetor.dat", &l);
it =0;
printf("\nGauss\n");
do{
  it++;
  dx = Gauss(M,m,n,v,p);
  printf("It:%d %8.4g\n", it, dx);
  for( i=0; i<m; i++) printf("%11.6g ", v[i]);
  puts("");
}while (dx > tolerance);

free(l1);

return 0;
}
