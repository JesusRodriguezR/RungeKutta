#include<stdio.h>
#include<omp.h>
#include<math.h>
#define MPI 3.14159265358979323846
#define MEULER 2.71828182846
void ec1();
void ec2();
void ec3();
void ec4();
double r1,r2,r3,r4;
long N=50000;
   FILE *fptr;
   FILE *fptr2;
   FILE *fptr3;
   FILE *fptr4;

void main()
{
    int i, nthreads;
    const double startTime = omp_get_wtime();
   
   
    

     
        
           
                ec1();
            
                ec2();
            
                ec3();
            
                ec4();
        

    
    const double endTime = omp_get_wtime();
    printf("tomo (%lf) segundos\n", (endTime - startTime));
    printf("Res: %f %f %f %f", r1,r2,r3,r4);
    fclose(fptr);
    fclose(fptr2);
    fclose(fptr3);
    fclose(fptr4);
}

void ec1()
{
    double h,t,w,ab;
    double w0=MPI/4,a=0,b=MPI, k1,k2,k3,k4;
    int i;
    fptr=fopen("f1.txt","w");
    w=w0;
    h=(b-a)/N;
    for(i=0;i<N;i++){
        t=a+(h*i);
        k1=h*(t*pow(MEULER,3.0*t)-2.0*w);
        k2=h*((t+h/2.0)*pow(MEULER,3.0*(t+h/2.0))-2.0*(w+k1/2.0));
        k3=h*((t+h/2.0)*pow(MEULER,3.0*(t+h/2.0))-2.0*(w+k2/2.0));
        k4=h*((t+h)*pow(MEULER,3.0*(t+h))-2.0*(w+k3));
        w=w+(k1+2.0*k2+2.0*k3+k4)/6.0;
        fprintf(fptr, "%f\t %f \t \n", t+h, w);
    }
    r1=w;
}

void ec2()
{
    double h,t,w,ab;
    double w0=MPI/4,a=0,b=MPI, k1,k2,k3,k4;
    int i;
    fptr2=fopen("f2.txt","w");
    w=w0;
    h=(b-a)/N;
    for(i=0;i<N;i++){
        t=a+(h*i);
        k1=h*(1+pow((t-w),2));
        k2=h*(1+pow((t+h/2.0)-w-(k1/2.0),2));
        k3=h*(1+pow((t+h/2.0)-w-(k2/2.0),2));
        k4=h*(1+pow((t+h-w-k3),2));
        w=w+(k1+2.0*k2+2.0*k3+k4)/6.0;
        fprintf(fptr2, "%f\t %f \t \n", t+h, w);
    }
    r2=w;
}

void ec3()
{
    double h,t,w,ab;
    double w0=MPI/4,a=0,b=MPI, k1,k2,k3,k4;
    int i;
    fptr3=fopen("f3.txt","w");
    w=w0;
    h=(b-a)/N;
    for(i=0;i<N;i++){
        t=a+(h*i);
        k1=h*(1+w/t);
        k2=h*(1+(w+k1/2.0)/(t+h/2.0));
        k3=h*(1+(w+k2/2.0)/(t+h/2.0));
        k4=h*(1+(w+k3)/(t+h));
        w=w+(k1+2.0*k2+2.0*k3+k4)/6.0;
        fprintf(fptr3, "%f\t %f \t \n", t+h, w);
    }
    r3=w;
}

void ec4()
{
    double h,t,w,ab;
    double w0=MPI/4,a=0,b=MPI, k1,k2,k3,k4;
    int i;
    fptr4=fopen("f4.txt","w");
    w=w0;
    h=(b-a)/N;
    for(i=0;i<N;i++){
        t=a+(h*i);
        k1=h*(cos(2*t*w)+sin(3*t*w));
        k2=h*(cos(2*(t+h/2.0)*(w+k1/2.0))+sin(3*(t+h/2.0)*(w+k1/2.0)));
        k3=h*(cos(2*(t+h/2.0)*(w+k2/2.0))+sin(3*(t+h/2.0)*(w+k2/2.0)));
        k4=h*(cos(2*(t+h)*(w+k3))+sin(3*(t+h)*(w+k3)));
        w=w+(k1+2.0*k2+2.0*k3+k4)/6.0;
        fprintf(fptr4, "%f\t %f \t \n", t+h, w);
    }
    r4=w;
}

