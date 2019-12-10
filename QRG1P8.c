#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define tags 0
#define tagr 0
//#include <cblas.h>
#include <math.h>
#include <time.h>
#include<sys/time.h>
#define mrows 20
#define ncols 10
#define N 100000

   double norm(double *M, int n);

int main(int argc,char* argv[])
{


printf("size of my matrix is %d and %d", mrows,ncols);

	printf("\n\n\n\n");
 //clock_t t1,t2;
 //struct timeval tv1,tv2; struct timezone tz;
	
	
	double *A, *Qt, *Q, *R, *r, *e, *e1, *z, *a1, *Aa, *AT;

	double no, no1;

	int i, i1, j, j1, k, k1, m, o1, dest, dest1, destw, source, rank=0, numtasks, offset, averow;

	int avecol, cols;

	int iw, jw, kw, ow, id, is;

		int q = 5; 

	A = (double *)malloc( (mrows*ncols) * sizeof (double));

    Qt = (double *)malloc( (mrows*ncols) * sizeof (double));

    Q = (double *)malloc( (mrows*ncols) * sizeof (double));

    R = (double *)malloc( (ncols*ncols) * sizeof (double));

	 r = (double *)malloc( (ncols*ncols) * sizeof (double));

	e = (double *)malloc( (mrows*ncols) * sizeof (double));

	e1 = (double *)malloc( (mrows*ncols) * sizeof (double));

	a1 = (double *)malloc( (mrows) * sizeof (double));
 
    z = (double *)malloc( (mrows) * sizeof (double));


	MPI_Status stat;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
/*initilaizing out matrix, a*/

for(i=0;i<mrows;i++)

		       for(j=0;j<ncols;j++)

			A[i+j*mrows]=i+j;


	A[0] =-1; A[8] = 2;   A[16] = 0;  A[24] = 1;   A[32] = 4;  A[40] = 1;  //A[48]=2;
	A[1] = 0; A[9] = 4;   A[17] = 1;  A[25] = -1;  A[33] = -2; A[41] = 3;  //A[49]=1;
	A[2] =1;  A[10] = -1; A[18] = 2;  A[26] = 4;   A[34]=-3;   A[42]=1;   // A[50]=0;
    A[3] =2;  A[11] = 0;  A[19] = 1;  A[27] = -3;  A[35] = 4;  A[43] = 5;  //A[51]=1;
    A[4] = 1; A[12] = 2;  A[20] = 0;  A[28] = 4;   A[36] = -3; A[44] = 1;  //A[52]=-2;
	A[5] =-2; A[13] = 1;  A[21] = 3;  A[29] = -1;  A[37]=1;    A[45]=2;    //A[53]=5;
	A[6] =2;  A[14] = 0;  A[22] = 4;  A[30] = -1;  A[38]=1;    A[46]=-3;  // A[54]=2;
	A[7] =3;  A[15] = -4; A[23] = 1;  A[31] = 2;   A[39]=4;    A[47]=-1;   //A[55]=2;

		
	printf("\n\n\n");

  //for(i=0;i<mrows;i++){

		//for(j=0;j<ncols;j++){

			//printf(" %lf ", A[i + j*mrows]);}

	//printf("\n");} 

	printf("\n\n\n\n");

	/*Defining appropriate data-type for the ghost regions of each process*/

			MPI_Datatype columntype;

            MPI_Type_vector(mrows, 1,1, MPI_DOUBLE, &columntype);
            MPI_Type_commit(&columntype);

        if( ncols % q==0 )

		avecol = ncols / q;

		else

			avecol = ncols / q + 1;

	if (rank == 0){



		//printf(" no. of tasks are: %d\n", numtasks);


        if( ncols % q==0 )

		avecol = ncols / q;

		else

			avecol = ncols / q + 1;

       

		offset = 0;

		for(j=0;j<q;j++){

			cols = avecol;

		if(j==0)
		{
			j=1;
			
			offset = cols*mrows;
		}

		dest = j;

			if((ncols%q)!=0 && j==q-1)

				cols = ncols % avecol;

		//printf("\n\n\n");

	//	printf("cols is %d", cols);

		//printf("\n\n\n");

			/*appointing neighbooring processors in a array*/

				MPI_Send(&cols,1,MPI_INT,dest,tags,MPI_COMM_WORLD);

				/*sending the necessary elements of the main matrix to each task*/

				for(k=0;k<cols;k++){

					MPI_Send(&offset,1,MPI_INT,dest,tags,MPI_COMM_WORLD);

				    MPI_Send(&A[offset],1,columntype,dest,tags,MPI_COMM_WORLD);

				if(k!=cols-1)

				offset+=mrows;

				}

			offset += mrows;
		}

	}


		//cols = avecol;
		//offset = cols;
		
		//for(i=0;i<mrows;i++)
	
			//a1[i] = A[i+mrows];

    //double no = cblas_dnrm2 ( mrows, a1, 1 );

    


/* workers as the rest of tasks to be sent*/

	if(rank != 0)

	{

	source = 0;
   
				MPI_Recv(&cols, 1, MPI_INT, source, tagr, MPI_COMM_WORLD, &stat);
				//MPI_Recv(&right, 1, MPI_INT, source, tagr, MPI_COMM_WORLD, &stat);
				//MPI_Recv(&left, 1, MPI_INT, source, tagr, MPI_COMM_WORLD, &stat);

for(k=0;k<cols;k++){
	
				MPI_Recv(&offset, 1, MPI_INT, source, tagr,MPI_COMM_WORLD, &stat);

	            //MPI_Recv(&offset1, 1, MPI_INT, source, tagr,MPI_COMM_WORLD, &stat);

	            MPI_Recv(&A[offset],1,columntype,source,tagr,MPI_COMM_WORLD,&stat);

			if(k!=cols){

				offset+=mrows;}

	//printf("\n\n\n");

	//printf("offset of worker is %d", offset);
 
	  //printf("\n\n\n");
			
			}
	


			//for(i1=0;i1<cols;i1++){

		//	MPI_Recv(&e1[i1*mrows],1,columntype,source,tagr,MPI_COMM_WORLD,&stat);}

	//}
	}

//to begin the necessary computation for QR factorization


	if(rank==0){


        if( ncols % q==0 )

		avecol = ncols / q;

		else

			avecol = ncols / q + 1;

		cols = avecol;
		//offset = cols;
		
		for(i=0;i<mrows;i++)
	
			a1[i] = A[i];

    no = norm(a1,mrows);

for(i=0;i<mrows;i++)

		e[i] = A[i]/no;

//r[0] = 0;

	for(i=0;i<mrows;i++)

		r[0] += e[i]*A[i];

	//for(i=0;i<mrows;i++){

		
		//	printf(" %lf ", e[i + mrows]);

		printf("\n\n");


		//printf("\n\n\n");

		//printf("cols is %d", cols);

		//printf("\n\n\n");

		for(j1=1;j1<cols;j1++)

		{
		
		for(i=0;i<mrows;i++)
		{
			z[i] = A[i+mrows*j1];
		}

		   for(o1=0;o1<j1;o1++){

			   for(k1=0;k1<mrows;k1++)
		    {
				r[o1+j1*ncols] += A[k1+j1*mrows]*e[k1+o1*mrows];
			}

              for(k1=0;k1<mrows;k1++)
			{
				z[k1]-= r[o1+j1*ncols]*e[k1+o1*mrows];
			}
		}

	 no1 = norm(z,mrows);

		for(k1=0;k1<mrows;k1++)
		{
			e[k1+j1*mrows] = z[k1]/no1;
		}

		for(k1=0;k1<mrows;k1++)
			 
			r[j1+j1*ncols]+=e[k1+j1*mrows]*A[k1+j1*mrows];
	}

	

			printf("\n\n\n");

	for(i=0;i<mrows;i++){

		for(j=0;j<cols;j++){

			printf(" %lf ", e[i + j*mrows]);}

	printf("\n\n\n");} 


    

		for(dest1=1;dest1<q;dest1++)

		{
			  //MPI_Send(&avecol,1,columntype,dest1,tags,MPI_COMM_WORLD);

			printf("\n\n\n");

			//printf("avecol is %d", cols);

			printf("\n\n\n");

			//printf("dest1 is %d", dest1);

			printf("\n\n\n");

				//MPI_Send(&dest1,1,columntype,dest1,tags,MPI_COMM_WORLD);

			for(i1=0;i1<cols;i1++){

	        MPI_Send(&e[i1*mrows],1,columntype,dest1,tags,MPI_COMM_WORLD);
			
			}
		}

}

	if(rank!=0)

	{

		//printf("\n\n\n");

		//printf("cols is %d", cols);

		//printf("\n\n\n");

		source = 0;

       // MPI_Recv(&avecol,1,columntype,source,tagr,MPI_COMM_WORLD,&stat);

		//MPI_Recv(&dest1,1,columntype,source,tagr,MPI_COMM_WORLD,&stat);

			for(i=0;i<avecol;i++){

		MPI_Recv(&e[i*mrows],1,columntype,source,tagr,MPI_COMM_WORLD,&stat);
		
		}

	}
	
if(rank!=0 && rank!=1){

printf("\n\n\n");

		printf("cols is %d", cols);

		printf("\n\n\n");

		//if(rank!=q-1){

    for(is=rank-1;is>0;is--){
			
//for(is=1;is<q;is++){

			for(i=0;i<avecol;i++){

			 MPI_Recv(&e[is*mrows*avecol+i*mrows],1,columntype,is,tags,MPI_COMM_WORLD, &stat);
				
				}
			}
	//}
	}

// sending the necessary data among processers

		if(rank!=0){

	//printf("\n\n\n");

		//printf("cols is %d", cols);

		//printf("\n\n\n");

		for(jw=0;jw<cols;jw++)

		{
		
		for(iw=0;iw<mrows;iw++)
		{
			z[iw] = A[iw+jw*mrows+mrows*rank*cols];
		}

		   for(ow=0;ow<jw+rank*cols;ow++){	

			   for(kw=0;kw<mrows;kw++)
					
		    {
				r[ow+jw*ncols+ncols*cols*rank] += A[kw+jw*mrows+mrows*rank*cols]*e[kw+ow*mrows];//+pw*mrows];
			}

              for(kw=0;kw<mrows;kw++)
			{
				z[kw]-= r[ow+jw*ncols+ncols*cols*rank]*e[kw+ow*mrows];
			}
		}

	 no1 = norm(z,mrows);

		for(kw=0;kw<mrows;kw++)
		{
			e[kw+(jw+rank*cols)*mrows] = z[kw]/no1;
		}

		for(kw=0;kw<mrows;kw++)
			 
			r[jw+rank*cols+jw*ncols+ncols*cols*rank]+=e[kw+jw*mrows+mrows*rank*cols]*A[kw+jw*mrows+mrows*rank*cols];
	}
	}

	

	//if(rank==2){

	//for(i=0;i<mrows;i++){

		//for(j=0;j<ncols;j++){

		//	printf(" %lf ", e[i + j*mrows]);}

//	printf("\n");} }


        if(rank!=q-1 && rank!=0) {

//printf("\n\n\n");

		//printf("cols is %d", cols);

		//printf("\n\n\n");

       for(id=rank+1;id<q;id++){
			
		for(i=0;i<avecol;i++){

			 MPI_Send(&e[rank*mrows*cols+i*mrows],1,columntype,id,tags,MPI_COMM_WORLD);
				
			}
		
	}
		}



	

		



		


	

      







	



































   // if(rank == 1){
	
	//for(i=0;i<mrows;i++){

		//for(j=0;j<ncols;j++){

			//printf(" %lf ", A[i + j*mrows]);}

	//printf("\n");} }





	MPI_Finalize();

	if(rank==2){

	for(i=0;i<mrows;i++){

		for(j=0;j<ncols;j++){

			printf(" %lf ", e[i + j*mrows]);}

	printf("\n");}} //}}



   


}

double norm(double *v, int n)

{

	int q;
	
	double S=0,norm1;

	for(q=0;q<n;q++){
	
	S+=v[q]*v[q];}

	norm1=sqrt(S);
	
	return norm1;
}

	




	

  











































		
		     


			

		    








	


			


 
