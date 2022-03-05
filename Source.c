#include "bisec_clustering_radius.h"

 //#define KMAX   128
# define R1 600
# define R2 200
//0.115811

int main(int argc, char** argv) {
    
    int i, j,kk, k, count, num_outdata, *outdata = NULL,nclusters;
    char *filename_path = NULL;
    int dim ,ndata,ndata2;
    double *data,tmp,*data2;
    data= (double *) calloc((ndata*dim), sizeof(double));
    // data2= (double *) calloc((ndata2*dim), sizeof(double));
    
  dim = atoi(argv[1]),  ndata = atoi(argv[2]),      ndata2=atoi(argv[3]);
  // dim = atoi(argv[2]),  ndata = atoi(argv[4]);
   

    
    for (i=0; i < argc; ++i) { printf("argv[%d] = '%s'\n", i, argv[i]); }
    
    printf("dim= %i   ndata= %i\n", dim,ndata);

    printf("-----------------------------------------------------------------------------\n");

    tmp = 1.0/RAND_MAX ; /*** Generating data[i] ***/
    for(i=0;i<(ndata*dim);i++)
    {
       data[i] = ((double)rand() / (double)RAND_MAX) * 1000.00;
        //data[i] = (double)rand()*tmp - 0.5;    /*** Generating data[i] in [-0.5, 05] ***/
        
        printf("\n data[%d]=%fl",i,data[i]);
    }

    /**********************************************************************/

  
    printf("call readFile\n");
//    filename_path="/Users/zohrehsafari/Desktop/ColorHistogram.dat";
//    filename_path ="/home/zosafari/profcode2/CASP.csv.bin";
//    filename_path="/home/zosafari/profcode2/ColorHistogram.dat";
//    filename_path="/home/zosafari/profcode2/HIGGS.csv";
//    filename_path ="/tmp/HIGGS.dat";


//    FILE * fp = fopen(filename_path,"rb");
//      if (fp == NULL) { printf("Failed to open binary file: %s.\n",filename_path); return 0; }
//      fread(data, sizeof(double), dim*ndata, fp);
//      fclose(fp);
//      printf("dim= %i   ndata= %i\n", dim,ndata);
    
//       for (i=0;i<ndata*dim;i++)
//            printf("data[%d]=%.10lf\n",i,data[i]);
//
    
//    #endif
     /**********************************************************************/

//  double *  datum  = (double *)calloc(dim, sizeof(double)) ;
//  double * buf  = (double *)calloc(ndata*dim, sizeof(double)) ;
//  int * cluster_assign = (int *)calloc(ndata, sizeof(int)) ;
//  double * cluster_center =(double *)calloc(ndata*dim, sizeof(double));
//  int  cluster_start[ndata], cluster_size[ndata];
//  kk = KMAX ;
//  double* cluster_ssd =(double *)calloc(ndata, sizeof(double));
//
//  double * cluster_radius =(double *)calloc(ndata*dim, sizeof(double));
//  double *cluster_center0 =(double *)calloc(dim, sizeof(double));
//  double *cluster_radius_pt =(double *)calloc(dim, sizeof(double));


//    i0 = 0,im = ndata ;
//    int iterat_limit=5;
//     bkmeans(ndata,5,  R2,  dim,  i0, im, data, // line of input
//               cluster_assign, datum,                                       // line of buffers
//                cluster_center, cluster_radius,
//                cluster_start, cluster_size, cluster_ssd);
//


    //Now begin building the structure.
    /*zz***************************************************************************************************zz*/
    
    
    struct GrowingTree3L *SS = (struct GrowingTree3L *)malloc(sizeof(struct GrowingTree3L));
    printf("Start growingtree_construc...\n");
    
    GrowingTree_construc(dim,R1,R2,ndata,data,SS);
    
    printf("\n \n call next files\n");
//   filename_path="/Users/zohrehsafari/Desktop/CASP.csv.bin";
//    FILE * fp1 = fopen(filename_path,"rb");
//    if (fp1 == NULL)  {printf("Failed to open binary file: %s.\n",filename_path); return 0;}
//    fread(data2, sizeof(double), dim*ndata2, fp1);
//    fclose(fp1);
//    printf(" dim= %i   ndata2= %i\n", dim,ndata2);
//
//       // for (i=0;i<ndata2*dim;i++)   printf("data2[%d]=%.10lf\n",i,data2[i]);

    
    printf("-----------------------------------------------------------------------------\n");
//    free(data);
//    ndata=ndata2;
//    data= (double *) calloc((ndata*dim), sizeof(double));
//
//    tmp = 1.0/RAND_MAX ; /*** Generating next data[i] ***/
//    for(i=0;i<(ndata*dim);i++)
//    {
//        data[i] = ((double)rand() / (double)RAND_MAX) * 1000.00;
//
//        printf("\n data2[%d]=%fl",i,data[i]);
//    }
//
    /**********************************************************************/
    //for further files:
    
    //for (i = 0; i < ndata2 * dim; i++)
    
    
//    GrowingTree_grow (dim,R1,R2,ndata,data,SS);
    // I delete KMAX as a passing parameter in the function, so I thinl it is better we have it and the use realloc. If we donot need KMAX, we can remove defining KMAX as a constant .
    
  
    
    
    
    
    /**********************************************************************/
    printf("\n\n \n Start searches...\n\n");
    
   
    
    
    
    free(data) ;
    return 0;
    
    
}

