/****** File:  SphereTree4L. c ******/
#ifndef SPHT
#define SPHT
#include "bisec_clustering_radius.h"



int GrowingTree_construc(int dim, double R1,double R2,int ndata, double *data, struct GrowingTree3L *ptr_root) {
    
    printf ("\n hello GrowingTree");

   int  i0,im,i,j,k,child_start,child_end,num_shift,startnewcluster,
   cluster_start[ndata],cluster_size[ndata],
   nclusters,node_indx,N1,N2,NC,k_split;
   int *cluster_assign;
   double tmp, dist_ss, *datum, *buf,*cluster_center, *current_center, *child_center ;
   int flagfirstchecking=0;
    
   datum  = (double *)calloc(dim, sizeof(double)) ;
   buf = (double *)calloc(ndata*dim, sizeof(double)) ;
   cluster_assign = (int *)calloc(ndata, sizeof(int)) ;
   current_center= (double *)calloc(dim, sizeof(double)) ;
   child_center  = (double *)calloc(dim, sizeof(double)) ;
   cluster_center =(double *)calloc(ndata*dim, sizeof(double));
   double * cluster_radius =(double *)calloc(ndata, sizeof(double));
    double * cluster_mean = (double *) calloc(dim, sizeof(double));
  
   ptr_root->dim = dim;
   ptr_root->NumData = ndata;
//int leafmax =ptr_root->leafmax; 

   i0=0; im=ndata;

/**************finding the radius of whole data set for threshold and easily pass to functions*****************/






/*** Compute centroid of the whole dataset ***/
    memset(cluster_center, 0.0, dim*sizeof(double)) ;
    for(i=i0; i<im; i++)
        for(j=0; j<dim; j++) cluster_center[j] += data[i*dim+j];
    for(j=0; j<dim; j++) cluster_center[j] /= (double)(im-i0) ;
    
     /*** Compute farthest pt to centroid. Store it in radius_pt[0] ***/
    memset(cluster_radius, 0.0, ndata*sizeof(double)) ;
   // ssd_initial = 0.0 ;
    for(i=i0; i<im; i++) {
        tmp = calc_dist_square(dim, data+i*dim, cluster_center) ;
       // ssd_initial += tmp ;
        if(cluster_radius[0] < tmp) {
            cluster_radius[0] = tmp ;
            //memcpy(cluster_radius_pt, data+i*dim, dim*sizeof(double)) ;
        }
}
              R1=(cluster_radius[0]/4);
//R1=log10(cluster_radius[0])/log10(2);              

printf("\nJust to make sure the values are correct:\n in Growing cluster_radius[0] is %fl and  R1 is %fl\n",cluster_radius[0], R1);                                                          




/*************************************************************************************************************/

   N1 = bkmeans(flagfirstchecking,ndata, 5, R1, dim, i0, im, data,cluster_assign,datum,cluster_center, cluster_radius, cluster_start, cluster_size);
   printf("\n now N1 is calculating....N1 is= %d\n",N1);

   ptr_root->N1 = N1 ;
   ptr_root->L1node_StartChild = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1node_NumChldrn  = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1node_StartDatum = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1node_DataSize   = (int *) calloc(N1, sizeof(int)) ;
   ptr_root->L1CentersRadii = (double *) calloc(N1*(dim+1), sizeof(double)) ;

   for(k=0; k<(ptr_root->N1); k++) {
      ptr_root->L1node_StartDatum[k] = cluster_start[k] ;
      ptr_root->L1node_DataSize[k] = cluster_size[k]  ;
      for(j=0; j<dim; j++) ptr_root->L1CentersRadii[k*(dim+1)+j] = cluster_center[k*dim+j];
      ptr_root->L1CentersRadii[k*(dim+1)+dim] = cluster_radius[k] ;
       ptr_root->L1node_NumChldrn[k]  = cluster_size[k];
      // ptr_root->L1node_StartChild[k] = cluster_start[k] ;   
      printf("GrWT L1 done  N1=%d, radius=%f, size=%d, start=%d\n", ptr_root->N1, cluster_radius[k], cluster_size[k], cluster_start[k]);

   }
    N1=ptr_root->N1;
    printf ("\nfinally N1 is= %d\n",N1);
printf("\n LET'S Go TO CHECK THE RADIUS AFTER KMEAN RUNNING\n");
 /*********************************************************************/


for (k=0;k<(ptr_root->N1);k++)
{
    if (ptr_root->L1CentersRadii[k*(dim+1)+dim] > R1 )

    { printf("cluster %d should be ready to be split as it has a radius greater than R1  ",k);
        // create a function for following statements:
       // flagfirstchecking=1; /*add it to the function SplitR_constru and bkmean*/
       flagfirstchecking=1;
         k_split=k;
        i0=ptr_root->L1node_StartDatum[k];
        im=i0+ptr_root->L1node_DataSize[k];
        startnewcluster=i0;
        
        
        NC = bkmeans(flagfirstchecking,ndata, 5, R1, dim, i0, im, data,cluster_assign,datum,cluster_center, cluster_radius, cluster_start, cluster_size);
        printf("Now cluster number %d converted to %d clusters \n");
        
         printf("\n ********************TO FIND SPACE FOR # OF NUM_SHIF**********\n");
        
        num_shift=NC-1;
        for(i= N1-1; i>=k_split+1 ; i--)
            { ptr_root->L1node_StartDatum[i+num_shift] = ptr_root->L1node_StartDatum[i] ;
            ptr_root->L1node_DataSize[i+num_shift]  = ptr_root->L1node_DataSize[i] ;
            ptr_root->L1CentersRadii[((i+num_shift)*(dim+1))+dim]= ptr_root->L1CentersRadii[(i*(dim+1))+dim] ;
            memcpy(ptr_root->L1CentersRadii + (i+num_shift)*dim, ptr_root->L1CentersRadii + i*dim, dim*sizeof(double)) ;
            ptr_root->L1node_NumChldrn[i+num_shift]=ptr_root->L1node_NumChldrn[i];
            }
        
        
        
    
        for(i=0; i<NC; i++)
              { ptr_root->L1node_NumChldrn[k_split]  = cluster_size[i] ;    printf("\n L1node_NumChldrn[%d] is %d\n",k_split, ptr_root->L1node_NumChldrn[k_split]);
                ptr_root->L1node_StartDatum[k_split] =startnewcluster+cluster_start[i] ;

                 ptr_root->L1node_DataSize[k_split]  = cluster_size[i] ;
                 printf("\n L1node_DataSize[%d] is %d\n",k_split, ptr_root->L1node_DataSize[k_split]);

                //if (ptr_root->L1node_DataSize[k_split] > leafmax) { leafmax = ptr_root->L1node_DataSize[k_split]; }
                for(j=0;j<dim;j++) ptr_root->L1CentersRadii[k_split*(dim+1)+j] =cluster_center[i*dim+j];
                ptr_root->L1CentersRadii[k_split*(dim+1)+dim] = cluster_radius[i] ;
                k_split++ ;}
            
           ptr_root->N1+=NC;
        
printf ("\n***** after changing cluster %d,  N1 is= %d\n",k,N1);
     
        
    }
    printf("YES : we are done , the Tree is going to be printed out\n");
    printf ("****************** Finally N1 with suitable R1 is= %d\n",N1);
    
    for (k=0;k<(ptr_root->N1);k++)
        printf("New GrWT L1 done  N1=%d, radius=%f, size=%d, start=%d\n", ptr_root->N1, cluster_radius[k], cluster_size[k], cluster_start[k]);
        
        
        
        
        
        
        
        
        
        
        
        

        /* Compute Mean for Split Location */
//        memset(cluster_mean, 0.0, dim*sizeof(double)) ;
//        for(k=i0; k<im; k++)
//            for(j=0; j<dim; j++) { /*printf(" \n data[i*dim+j]= %fl",data[i*dim+j]);*/  cluster_mean[j] += data[k*dim+j];}
//        for(j=0; j<dim; j++) cluster_mean[j] /= (double)(im-i0) ;
//        printf ("\n the mean while calculate outside split is");
//        for (j=0;j<dim;j++)  printf(" \n cluster_mean[%d]= %fl",j,cluster_mean[j]);
//
//
//        printf("\n ************** we are here to call split\n");
//
//        SplitR_construc(R1,i0,im, k_split,data,ndata,dim,cluster_assign,ptr_root,flagfirstchecking,cluster_mean);
    

}



#if 0
   FILE *L1_radius;
   FILE *L1_size;
   FILE *L1_r_s;
   L1_radius = fopen("SST_L1_radius.xvg", "w");
   L1_size   = fopen("SST_L1_size.xvg", "w");
   L1_r_s    = fopen("SST_L1_r_s.xvg", "w");

   for(k=0; k<ptr_root->N1; k++) {
      fprintf(L1_radius, "%d   %f\n", k, ptr_root->L1CentersRadii[k*(dim+1)+dim] );
      fprintf(L1_size, "%d    %d\n", k, ptr_root->L1node_DataSize[k]);
      fprintf(L1_r_s, "%d    %f\n", k, ptr_root->L1CentersRadii[k*(DIM+1)+DIM]/ptr_root->L1node_DataSize[k]);
   }

   fclose(L1_radius);
   fclose(L1_size);
   fclose(L1_r_s);


//   FILE *L2_radius;
//   FILE *L2_size;
//   FILE *L2_r_s;
//   L2_radius = fopen("LP_L2_radius.xvg", "w");
//   L2_size   = fopen("LP_L2_size.xvg", "w");
//   L2_r_s    = fopen("LP_L2_r_s.xvg", "w");

//   for(k=0; k<ptr_root->N2; k++) {
//      fprintf(L2_radius, "%d   %f\n", k, ptr_root->L2CentersRadii[k*(DIM+1)+DIM] );                         // correct it
//      fprintf(L2_size, "%d    %d\n", k, ptr_root->L2node_DataSize[k]);
//      fprintf(L2_r_s, "%d    %f\n", k, ptr_root->L2CentersRadii[k*(DIM+1)+DIM]/ptr_root->L2node_DataSize[k]);
//   }

//   fclose(L2_radius);
//   fclose(L2_size);
//   fclose(L2_r_

# endif

   /*** free memory: buf, cluster_assign ***/
   free(datum);
   free(buf) ;
   free(cluster_assign);
   free(cluster_center);
   free(current_center);
   free(child_center) ;
  free(cluster_radius);
free(cluster_mean);
return 0;
   /*** free_memory: ptr_root->ChldrnCentersRadii, ptr_root->ChildNodeArrays  in main() ***/

} /********************************** End of function SStree_construc() ****************************/


int SplitR_construc(double R1,int i0,int im, int k_split, double *data,int ndata, int dim, int * cluster_assign, struct GrowingTree3L *ptr_root,int flagfirstchecking,double * cluster_mean)
{        

    printf("\n we are in split function ");

    int k,num_shift, i, j, maxVardim, size[2],start[2], start0, start1, end0, end1, N = 0;
    double *datum, *center, tmp = 0.0, radius[2], *variance, maxVar = 0.0, maxMean = 0.0;

    datum = (double *) calloc(dim, sizeof(double));
    variance = (double *) calloc(dim, sizeof(double));
    center = (double *) calloc(2 * dim, sizeof(double));


    // NC = bkmeans(ndata, 5, R1, dim, i0, im, data,cluster_assign,datum,cluster_center, cluster_radius, cluster_start, cluster_size);
    // I should apply two-means, only split the cluster into 2 clusters, ...yes it would be better.


    /*Compute Mean for Split Location*/
    memset(cluster_mean, 0.0, dim * sizeof(double));
    for (i = i0; i < im; i++)
        for (j = 0; j < dim; j++) { cluster_mean[j] += data[i * dim + j]; }
    for (j = 0; j < dim; j++) cluster_mean[j] /= (double) (im - i0);
    printf("\n the mean while calculate inside split is");
    for (j = 0; j < dim; j++) printf(" \n cluster_mean[%d]= %fl", j, cluster_mean[j]);



    /***find dimension with the highest variance***/
    for (i = i0; i < im; i++) {
        for (j = 0; j < dim; j++) {
            tmp = data[i * dim + j] - cluster_mean[j];
            variance[j] += tmp * tmp;
        }
    }
    for (j = 0; j < dim; j++) {
        variance[j] /= (double) (im - i0);
        if (variance[j] > maxVar) { maxVardim = j; /*maxVar=variance[j];*/}
    }
    maxMean = cluster_mean[maxVardim];
    printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ maxVardim= %d and maxMean=%fl  (maxMean=cluster_mean[maxVardim]; )@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n",
           maxVardim, maxMean);



    /***** clustering into 2parts   *****/
    printf("\n clustering into 2parts  \n  ");
    size[0] = 0;
    size[1] = 0;
    for (i = i0; i < im; i++)   /* to assign these data and find center and size of each part*/
    {
        memcpy(datum, data + i * dim, dim * sizeof(double));
        // printf("\n LETS SEE DATA POINT ::::   data[i*dim+maxVardim]=data[%d*%d+%d]=%fl\n",i,dim,maxVardim,data[i*dim+maxVardim]);
        k = data[i * dim + maxVardim] < maxMean ? 0 : 1; // printf("\n K is =%d\n",k);
        // if (data[i*dim+maxVardim] < maxMean)  k=0 ;  else k=1;    printf("\n K is =%d\n",k);

        cluster_assign[i] = k;  //printf("\n cluster_assign[%d]=%d\n",i,k);
        size[k]++;
        for (j = 0; j < dim; j++) center[k * dim + j] += datum[j];
    }


//    printf("\n ::: ::: ::: :::data array ");  for(k=0;k<ndata*dim;k++)   printf("\n data[%d]=%fl",k,data[k]);

/***update caenter of split cluster***/
    if (size[0] > 0)
        for (j = 0; j < dim; j++) center[j] /= size[0];
    printf("\n size[0]=%d\n ", size[0]);
    if (size[1] > 0)
        for (j = 0; j < dim; j++) center[dim + j] /= size[1];
    printf("\n size[1]=%d\n ", size[1]);

    /****** Data Re-ordering ******/
    start[0] = i0;
    start[1] = i0 + size[0];
    start0 = start[0];
    end0 = start0 + size[0];
    start1 = start[1];
    end1 = start1 + size[1];
    radius[0] = 0.0;
    radius[1] = 0.0;

    while (start0 < end0) { /*** end0 and end1 MUST NOT EXCEED im ***/
        while ((start0 < end0) && (cluster_assign[start0] == 0)) start0++; //printf("\nstart0=%d\n",start0);
        while ((start1 < end1) && (cluster_assign[start1] == 1)) start1++; //printf("\nstart1=%d\n",start1);
        if ((start0 < end0) && (start1 < end1)) {
            memcpy(datum, data + start0 * dim, dim * sizeof(double));
            memcpy(data + start0 * dim, data + start1 * dim, dim * sizeof(double));
            memcpy(data + start1 * dim, datum, dim * sizeof(double));
            start0++;
            start1++;
        } else if ((start0 < end0) && (start1 == end1)) {
            printf("Error: start0<end0 && start1==end1\n");    /* I have to  consider the update of this condition here as it may happen without return,because it goes out not only split but also insertion function*/ /*return ++ndata;*/     printf(
                    "\nstart0=%d and start1=%d and end0=%d and end1=%d\n", start0, start1, end0, end1);
        }

        else if ((start0 == end0) && (start1 < end1)) {
            printf("Error: start0==end0 && start1<end1\n");  /* I have to  consider the update of this condition here as it may happen*/ /*return ++ndata;*/ printf(
                    "\nstart0=%d and start1=%d and end0=%d and end1=%d\n", start0, start1, end0, end1);
        }
        else { ; }
    }  /*** End of Data Re-ordering ***/

/****update cluster_assign with the real one****/ //
    i0 = start[0];
    end0 = start0 + size[0];
    for (i = i0; i < end0; i++) { cluster_assign[i] = k_split; }
    /*PLEASE CHECK if my guess is true that cluster_assign is only based on indexes of level 1-clusters*/
    i0 = start[1];
    end1 = start1 + size[1];
     for (i = i0; i < end1; i++) { cluster_assign[i] = k_split + 1; }



/**update radius of split node**/
    printf("\n update radius of cluster \n ");
    for (k = 0; k < 2; k++) {
        radius[k] = 0.0;
        for (j = start[k]; j < start[k] + size[k]; j++) {
            tmp = calc_dist_square(dim, data + j * dim, center + k * dim);
            if (radius[k] < tmp)   radius[k] = tmp;
        }
        printf("the radius for k=%d is %fl\n",k,radius[k]);
    }




/**shift and save new cluster**/
    printf("\n save new cluster\n ");
    for (k = ptr_root->N1-1; k > k_split; k--) {
        ptr_root->L1node_StartChild[k + 1] = ptr_root->L1node_StartChild[k];
        ptr_root->L1node_StartDatum[k + 1] = ptr_root->L1node_StartDatum[k];
        ptr_root->L1node_DataSize[k + 1] = ptr_root->L1node_DataSize[k];
        ptr_root->L1node_NumChldrn[k + 1] = ptr_root->L1node_NumChldrn[k];
        ptr_root->L1CentersRadii[((k + 1) * (dim + 1)) + dim] = ptr_root->L1CentersRadii[(k * (dim + 1)) + dim];
        memcpy(ptr_root->L1CentersRadii + (k + 1) * dim, ptr_root->L1CentersRadii + k * dim, dim * sizeof(double));

    }


    k = k_split;
    ptr_root->L1node_StartDatum[k] = start[0];
    ptr_root->L1node_StartDatum[k + 1] = start[1];

    ptr_root->L1node_DataSize[k] = size[0];
    ptr_root->L1node_DataSize[k + 1] = size[1];
    ptr_root->L1node_NumChldrn[k] = size[0];
    ptr_root->L1node_NumChldrn[k + 1] = size[1];
    ptr_root->L1node_StartChild[k] = start[0];
    ptr_root->L1node_StartChild[k + 1] = start[1];


    ptr_root->L1CentersRadii[(k * (dim + 1)) + dim] = radius[0];
    ptr_root->L1CentersRadii[((k + 1) * (dim + 1)) + dim] = radius[1];
    memcpy(ptr_root->L1CentersRadii + k * dim, center, dim * sizeof(double));
    memcpy(ptr_root->L1CentersRadii + (k + 1) * dim, center + dim, dim * sizeof(double));

    ptr_root->N1++;


//check the radius of this new 2clusters:
for (i=k_split ; i <= k_split+1 ; i++)
 {
     if (radius[i] > R1)
   {
        printf("cluster %d needs to be split\n ",i);
        i0=ptr_root->L1node_StartDatum[i];
        im=i0+ptr_root->L1node_DataSize[i];
        k_split=i;
        SplitR_construc(R1,i0,im, k_split,data,ndata,dim,cluster_assign,ptr_root,flagfirstchecking,cluster_mean);

    }

  }

return 0;
}























/*********************************************************************************************************/

int GrowingTree_grow (int dim,double R1,double R2,int ndata,double *data, struct GrowingTree3L *ptr_root)
{
   double *dist, *center, *datum,mindist,mindist2,tmp=0.0;
       int i,k,j,index_close1_cluster,index_close2_cluster, Index_Outer=0,Index_Inner=0;
       center= (double *) calloc(dim, sizeof(double)) ;
       datum= (double *) calloc(dim, sizeof(double)) ;
       dist= (double *) calloc(dim, sizeof(double)) ;                             // the distances of 2 points which have the same dimmension, so I think the size should be dim   ???
       double * StoreBOuterArray = (double *)malloc(sizeof(double)*ndata*dim);    //the size should be changed  use realloc  ....> .... bezarim baad realloc konim
       double * StoreBInnerArray = (double *)malloc(sizeof(double)*ndata*dim);    //the size should be changed
       int * AssignInnerBigCluster = (int *)malloc(sizeof(int)*ndata);            //...size
       int * BigClusterAdd=(int *)calloc(ptr_root->N1, sizeof(int)) ;
       int * BigClusterOuterAdd=(int *)calloc(ptr_root->N1, sizeof(int)) ;
      int * ChangeSizeSmallCluster[ptr_root->N1];
       for (i=0;i<ptr_root->N1;i++)
     ChangeSizeSmallCluster[i]=(int *)calloc( ptr_root->N2, sizeof(int)) ;   //??
   // int * ChangeSizeSmallCluster[i]=(int *)calloc(ptr_root->N1* ptr_root->N2, sizeof(int)) ;
       int n,t=0;
       int StartInsideBigC[ptr_root->N1];     // keeping start of points which their dist>R2 for each cluster
       memset(StartInsideBigC, 0, ptr_root->N1*sizeof(int)) ;
   
/******* compute distances of all data to every clusters*******/
   
for(i=0; i<ndata ;i++)
{
   for(k=0; k<ptr_root->N1; k++)
          {
           for(j=0; j<dim; j++) center[j] = ptr_root->L1CentersRadii[k*(dim+1)+j] ;
           dist[k] = calc_dist_square(dim, data+i*dim, center);
          }
    mindist=dist[0];
      for (k=1;k<ptr_root->N1;k++)
      {
          if (dist[k]<mindist) {mindist=dist[k];}
         index_close1_cluster=k;
      }
   /***** distinguish if it will be in level 1 or 2 , it is also possible to make NEW cluaters in level 1 or 2********/
    if (mindist > R1)
                   {
                   printf("\n mindist > R1 so we should create new cluster but first store in the array\n");
                       //if (Index_Outer>=....)   StoreBOuterArray = (double *)realloc(sizeof(double)*....+1*dim);                               
                        for(j=0; j<dim; j++) StoreBOuterArray[j+(Index_Outer*dim)]=data[(i*dim)+j];
                       Index_Outer+=1;
                   }
   
    else if (mindist <= R1)
                  {
   /**calculate distances to inside clusters**/
                       for(k=0;k<ptr_root->L1node_NumChldrn[index_close1_cluster];k++)
                            {
                               for(j=0; j<dim; j++) center[j] = ptr_root->L2CentersRadii[index_close1_cluster][k*(dim+1)+j] ;
                               dist[k] = calc_dist_square(dim, data+i*dim, center);
                            }
                       mindist2=dist[0];
                       for (k=1;k<ptr_root->L1node_NumChldrn[index_close1_cluster];k++)
                            {
                                if (dist[k]<mindist2) { mindist2=dist[k];}
                            index_close2_cluster=k;
                            }
            if (mindist2 > R2)
            /****** we should create another inner cluster and also to figure out that cluster,so we need another array to assign****/
                            {
                              //   if (Index_Inner>=.....)   StoreBInnerArray = (double *)realloc(sizeof(double)*.....+1*dim);                       
                               memcpy(datum, data+i*dim, dim*sizeof(double));
                               for(j=0; j<dim; j++)    StoreBInnerArray[j+(Index_Inner*dim)]=data[(i*dim)+j];
                               AssignInnerBigCluster[Index_Inner]=index_close1_cluster;
                               Index_Inner++;
                               BigClusterAdd[index_close1_cluster]++;
                               BigClusterOuterAdd[index_close1_cluster]++ ;
                               for(j=0;j<dim;j++)  ptr_root->L1CentersRadii[index_close1_cluster*dim+j]+= datum[j];
                            }
           else if (mindist2 <= R2)
                            {
                               memcpy(datum, data+i*dim, dim*sizeof(double));
                               ptr_root->L2node_DataSize[index_close1_cluster][index_close2_cluster]+=1;
                               ChangeSizeSmallCluster[index_close1_cluster][index_close2_cluster]++;
                               for(j=0;j<dim;j++)  ptr_root->L2CentersRadii[index_close1_cluster][index_close2_cluster*dim+j]+= datum[j];
                               BigClusterAdd[index_close1_cluster]++;
                             }
                   }
}
   /*******re-adjust size, start of small clusters AND size,start,startchild & numchild  of bigcluster and also reajdust center of small cluaters and big ones.********/
   
 int position=0;
 int SumChangeSizeSmallClusters=0;
 for(i=0;i<ptr_root->N1;i++)
 {
   for(k=0; k<ptr_root->L1node_NumChldrn[i]; k++)
     {
         ptr_root->L2node_DataSize[i][k] += ChangeSizeSmallCluster[i][k];
         SumChangeSizeSmallClusters+= ChangeSizeSmallCluster[i][k];
         if ((i!=0)&&(k==0)){ ptr_root->L2node_StartDatum[i][k] = position;}
             else { for(j=k; k>0; j--)   ptr_root->L2node_StartDatum[i][k] +=position+ChangeSizeSmallCluster[i][k-1];}
       
         if (ptr_root->L2node_DataSize[i][k]>0)   for (j=0;j<dim;j++)  ptr_root->L2CentersRadii[i][k*dim+j]/= (double)ptr_root->L2node_DataSize[i][k];        // re-adjust center of big clusters
     }
     
     StartInsideBigC[i]=ptr_root->L1node_StartDatum[i]+ptr_root->L1node_DataSize[i]+SumChangeSizeSmallClusters;       //later I'll need the start for points which their dists>R2
     ptr_root->L1node_StartDatum[i+1]= StartInsideBigC[i]+BigClusterOuterAdd[i];
     ptr_root->L1node_DataSize[i]+=BigClusterAdd[i];
     position=ptr_root->L1node_StartDatum[i]+ptr_root->L1node_DataSize[i];
     ptr_root->L1node_StartChild[i+1]+=BigClusterAdd[i];
     if (ptr_root->L1node_DataSize[i]>0)   for (j=0;j<dim;j++)  ptr_root->L1CentersRadii[i*dim+j]/= (double)ptr_root->L1node_DataSize[i];               // re-adjust center of big clusters
     SumChangeSizeSmallClusters=0;
     
 }


/***Re-adjust radius of small clusters***/
  for(i=0;i<ptr_root->N1;i++)
   {
     for(k=0; k<ptr_root->N2; k++)
     {
      if (ChangeSizeSmallCluster[i][k]!=0)
      {
          for(j=0; j<dim; j++) center[j] = ptr_root->L2CentersRadii[i][k*(dim+1)+j] ;
             for(j=ptr_root->L2node_StartDatum[i][k]; j<((ptr_root->L2node_DataSize[i][k])-1) ; j++)
               {
                 tmp = calc_dist_square(dim, data+j*dim, center) ;        //data+j*dim  is not true
                 if(( ptr_root->L2CentersRadii[i][k]< tmp) && (tmp<=R2))            /*  ????     radius should be kept as a value         */
                      ptr_root->L2CentersRadii[i][k]  = tmp ;
                }
      }
     }
   }
   
clustering_innerpoints(R2,StoreBInnerArray,AssignInnerBigCluster,dim,ptr_root,Index_Inner,BigClusterOuterAdd,StartInsideBigC);
   
clustering_outerpoints(R1,R2,StoreBOuterArray,AssignInnerBigCluster,dim,ptr_root,ndata,Index_Outer);

return 0;}
/*************************************clustering innerpoints****************************/

/**** first investigating insidepoints  ****/
int clustering_innerpoints(double R2,double *StoreBInnerArray,int *AssignInnerBigCluster,int dim,struct GrowingTree3L *ptr_root,int Index_Inner,int *BigClusterOuterAdd,int *StartInsideBigC)
{
   int  N2,n,position,t ,k,count_ch,j,i0,im,i,N2_1,shift_index=0;
   
   int *cluster_start, *cluster_size,child_end,child_start;
   
   int *cluster_assign;
   double tmp, dist_ss, *datum, *buf,*cluster_center,*current_center, *child_center, *cluster_radius ; //child center va current center taarif kon
   
   datum  = (double *)calloc(dim, sizeof(double)) ;
   buf  = (double *)calloc(Index_Inner*dim, sizeof(double)) ;
   cluster_assign = (int *)calloc(Index_Inner, sizeof(int)) ;
   cluster_center =(double *)calloc(Index_Inner*dim, sizeof(double));   // index_inner can also be very large,so..                           ??
   cluster_size =(int *)calloc(Index_Inner, sizeof(int));
   cluster_start =(int *)calloc(Index_Inner, sizeof(int));
   cluster_radius = (double *)calloc(Index_Inner, sizeof(double)) ;
   
   n=0;
   int start0,end0,start1,start[ptr_root->N1];           // a temporary kepping start
   memset(start,0,ptr_root->N1*sizeof(int));
   
   
   
   /***** order the inside points    *****/
  position=0;
  for(k=0; k<ptr_root->N1; k++)
   {
    if (BigClusterOuterAdd[k]!=0)  start[k]=position;             //if it is Null,what I should hold for that cluster                        ???
       
        position+=BigClusterOuterAdd[k];
   }
   
  for(k=0; k<((ptr_root->N1)-1); k++)
   {
     if (BigClusterOuterAdd[k]!=0)
     {
         start0= start[k] ;
         end0  = start0 + BigClusterOuterAdd[k];
         if (BigClusterOuterAdd[k+1]!=0)  start1= start[k+1];
         else {  while (BigClusterOuterAdd[k+1]==0) k++;   start1= start[k];}
         while(start0 < end0)
         {
             while((start0<end0)&&(AssignInnerBigCluster[start0]==k)) start0++ ;
             if(start0<end0)
             {
                 while((start1<Index_Inner)&&(AssignInnerBigCluster[start1]!=k)) start1++ ;
                 if(start1==Index_Inner) { printf("\nError: start1==index_Inner.\n");   return 0; }
                 memcpy(datum, StoreBInnerArray+start0*dim, dim*sizeof(double)) ;
                 memcpy(StoreBInnerArray+start0*dim, StoreBInnerArray+start1*dim, dim*sizeof(double)) ;
                 memcpy(StoreBInnerArray+start1*dim, datum, dim*sizeof(double)) ;
                 AssignInnerBigCluster[start1] = AssignInnerBigCluster[start0] ;
                 AssignInnerBigCluster[start0] = k;
                 start0++;   start1++;
             }
         }
         
     }
  }/****** End of loop for(k=0; ******/
         
  for(k=0; k<ptr_root->N1; k++)
  {     n=BigClusterOuterAdd[k];
         if (n!=0)
       {
          i0= start[k];
          im= i0 + BigClusterOuterAdd[k] ;
          N2_1 = bkmeans (n,5, R2,dim, i0, im,StoreBInnerArray,cluster_assign, datum,cluster_center, cluster_radius, cluster_start, cluster_size);
           t=ptr_root->L1node_NumChldrn[k]+N2_1;
           ptr_root->L2node_StartDatum[k]= (int *) realloc(ptr_root->L2node_StartDatum[k], t*sizeof(int)) ;
           ptr_root->L2node_DataSize[k]  = (int *) realloc(ptr_root->L2node_DataSize[k], t*sizeof(int)) ;
           ptr_root->L2CentersRadii[k]= (double *) realloc( ptr_root->L2CentersRadii[k],t*(dim+1)*sizeof(double)) ;
           
           count_ch= ptr_root->L1node_NumChldrn[k];
           ptr_root->L1node_NumChldrn[k]+= N2_1;
           ptr_root->N2 += N2_1;
      
          for(i=1; i<N2_1; i++)
           {
              ptr_root->L2node_StartDatum[k][count_ch] = cluster_start[i]+ StartInsideBigC[k];       //???  RECHECK PLEASEEEEEE
              ptr_root->L2node_DataSize[k][count_ch] = cluster_size[i];
              for(j=0;j<dim;j++) ptr_root->L2CentersRadii[k][count_ch*(dim+1)+j] =  cluster_center[i*dim+j];
              ptr_root->L2CentersRadii[k][count_ch*(dim+1)+dim] = cluster_radius[i] ;
              count_ch++ ;
           }
               
        }
       
     }
   N2= ptr_root->N2;
 
   return 0;  //return N2;
}
   
   
/*************************************clustering outerpoints****************************/

/****  investigating outsidepoints  ****/
   
int clustering_outerpoints(int R1,int R2,double *StoreBOuterArray,int *AssignInnerBigCluster ,int dim,struct GrowingTree3L *ptr_root,int ndata,int Index_Outer)
{
   int  node_indx,count_ch,l,j,i0,im,i,k,N1_1,N1_2,shift_index=0;
   int *cluster_start, *cluster_size,child_end,child_start;
   int *cluster_assign;
   double tmp, dist_ss, *datum, *buf,*cluster_center, cluster_radius[Index_Outer],*current_center, *child_center ;
   int n,N1,N2,StartOutsideBigC;
   
   datum  = (double *)calloc(dim, sizeof(double)) ;
   buf  = (double *)calloc(Index_Outer*dim, sizeof(double)) ;
   cluster_assign = (int *)calloc(Index_Outer, sizeof(int)) ;
   cluster_center =(double *)calloc(Index_Outer*dim, sizeof(double));
   cluster_size =(int *)calloc(Index_Outer, sizeof(int));
   cluster_start =(int *)calloc(Index_Outer, sizeof(int));
   
   
   count_ch=0;
   i0= StoreBOuterArray[0];
   im=i0+Index_Outer-1;
   StartOutsideBigC= ptr_root->L1node_StartDatum[(ptr_root->N1)-1]+ptr_root->L1node_DataSize[(ptr_root->N1)-1];
   
   N1_1=bkmeans(Index_Outer,5, R1,dim, i0, im, StoreBOuterArray,cluster_assign,datum,cluster_center, cluster_radius, cluster_start, cluster_size);
   
   count_ch=ptr_root->N1;
   
   for(k=0; k<(N1_1); k++)
      {
        ptr_root->L1node_DataSize[count_ch] = cluster_size[k];
        ptr_root->L1node_StartDatum[count_ch] = cluster_start[k]+StartOutsideBigC ;
        for(j=0; j<dim; j++) ptr_root->L1CentersRadii[count_ch*(dim+1)+j] = cluster_center[k*dim+j];
        ptr_root->L1CentersRadii[count_ch*(dim+1)+dim] = cluster_radius[k] ;
        count_ch++;
      }
   node_indx=0;
   for(k=0; k<N1_1; k++)
     {
       i0= ptr_root->L1node_StartDatum[k];   im= i0 + ptr_root->L1node_DataSize[k] ;
         n=ptr_root->L1node_DataSize[k];                                            //RECHECK PLEASE
       N1_2 = bkmeans(n,5, R1,dim, i0, im, StoreBOuterArray,cluster_assign,datum,cluster_center, cluster_radius, cluster_start, cluster_size) ;

              ptr_root->L1node_NumChldrn[k]  = N1_2;
              ptr_root->L1node_StartChild[k] = StartOutsideBigC +node_indx ;

               for(i=0; i<N1_2; i++)
                 {
                 ptr_root->L2node_StartDatum[k+N1][i] = cluster_start[i] + StartOutsideBigC ;
                 ptr_root->L2node_DataSize[k+N1][i]   = cluster_size[i] ;
          
                 for(j=0;j<dim;j++) ptr_root->L2CentersRadii[k][i*(dim+1)+j] = cluster_center[i*dim+j];
            
                ptr_root->L2CentersRadii[k][i*(dim+1)+dim] = cluster_radius[i] ;
                node_indx++ ;
                 }
    }
   ptr_root->N1+=N1_1;
   N1=ptr_root->N1;
   ptr_root->N2+=node_indx;
   N2=ptr_root->N2;
   


/*** Re-adjust the radius of each L1 node ***/
   
   for(i=0; i<ptr_root->N1; i++)
    {
       for(j=0; j<dim/*DIM*/; j++) current_center[j] = ptr_root->L1CentersRadii[i*(dim+1)+j] ;

       child_start = ptr_root->L1node_StartChild[i] ;
       child_end = child_start + ptr_root->L1node_NumChldrn[i] ;
       for(k=child_start; k<child_end; k++)
         {
           for(j=0; j<dim/*DIM*/; j++) child_center[j]= ptr_root->L2CentersRadii[i][k*(dim+1)+j] ;
           dist_ss = calc_dist_square(dim, current_center, child_center) ;
           dist_ss=sqrt(dist_ss);
           tmp = ptr_root->L2CentersRadii[i][k*(dim+1)+dim] + dist_ss ;
           if((tmp<=R1)  && (ptr_root->L1CentersRadii[i*(dim+1)+dim] < tmp) )
               ptr_root->L1CentersRadii[i*(dim+1)+dim] = tmp ;
           else printf("\n ERROR \n");
         }
    }

   /*** free memory: buf, cluster_assign ***/
   free(datum);
   free(buf) ;
   free(cluster_assign);
   free(cluster_center);
   free(current_center);
   free(child_center) ;

   return 0;
}



/********************************************* End of function Growing tree grow()***************************/




#endif
