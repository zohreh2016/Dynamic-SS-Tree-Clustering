#if !defined BISEC_CLUSTERING_RADIUS_H
#define BISEC_CLUSTERING_RADIUS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include "math.h"
#include "string.h"
//#include "SphereTree4L.h"



#define min(x,y)        ((x)>(y) ? (y) : (x))


/****************************************************************
 o                                       L0
 
 o              o                o                o            L1
 
 o o o          o o o o         o o o o        o o o o o        L2
 
 ooo oo ooo     ooo ooo oo oo   oo oo oo oo    ooo oo oo oo oo    L3=Leaves=data
 
 In this example, NumLevels=4, N1=4, N2=16, N3=37
 L1node_NumChldrn[4]:  3,        4,         4,        5
 L2node_NumChldrn[16]:3,2,3,  3,3,2,2,   2,2,2,2,  3,2,2,2,2
 ****************************************************************/

struct GrowingTree3L {    /* Number of levels is 2 or 3: root is at L0        */
    
    int dim,       /* dimension of each datum                     */
    NumData,   /* total num of data items                     */
    //NumLevels, /* num of levels of the tree                   */
    N1,        /* N1: num of L1 clusters. N0=1 and is omitted */
    N2;       /* N2: total num of L2 clusters.               */
    //leafmax;   /* leafmax: maximum num data points in any cluster. */
    
    //double R1,R2;   /* threshold radii of L1, L2 clusters respectively  */
    
    int *L1node_NumChldrn,  /* L1…[N1]: num of children of each L1 node.
                             Total number of L1 children = N2 =
                             L1NumChldrn[0]+L1NumChldrn[1]+...+L1NumChldrn[N1-1] */
       *L1node_StartChild,  /* L1…[N1]: L2 cluster index of the
                          start child of each L1 node                       */
       *L1node_DataSize,   /* L1…[N1]: num of data of each L1 cluster */
       *L1node_StartDatum, /* L1…[N1]: start datum of each L1 cluster */
    
      **L2node_DataSize,   /* L2…[N2][varies]: num of data of each L2 cluster   */
      **L2node_StartDatum, /* L2…[N2][varies]: start datum of each L2 cluster   */
    
     /*** L2node_DataSize[i][j]: num of data in j-th L2 cluster  ***
     *** which is a child of the i-th L1 cluster                ***/
    
    //////zz
    **L2node_NumChldrn, /* L2…[N2][varies]: num of children of each L1 node.
                            Total number of L1 children = N3 =
                            L1NumChldrn[0]+L1NumChldrn[1]+...+L1NumChldrn[N2-1] */
   
    
    **L2node_StartChild;/* L2…[N2][varies]: L2 cluster index of the
                        start child of each L2 node                         */

    ///////zz
    double *L1CentersRadii, /* L1centerRadii[N1*(DIM+1)]               */
    **L2CentersRadii, /* L2centerRadii[N2*(DIM+1)]               */    /* L1centerRadii[N1][varies*(DIM+1)] */
    *data;
}; /**********************  end of struct LPT4L  **********************/




int GrowingTree_construc (int dim, double R1,double R2, int ndata, double *data,struct GrowingTree3L *ptr_root);

int GrowingTree_grow (int dim,double R1,double R2,int ndata,double *data, struct GrowingTree3L *ptr_root);

int allneighbors_search(int dim, int ndata, double *data, struct GrowingTree3L *ptr_root,
                        double *query, double delta, int *num_outdata, int **outdata);
/* Returns the number of points to which the distances from the query are calculated */


int clustering_innerpoints(double R2,double *StoreBInnerArray,int *AssignInnerBigCluster,int dim,struct GrowingTree3L *ptr_root,int Index_Inner,int *BigClusterOuterAdd,int *StartInsideBigC);
int clustering_outerpoints(int R1,int R2,double *StoreBOuterArray,int *AssignInnerBigCluster ,int dim,struct GrowingTree3L *ptr_root,int ndata,int Index_Outer);
int SplitR_construc(double R1,int i0,int im, int k_split, double *data,int ndata, int dim, int * cluster_assign, struct GrowingTree3L *ptr_root,int flagfirstchecking,double * cluster_mean);



/******************** Functions for clustering *********************/

 double calc_dist_square(int dim, double *datum1, double *datum2);




#endif
