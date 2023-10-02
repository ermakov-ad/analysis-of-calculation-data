#include "math_utils.h"
 
using namespace std;
 
 int solve_system( struct Parameters2d *params2d, double **A, double *b, const int n, const double eps, double *x )
{ 
    int i,j,m=0;
    double norma;
    double w = 1.9;
    double a;
    //double xn[n];
    double *xn; // array with the information about permutations performed during LU-decomposition process
    if (BOUN_TYPE == 0 ){
    get_memory_for_1D_double_array( params2d->cells_number_x * params2d->cells_number_y * params2d->cells_number_z, &xn );
    double S, Sn;
 
    for( i = 0; i < params2d->cells_number_x * params2d->cells_number_y * params2d->cells_number_z; i++){          
            xn[i] = x[i];
            //x[i] = 0.0;
    }
    do{
        m++;
        norma = 0;
        S = 0;
        Sn = 0;
        for(i = 0; i < n; i++){
            x[i] = b[i];
            //for(j = 0; j < n; j++){                  
                //if(i!=j)
                //    x[i] = x[i] - A[i][j] * x[j];
                //printf("%lf\n", x[i]);
            if (i == 0)
                x[i] = x[i]  - A[i][4] * x[i + 1]
                - A[i][5] * x[i + params2d->cells_number_x] - A[i][6] * x[i + params2d->cells_number_x * params2d->cells_number_y];
            else
            if (i >=1 && i < params2d->cells_number_x )
                x[i] = x[i]   - A[i][2] * x[i - 1] - A[i][4] * x[i + 1]
                - A[i][5] * x[i + params2d->cells_number_x] - A[i][6] * x[i + params2d->cells_number_x * params2d->cells_number_y];
            else
            if (i >= params2d->cells_number_x && i < params2d->cells_number_x * params2d->cells_number_y)
                x[i] = x[i] - A[i][1] * x[i - params2d->cells_number_x] - A[i][2] * x[i - 1] - A[i][4] * x[i + 1]
                - A[i][5] * x[i + params2d->cells_number_x] - A[i][6] * x[i + params2d->cells_number_x * params2d->cells_number_y];
            else
            if ( i == n - 1)
                x[i] = x[i] - A[i][0] * x[i - params2d->cells_number_x * params2d->cells_number_y] - A[i][1] * x[i - params2d->cells_number_x] - A[i][2] * x[i - 1];
            else
            if ( i < n - 1 && i > n - params2d->cells_number_x - 1)
                x[i] = x[i] - A[i][0] * x[i - params2d->cells_number_x * params2d->cells_number_y] - A[i][1] * x[i - params2d->cells_number_x] - A[i][2] * x[i - 1] - A[i][4] * x[i + 1];
            else
            if ( i > n - params2d->cells_number_x * params2d->cells_number_y - 1 )
                x[i] = x[i] - A[i][0] * x[i - params2d->cells_number_x * params2d->cells_number_y] - A[i][1] * x[i - params2d->cells_number_x] - A[i][2] * x[i - 1] - A[i][4] * x[i + 1]
                - A[i][5] * x[i + params2d->cells_number_x];
            else      
                x[i] = x[i] - A[i][0] * x[i - params2d->cells_number_x * params2d->cells_number_y] - A[i][1] * x[i - params2d->cells_number_x] - A[i][2] * x[i - 1] - A[i][4] * x[i + 1]
                - A[i][5] * x[i + params2d->cells_number_x] - A[i][6] * x[i + params2d->cells_number_x * params2d->cells_number_y];
                //printf("%lf %d\n", x[i], i);

            //}
            //printf("before%lf %lf\n", x[i], A[i][i] );
            x[i] /= A[i][3];                  
            x[i] = w * x[i] + ( 1 - w ) * xn[i];
            S += fabs(x[i]);
            Sn += fabs(xn[i]);
            //printf("after%lf %lf\n", x[i], xn[i] );               
            //if( fabs( (x[i] - xn[i] + 1.e-10) / (x[i] ) ) > norma )
            //    norma = fabs( (x[i] - xn[i] + 1.e-10) / ( x[i] ) );
            xn[i] = x[i];
            
        }
        norma = fabs(S-Sn)/(S+1.e-10);
        //printf("norm = %lf \n", norma);
        
    }
    while(norma > eps||m<5);
    }
    else
    if (BOUN_TYPE == 1){
    get_memory_for_1D_double_array( params2d->cells_number_x * params2d->cells_number_y * params2d->cells_number_z, &xn );
    double S, Sn;
 
    for( i = 0; i < params2d->cells_number_x * params2d->cells_number_y * params2d->cells_number_z; i++){          
            //xn[i] = 0.0;//x[i];
        xn[i] = x[i];
            //x[i] = 0.0;
    }
    do{
        m++;
        norma = 0;
        S = 0;
        Sn = 0;
                for ( int i = 1; i < params2d->cells_number_x - 1; i++ )
		    for ( int j = 1; j <  params2d->cells_number_y - 1; j++ )
		        for ( int k = 1; k < params2d->cells_number_z - 1; k++ ){
                            int cntr = k * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i;
            x[cntr] = b[cntr];
            //for(j = 0; j < n; j++){                  
                //if(i!=j)
                //    x[i] = x[i] - A[i][j] * x[j];
                //printf("%lf\n", x[i]);
            //zyxcxyz//
            if ( i == 1 )
                x[cntr] -= ( A[cntr][2]
            * x[k * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + params2d->cells_number_x - 2] + A[cntr][4] * 
                x[k * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i + 1] );
            else
                if ( i == params2d->cells_number_x - 2 )
                x[cntr] -= ( A[cntr][2]
            * x[k * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i - 1] + A[cntr][4] * 
                x[k * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + 1] );
                else
                x[cntr] -= ( A[cntr][2]
            * x[k * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i - 1] + A[cntr][4] * 
                x[k * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i + 1] );

            if ( j == 1 )
                x[cntr] -= ( A[cntr][1]
            * x[k * params2d->cells_number_x * params2d->cells_number_y + ( params2d->cells_number_y - 2 ) * params2d->cells_number_x + i] + A[cntr][5] * 
                x[k * params2d->cells_number_x * params2d->cells_number_y + ( j + 1 ) * params2d->cells_number_x + i] );
            else
                if ( j == params2d->cells_number_y - 2 )
                x[cntr] -= ( A[cntr][1]
            * x[k * params2d->cells_number_x * params2d->cells_number_y + ( j - 1 ) * params2d->cells_number_x + i] + A[cntr][5] * 
                x[k * params2d->cells_number_x * params2d->cells_number_y + 1 * params2d->cells_number_x + i] );
                else
                x[cntr] -= ( A[cntr][1]
            * x[k * params2d->cells_number_x * params2d->cells_number_y + ( j - 1) * params2d->cells_number_x + i ] + A[cntr][5] * 
                x[k * params2d->cells_number_x * params2d->cells_number_y + ( j + 1 ) * params2d->cells_number_x + i ] );

            if ( k == 1 )
                x[cntr] -= ( A[cntr][0]
            * x[( params2d->cells_number_z - 2 ) * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i] + A[cntr][6] * 
                x[( k + 1 ) * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i] );
            else
                if ( k == params2d->cells_number_z - 2 )
                x[cntr] -= ( A[cntr][0]
            * x[( k - 1 ) * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i] + A[cntr][6] * 
                x[1 * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i] );
                else
                x[cntr] -= ( A[cntr][0]
            * x[( k - 1 ) * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i] + A[cntr][6] * 
                x[( k + 1 ) * params2d->cells_number_x * params2d->cells_number_y + j * params2d->cells_number_x + i] );


            //}
            //printf("before%lf %lf\n", x[cntr], A[cntr][3] );
            x[cntr] /= A[cntr][3];                  
            x[cntr] = w * x[cntr] + ( 1 - w ) * xn[cntr];
            S += fabs(x[cntr]);
            Sn += fabs(xn[cntr]);
            //printf("after%lf %lf %d %d %d %d\n", x[cntr], xn[cntr], m, i , j , k );               
            //if( fabs( (x[i] - xn[i] + 1.e-10) / (x[i] ) ) > norma )
            //    norma = fabs( (x[i] - xn[i] + 1.e-10) / ( x[i] ) );
            xn[cntr] = x[cntr];
            
        }
        norma = fabs(S-Sn)/(S+1.e-10);
        //printf("norm = %lf \n", norma);
        
    }
    while(norma > eps||m<5);
    }
    //double max = -1.e20;
    //double min = 1.e20;
    //for( i = 0; i < params2d->cells_number_x * params2d->cells_number_y * params2d->cells_number_z; i++){          
    //        if (x[i] < min && fabs(x[i]) > 0.001)
    //            min = x[i];
    //        if (x[i] > max)
    //            max = x[i];
    //}
    //printf("%lf %lf\n", min, max);
    //for( i = 0; i < params2d->cells_number_x * params2d->cells_number_y * params2d->cells_number_z; i++){  
    //    if (max <= 0)
    //        x[i] = x[i] + ( max - min ) / 2;
    //    if (min > 0.01)
    //        x[i] = x[i] - ( max - min ) / 2;
    //}
    printf("norm = %lf\n", norma);
    free( xn );
    return 0;
  }