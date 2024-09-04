
#include <fstream>
#include <cmath>
using namespace std;


#include "Simplex.h"
#include "Base/Random.h"
#include "Base/Estring.h"
#include "Base/Matrix.h"
// Notice, val is not passed by reference
double dabs(double val)
{
    val = pow(val * val, 0.5);
    return val;
}

// void Amoeba(double *Values[], int ndim, CharString & OutputFileName,
//     CharString & InputFileName, double (*funk) (CharString & InFileName,
//         double val), double Limits[], double Noise, double ftol)
// {

//     /**************\
//     * SETUP MATRIX *
//     \**************/

//     double *y = new double[ndim + 1];

//     Matrix < double >p(ndim + 1, ndim + 1);

//     // Spreading Factor of 5%
//     double *Spread = new double[ndim];
//     double minspread = 1000.;

//     for(int i = 0; i < ndim; ++i)
//     {
//         p[0][i] = *Values[i];
//         Spread[i] = dabs(0.05 * p[0][i]);
//         if(Spread[i] < minspread)
//             minspread = Spread[i];
//     }
//     y[0] = funk(InputFileName, Noise);
//     Random inNoise(-minspread, minspread);

//     for(int j = 1; j < ndim + 1; ++j)
//     {
//         for(int i = 0; i < ndim; ++i)
//         {
//             p[j][i] =
//                 p[0][i] + Spread[i] * (double)(ndim / 2 -
//                     j) / (double)ndim + inNoise.Draw();
//             *Values[i] = p[j][i];
//             if(p[j][i] < Limits[i])
//                 p[j][i] = Limits[i];
//             if(p[j][i] > Limits[i + ndim])
//                 p[j][i] = Limits[i + ndim];
//         }
//         y[j] = funk(InputFileName, Noise);
//     }

//     /*  ***************
//     * SIMPLEX PROGRAM
//     * *************** */

//     int mpts = ndim + 1;
//     int iter = 0;
//     int itmax = 10000;

//     double *pbar = new double[ndim];
//     double *prr = new double[ndim];
//     double *pr = new double[ndim];
//     double alpha = 1.0;
//     double beta = 0.5;
//     double gamma = 2.0;
//     double ypr, yprr;
//     streampos pos;
//     fstream FileName(OutputFileName.Ec_str(), ios::out);

//     pos = FileName.tellp();
//     double rtol = 1.0e10;

//     while(iter < 10 || (rtol > ftol && iter < itmax))
//     {
//         ++iter;
//         FileName.seekg(pos);
//         int ilo = 1;
//         int ihi, inhi;

//         if(y[0] > y[1])
//         {
//             ihi = 0;
//             inhi = 1;
//         }
//         else
//         {
//             ihi = 1;
//             inhi = 0;
//         }
//         for(int i = 0; i < mpts; ++i)
//         {
//             if(y[i] < y[ilo])
//                 ilo = i;
//             if(y[i] > y[ihi])
//             {
//                 inhi = ihi;
//                 ihi = i;
//             }
//             else if(y[i] > y[inhi])
//             {
//                 if(i != ihi)
//                     inhi = i;
//             }
//         }
//         for(int jj = 0; jj < ndim; ++jj)
//         {
//             pbar[jj] = 0.;
//         }
//         for(int ii = 0; ii < mpts; ++ii)
//         {
//             if(ii != ihi)
//             {
//                 for(int jj = 0; jj < ndim; ++jj)
//                 {
//                     pbar[jj] = pbar[jj] + p[ii][jj];
//                 }
//             }
//         }
//         for(int j2 = 0; j2 < ndim; ++j2)
//         {
//             pbar[j2] = pbar[j2] / (double)ndim;
//             pr[j2] = (1.0 + alpha) * pbar[j2] - alpha * p[ihi][j2];
//             if(pr[j2] < Limits[j2])
//             {
//                 pr[j2] = Limits[j2];
//             }
//             if(pr[j2] > Limits[j2 + ndim])
//             {
//                 pr[j2] = Limits[j2 + ndim];
//             }
//             *Values[j2] = pr[j2];
//         }
//         ypr = funk(InputFileName, Noise);
//         if(ypr < y[ilo])
//         {
//             for(int j3 = 0; j3 < ndim; ++j3)
//             {
//                 prr[j3] = gamma * pr[j3] + (1.0 - gamma) * pbar[j3];
//                 if(prr[j3] < Limits[j3])
//                     prr[j3] = Limits[j3];
//                 if(prr[j3] > Limits[j3 + ndim])
//                     prr[j3] = Limits[j3 + ndim];
//                 *Values[j3] = prr[j3];
//             }
//             double yprr = funk(InputFileName, Noise);

//             if(yprr < y[ilo])
//             {
//                 for(int j4 = 0; j4 < ndim; ++j4)
//                 {
//                     p[ihi][j4] = prr[j4];
//                 }
//                 y[ihi] = yprr;
//             }
//             else
//             {
//                 for(int j5 = 0; j5 < ndim; ++j5)
//                 {
//                     p[ihi][j5] = pr[j5];
//                 }
//                 y[ihi] = ypr;
//             }
//         }
//         else if(ypr > y[inhi])
//         {
//             if(ypr < y[ihi])
//             {
//                 for(int j6 = 0; j6 < ndim; ++j6)
//                 {
//                     p[ihi][j6] = pr[j6];
//                 }
//                 y[ihi] = ypr;
//             }
//             for(int j7 = 0; j7 < ndim; ++j7)
//             {
//                 prr[j7] = beta * p[ihi][j7] + (1.0 - beta) * pbar[j7];
//                 if(prr[j7] < Limits[j7])
//                 {
//                     prr[j7] = Limits[j7];
//                 }
//                 if(prr[j7] > Limits[j7 + ndim])
//                 {
//                     prr[j7] = Limits[j7 + ndim];
//                 }
//                 *Values[j7] = prr[j7];
//             }
//             yprr = funk(InputFileName, Noise);
//             if(yprr < y[ihi])
//             {
//                 for(int j = 0; j < ndim; ++j)
//                 {
//                     p[ihi][j] = prr[j];
//                 }
//                 y[ihi] = yprr;
//             }
//             else
//             {
//                 for(int i = 0; i < mpts; ++i)
//                 {
//                     if(i != ilo)
//                     {
//                         for(int j = 0; j < ndim; ++j)
//                         {
//                             pr[j] = 0.5 * (p[i][j] + p[ilo][j]);
//                             if(pr[j] < Limits[j])
//                                 pr[j] = Limits[j];
//                             if(pr[j] > Limits[j + ndim])
//                                 pr[j] = Limits[j + ndim];
//                             *Values[j] = pr[j];
//                             p[i][j] = pr[j];
//                         }
//                         y[i] = funk(InputFileName, Noise);
//                     }
//                 }
//             }
//         }
//         else
//         {
//             if(true)
//             {
//                 for(int j = 0; j < ndim; ++j)
//                 {
//                     p[ihi][j] = pr[j];
//                 }
//                 y[ihi] = ypr;
//             }
//         }
//         for(int j = 0; j < ndim; ++j)
//         {
//             *Values[j] = p[ilo][j];     // Keep values to low ones
//         }

//         for(int i2 = 0; i2 < ndim; i2 += 5)
//         {
//             for(int j = 0; j < 5 && i2 + j < ndim; ++j)
//             {
//                 cout << i2 + j << " : " << *Values[i2 + j] << " \t";
//                 FileName << *Values[i2 + j] << endl;
//             }
//             cout << endl;
//         }
//         rtol = 2. * dabs(y[ihi] - y[ilo]) / (dabs(y[ihi]) + dabs(y[ilo]));
//         cout << "\tSum of Squares : " << y[ilo] << " \t Tolerance : "
//             << rtol << " \tIteration : " << iter << endl;
        
//         if(iter == itmax)
//         {
//             cout << "Maximum no. of iterations exceeded" << endl;
//             return;
//         }
//     }
    
//     cout << "Done." << endl;
//     delete[]y;
//     delete[]pbar;
//     delete[]prr;
//     delete[]pr;
//     delete[]Spread;
//     return;
// }
