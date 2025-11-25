#include "mex.h"
#include "matrix.h"
#include <algorithm>
#include <cmath>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check inputs
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("mpmsc:nrhs", "Three inputs required: X, r, w");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("mpmsc:nlhs", "One output required");
    }
    
    // Get input matrix X (m x n) - rows are spectra
    const double* X = mxGetPr(prhs[0]);
    int m = mxGetM(prhs[0]);  // number of spectra
    int n = mxGetN(prhs[0]);  // number of channels
    
    // Get reference vector r (length n) - must be row vector
    const double* r = mxGetPr(prhs[1]);
    int r_len = std::max(mxGetM(prhs[1]), mxGetN(prhs[1]));
    
    if (r_len != n) {
        mexErrMsgIdAndTxt("mpmsc:dimMismatch", 
                         "Reference vector length must match number of columns in X");
    }
    
    // Get window size w
    int w = (int)mxGetScalar(prhs[2]);
    
    // Create output matrix
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    double* Z = mxGetPr(plhs[0]);
    
    // Allocate temporary arrays for moving statistics
    double* Xmean = new double[m * n];
    double* rmean = new double[n];
    double* rsum = new double[n];
    double* esum = new double[n];
    double* Denom = new double[n];
    double* Numer = new double[m * n];
    double* Slope = new double[m * n];
    
    // Temporary window buffers
    double* window_buf = new double[w];
    
    // ------------------------------------------------
    // PRE-COMPUTE MOVING SUMS AND MEANS
    // ------------------------------------------------
    
    // Compute Xmean = movmean(X, w, 2) - moving mean along columns (dim 2)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int start = std::max(0, j - (w-1)/2);
            int end = std::min(n - 1, j + (w-1)/2);
            int win_size = end - start + 1;
            
            double sum = 0.0;
            for (int k = start; k <= end; k++) {
                sum += X[i + k * m];
            }
            Xmean[i + j * m] = sum / win_size;
        }
    }
    
    // Compute rmean = movmean(r, w)
    for (int j = 0; j < n; j++) {
        int start = std::max(0, j - (w-1)/2);
        int end = std::min(n - 1, j + (w-1)/2);
        int win_size = end - start + 1;
        
        double sum = 0.0;
        for (int k = start; k <= end; k++) {
            sum += r[k];
        }
        rmean[j] = sum / win_size;
    }
    
    // Compute rsum = movsum(r, w)
    for (int j = 0; j < n; j++) {
        int start = std::max(0, j - (w-1)/2);
        int end = std::min(n - 1, j + (w-1)/2);
        
        double sum = 0.0;
        for (int k = start; k <= end; k++) {
            sum += r[k];
        }
        rsum[j] = sum;
    }
    
    // Compute esum = movsum(ones(1,n), w)
    for (int j = 0; j < n; j++) {
        int start = std::max(0, j - (w-1)/2);
        int end = std::min(n - 1, j + (w-1)/2);
        esum[j] = (double)(end - start + 1);
    }
    
    // ------------------------------------------------
    // COMPUTE THE MOVING SLOPE
    // ------------------------------------------------
    
    // Compute Denom = movsum(r.*r, w) - 2*rmean.*rsum + rmean.*(rmean.*esum)
    for (int j = 0; j < n; j++) {
        int start = std::max(0, j - (w-1)/2);
        int end = std::min(n - 1, j + (w-1)/2);
        
        double sum_r_sq = 0.0;
        for (int k = start; k <= end; k++) {
            sum_r_sq += r[k] * r[k];
        }
        Denom[j] = sum_r_sq - 2.0 * rmean[j] * rsum[j] + rmean[j] * rmean[j] * esum[j];
    }
    
    // Compute Numer = movsum(X.*r, w, 2) - movsum(X,w,2).*rmean - Xmean.*rsum + Xmean.*(rmean.*esum)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int start = std::max(0, j - (w-1)/2);
            int end = std::min(n - 1, j + (w-1)/2);
            
            // movsum(X.*r, w, 2)
            double sum_xr = 0.0;
            for (int k = start; k <= end; k++) {
                sum_xr += X[i + k * m] * r[k];
            }
            
            // movsum(X, w, 2)
            double sum_x = 0.0;
            for (int k = start; k <= end; k++) {
                sum_x += X[i + k * m];
            }
            
            Numer[i + j * m] = sum_xr - sum_x * rmean[j] - Xmean[i + j * m] * rsum[j] + 
                               Xmean[i + j * m] * rmean[j] * esum[j];
        }
    }
    
    // Compute Slope = Denom ./ Numer (element-wise division)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (std::abs(Numer[i + j * m]) < 1e-10) {
                Slope[i + j * m] = 0.0;  // Avoid division by zero
            } else {
                Slope[i + j * m] = Denom[j] / Numer[i + j * m];
            }
        }
    }
    
    // ------------------------------------------------
    // CORRECTED SPECTRA
    // ------------------------------------------------
    // Z = (X - Xmean) .* Slope + rmean
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            Z[i + j * m] = (X[i + j * m] - Xmean[i + j * m]) * Slope[i + j * m] + rmean[j];
        }
    }
    
    // Clean up
    delete[] Xmean;
    delete[] rmean;
    delete[] rsum;
    delete[] esum;
    delete[] Denom;
    delete[] Numer;
    delete[] Slope;
    delete[] window_buf;
}
