#include "mex.h"
#include "matrix.h"
#include <algorithm>
#include <cmath>

// Helper function: Compute mean of a vector
double computeMean(const double* vec, int len) {
    double sum = 0.0;
    for (int i = 0; i < len; i++) {
        sum += vec[i];
    }
    return sum / len;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check number of inputs and outputs
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("pmsc2_helland:nrhs", "Three inputs required: X, r, v");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("pmsc2_helland:nlhs", "One output required");
    }
    
    // Input matrix X (m x n) 
    const double* X = mxGetPr(prhs[0]);
    int m = mxGetM(prhs[0]);  // number of spectra
    int n = mxGetN(prhs[0]);  // number of channels (e.g., wavelengths)
    
    // Reference vector r 
    const double* r = mxGetPr(prhs[1]);
    int r_len = std::max(mxGetM(prhs[1]), mxGetN(prhs[1]));
    // Check length of reference mean
    if (r_len != n) {
        mexErrMsgIdAndTxt("pmsc2_helland:dimMismatch", 
                         "Reference vector length must match number of columns in X");
    }
    
    // Window half-width v
    int v = (int)mxGetScalar(prhs[2]);
    
    // Output matrix
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    double* Z = mxGetPr(plhs[0]);
    
    // Initialize Z: Copy X to Z 
    for (int i = 0; i < m * n; i++) {
        Z[i] = X[i];
    }
    
    // Temporary storage for intermediate calculations
    double* Xband = new double[m * (2*v + 1)];  
    double* rband = new double[2*v + 1];
    double* rc    = new double[2*v + 1];
    double* Xmean = new double[m];
    double* Xc = new double[m * (2*v + 1)];
    
    // Cycle though each channel j
    for (int j = 0; j < n; j++) {
        int start = std::max(0, j - v);
        int end = std::min(n - 1, j + v);
        int band_size = end - start + 1;
        
        // Extract data window or band: Xband and rband
        for (int col = 0; col < band_size; col++) {
            rband[col] = r[start + col];
            for (int row = 0; row < m; row++) {
                Xband[row + col * m] = X[row + (start + col) * m];
            }
        }
        
        // rmean: mean of reference spectrum on band
        double rmean = computeMean(rband, band_size);
        
        // Xmean: mean of each row of X on band
        for (int i = 0; i < m; i++) {
            double sum = 0.0;
            for (int col = 0; col < band_size; col++) {
                sum += Xband[i + col * m];
            }
            Xmean[i] = sum / band_size;
        }
        
        // Xc: Mean-center each row of X on band
        for (int col = 0; col < band_size; col++) {
            for (int i = 0; i < m; i++) {
                Xc[i + col * m] = Xband[i + col * m] - Xmean[i];
            }
        }
        
        // rc: Mean-center reference spectrum on band
        for (int k = 0; k < band_size; k++) {
            rc[k] = rband[k] - rmean;
        }
        
        // Inner product: rc' * rc 
        double rc_dot_rc = 0.0;
        for (int k = 0; k < band_size; k++) {
            rc_dot_rc += rc[k] * rc[k];
        }
        
        // If reference has no variance, set to rmean
        if (std::abs(rc_dot_rc) < 1e-10) {
            for (int i = 0; i < m; i++) {
                Z[i + j * m] = rmean;
            }
            continue;
        }
        
        // b = (Xc * rc) / (rc' * rc) for each row
        double* b = new double[m];
        for (int i = 0; i < m; i++) {
            double dot = 0.0;
            for (int k = 0; k < band_size; k++) {
                dot += Xc[i + k * m] * rc[k];
            }
            b[i] = dot / rc_dot_rc;
        }
        
        // Find which column in band corresponds to j
        int j_in_band = j - start;
        
        // Apply correction: Z[:,j] = rmean + Xc[:,j_in_band] ./ b
        for (int i = 0; i < m; i++) {
            if (std::abs(b[i]) < 1e-10) {
                Z[i + j * m] = rmean;
            } else {
                Z[i + j * m] = rmean + Xc[i + j_in_band * m] / b[i];
            }
        }
        
        delete[] b;
    }
    
    // Clean up
    delete[] Xband;
    delete[] rband;
    delete[] rc;
    delete[] Xmean;
    delete[] Xc;
}
