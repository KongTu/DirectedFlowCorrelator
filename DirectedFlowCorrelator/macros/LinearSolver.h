#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"

const int MAX = 4;
double determ(double a[MAX][MAX],int n) {
  double det=0;
  int p, h, k, i, j;
  double temp[MAX][MAX];
  if(n==1) {
    return a[0][0];
  } else if(n==2) {
    det=(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
    return det;
  } else {
    for(p=0;p<n;p++) {
      h = 0;
      k = 0;
      for(i=1;i<n;i++) {
        for( j=0;j<n;j++) {
          if(j==p) {
            continue;
          }
          temp[h][k] = a[i][j];
          k++;
          if(k==n-1) {
            h++;
            k = 0;
          }
        }
      }
      det=det+a[0][p]*pow(-1,p)*determ(temp,n-1);
    }
    return det;
  }
}
const int MAX1=3;
double determ3(double a[MAX1][MAX1],int n) {
  double det=0;
  int p, h, k, i, j;
  double temp[MAX1][MAX1];
  if(n==1) {
    return a[0][0];
  } else if(n==2) {
    det=(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
    return det;
  } else {
    for(p=0;p<n;p++) {
      h = 0;
      k = 0;
      for(i=1;i<n;i++) {
        for( j=0;j<n;j++) {
          if(j==p) {
            continue;
          }
          temp[h][k] = a[i][j];
          k++;
          if(k==n-1) {
            h++;
            k = 0;
          }
        }
      }
      det=det+a[0][p]*pow(-1,p)*determ3(temp,n-1);
    }
    return det;
  }
}

// double determ3(double a[4][3][3],int n,int index) {
//   double det=0;
//   int p, h, k, i, j;
//   double temp[4][3][3];
//   if(n==1) {
//     return a[index][0][0];
//   } else if(n==2) {
//     det=(a[index][0][0]*a[index][1][1]-a[index][0][1]*a[index][1][0]);
//     return det;
//   } else {
//     for(p=0;p<n;p++) {
//       h = 0;
//       k = 0;
//       for(i=1;i<n;i++) {
//         for( j=0;j<n;j++) {
//           if(j==p) {
//             continue;
//           }
//           temp[index][h][k] = a[index][i][j];
//           k++;
//           if(k==n-1) {
//             h++;
//             k = 0;
//           }
//         }
//       }
//       det=det+a[index][0][p]*pow(-1,p)*determ3(temp,n-1,index);
//     }
//     return det;
//   }
// }

vector<double> GetDmesonSigBkg( vector<double> frac_mc, vector<double> v1_data, vector<double> v1_data_error ){

	double A[4];
	double B[4];
	double C[4];
	double D[4];

	A[0] = frac_mc[0];//f_obs_sig_d0
	A[1] = frac_mc[3];//f_obs_swap_d0bar
	A[2] = frac_mc[4];//f_sb_sig_d0
	A[3] = frac_mc[7];//f_sb_swap_d0bar

	B[0] = frac_mc[1];//f_obs_swap_d0
	B[1] = frac_mc[2];//f_obs_sig_d0bar
	B[2] = frac_mc[5];//f_sb_swap_d0
	B[3] = frac_mc[6];//f_sb_sig_d0bar

	C[0] = 1-A[0]-B[0];
	C[1] = 0.0;
	C[2] = 1-A[2]-B[2];
	C[3] = 0.0;

	D[0] = 0.0;
	D[1] = 1-A[1]-B[1];
	D[2] = 0.0;
	D[3] = 1-A[3]-B[3];

	double Y[4];

	Y[0] = v1_data[0];//v1_obs_d0
	Y[1] = v1_data[1];//v1_obs_d0bar
	Y[2] = v1_data[2];//v1_side_d0
	Y[3] = v1_data[3];//v1_side_d0bar

	double Y_e[4];
	Y_e[0] = v1_data_error[0];//error_v1_obs_d0
	Y_e[1] = v1_data_error[1];//error_v1_obs_d0bar
	Y_e[2] = v1_data_error[2];//error_v1_side_d0
	Y_e[3] = v1_data_error[3];//error_v1_side_d0bar

	double DD[4][4];
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){

			if(i == 0 ) DD[j][i] = A[j];
			if(i == 1 ) DD[j][i] = B[j];
			if(i == 2 ) DD[j][i] = C[j];
			if(i == 3 ) DD[j][i] = D[j];
		}
	}

	double D0[4][4];
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){

			if(i == 0 ) D0[j][i] = Y[j];
			if(i == 1 ) D0[j][i] = B[j];
			if(i == 2 ) D0[j][i] = C[j];
			if(i == 3 ) D0[j][i] = D[j];
		}
	}
	double D1[4][4];
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){

			if(i == 0 ) D1[j][i] = A[j];
			if(i == 1 ) D1[j][i] = Y[j];
			if(i == 2 ) D1[j][i] = C[j];
			if(i == 3 ) D1[j][i] = D[j];
		}
	}
	double D2[4][4];
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){

			if(i == 0 ) D2[j][i] = A[j];
			if(i == 1 ) D2[j][i] = B[j];
			if(i == 2 ) D2[j][i] = Y[j];
			if(i == 3 ) D2[j][i] = D[j];
		}
	}

	double D3[4][4];
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){

			if(i == 0 ) D3[j][i] = A[j];
			if(i == 1 ) D3[j][i] = B[j];
			if(i == 2 ) D3[j][i] = C[j];
			if(i == 3 ) D3[j][i] = Y[j];
		}
	}

	double DDD =  determ(DD, 4);
	double D00 =  determ(D0, 4);
	double D11 =  determ(D1, 4);
	double D22 =  determ(D2, 4);
	double D33 =  determ(D3, 4);

	vector<double> tmp;
	tmp.push_back( D00/DDD );//sig_d0
	tmp.push_back( D11/DDD );//sig_d0bar
	tmp.push_back( D22/DDD );//bkg_d0
	tmp.push_back( D33/DDD );//bkg_d0bar

// Error for all:

	double Y1_assoc1[3][3];
    const int index1[3] = {1,2,3};
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
	
			if(i == 0) Y1_assoc1[j][i] = B[ index1[j] ];
			if(i == 1) Y1_assoc1[j][i] = C[ index1[j] ];
			if(i == 2) Y1_assoc1[j][i] = D[ index1[j] ];
		}
	}

	double Y1_assoc2[3][3];
	const int index2[3] = {0,2,3};
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){

			if(i == 0) Y1_assoc2[j][i] = B[ index2[j] ];
			if(i == 1) Y1_assoc2[j][i] = C[ index2[j] ];
			if(i == 2) Y1_assoc2[j][i] = D[ index2[j] ];
		}
	}
	double Y1_assoc3[3][3];
	const int index3[3] = {0,1,3};
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){

			if(i == 0) Y1_assoc3[j][i] = B[ index3[j] ];
			if(i == 1) Y1_assoc3[j][i] = C[ index3[j] ];
			if(i == 2) Y1_assoc3[j][i] = D[ index3[j] ];
		}
	}
	double Y1_assoc4[3][3];
	const int index4[3] = {0,1,2};
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){

			if(i == 0) Y1_assoc4[j][i] = B[ index4[j] ];
			if(i == 1) Y1_assoc4[j][i] = C[ index4[j] ];
			if(i == 2) Y1_assoc4[j][i] = D[ index4[j] ];
		}
	}
	
	double Y2_assoc1[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
		
			if(i == 0) Y2_assoc1[j][i] = A[ index1[j] ];
			if(i == 1) Y2_assoc1[j][i] = C[ index1[j] ];
			if(i == 2) Y2_assoc1[j][i] = D[ index1[j] ];
		}
	}
	double Y2_assoc2[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
		
			if(i == 0) Y2_assoc2[j][i] = A[ index2[j] ];
			if(i == 1) Y2_assoc2[j][i] = C[ index2[j] ];
			if(i == 2) Y2_assoc2[j][i] = D[ index2[j] ];
		}
	}
	double Y2_assoc3[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
		
			if(i == 0) Y2_assoc3[j][i] = A[ index3[j] ];
			if(i == 1) Y2_assoc3[j][i] = C[ index3[j] ];
			if(i == 2) Y2_assoc3[j][i] = D[ index3[j] ];
		}
	}
	double Y2_assoc4[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
					
			if(i == 0) Y2_assoc4[j][i] = A[ index4[j] ];
			if(i == 1) Y2_assoc4[j][i] = C[ index4[j] ];
			if(i == 2) Y2_assoc4[j][i] = D[ index4[j] ];
		}
	}
	
	double Y3_assoc1[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
					
			if(i == 0) Y3_assoc1[j][i] = A[ index1[j] ];
			if(i == 1) Y3_assoc1[j][i] = B[ index1[j] ];
			if(i == 2) Y3_assoc1[j][i] = D[ index1[j] ];
		}
	}
	double Y3_assoc2[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
		
			if(i == 0) Y3_assoc2[j][i] = A[ index2[j] ];
			if(i == 1) Y3_assoc2[j][i] = B[ index2[j] ];
			if(i == 2) Y3_assoc2[j][i] = D[ index2[j] ];
		}
	}
	double Y3_assoc3[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
					
			if(i == 0) Y3_assoc3[j][i] = A[ index3[j] ];
			if(i == 1) Y3_assoc3[j][i] = B[ index3[j] ];
			if(i == 2) Y3_assoc3[j][i] = D[ index3[j] ];
		}
	}
	double Y3_assoc4[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){

			if(i == 0) Y3_assoc4[j][i] = A[ index4[j] ];
			if(i == 1) Y3_assoc4[j][i] = B[ index4[j] ];
			if(i == 2) Y3_assoc4[j][i] = D[ index4[j] ];
		}
	}

	double Y4_assoc1[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){

			if(i == 0) Y4_assoc1[j][i] = A[ index1[j] ];
			if(i == 1) Y4_assoc1[j][i] = B[ index1[j] ];
			if(i == 2) Y4_assoc1[j][i] = C[ index1[j] ];
		}
	}
	double Y4_assoc2[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){

			if(i == 0) Y4_assoc2[j][i] = A[ index2[j] ];
			if(i == 1) Y4_assoc2[j][i] = B[ index2[j] ];
			if(i == 2) Y4_assoc2[j][i] = C[ index2[j] ];
		}
	}
	double Y4_assoc3[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
					
			if(i == 0) Y4_assoc3[j][i] = A[ index3[j] ];
			if(i == 1) Y4_assoc3[j][i] = B[ index3[j] ];
			if(i == 2) Y4_assoc3[j][i] = C[ index3[j] ];
		}
	}
	double Y4_assoc4[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){

			if(i == 0) Y4_assoc4[j][i] = A[ index4[j] ];
			if(i == 1) Y4_assoc4[j][i] = B[ index4[j] ];
			if(i == 2) Y4_assoc4[j][i] = C[ index4[j] ];
		}
	}


//X1 error

	double D00_e = determ3(Y1_assoc1,3);
	double D11_e = determ3(Y1_assoc2,3);
	double D22_e = determ3(Y1_assoc3,3);
	double D33_e = determ3(Y1_assoc4,3);

	D00_e = D00_e/DDD;
	D11_e = D11_e/DDD;
	D22_e = D22_e/DDD;
	D33_e = D33_e/DDD;

	double final_error_1 = sqrt(D00_e*D00_e*Y_e[0]*Y_e[0] + D11_e*D11_e*Y_e[1]*Y_e[1] + D22_e*D22_e*Y_e[2]*Y_e[2] +D33_e*D33_e*Y_e[3]*Y_e[3]);

	tmp.push_back(final_error_1);

//X2 error

	D00_e = determ3(Y2_assoc1,3);
	D11_e = determ3(Y2_assoc2,3);
	D22_e = determ3(Y2_assoc3,3);
	D33_e = determ3(Y2_assoc4,3);

	D00_e = D00_e/DDD;
	D11_e = D11_e/DDD;
	D22_e = D22_e/DDD;
	D33_e = D33_e/DDD;

	double final_error_2 = sqrt(D00_e*D00_e*Y_e[0]*Y_e[0] + D11_e*D11_e*Y_e[1]*Y_e[1] + D22_e*D22_e*Y_e[2]*Y_e[2] +D33_e*D33_e*Y_e[3]*Y_e[3]);

	tmp.push_back(final_error_2);

//X3 error

	D00_e = determ3(Y3_assoc1,3);
	D11_e = determ3(Y3_assoc2,3);
	D22_e = determ3(Y3_assoc3,3);
	D33_e = determ3(Y3_assoc4,3);

	D00_e = D00_e/DDD;
	D11_e = D11_e/DDD;
	D22_e = D22_e/DDD;
	D33_e = D33_e/DDD;

	double final_error_3 = sqrt(D00_e*D00_e*Y_e[0]*Y_e[0] + D11_e*D11_e*Y_e[1]*Y_e[1] + D22_e*D22_e*Y_e[2]*Y_e[2] +D33_e*D33_e*Y_e[3]*Y_e[3]);

	tmp.push_back(final_error_3);

//X4 error

	D00_e = determ3(Y4_assoc1,3);
	D11_e = determ3(Y4_assoc2,3);
	D22_e = determ3(Y4_assoc3,3);
	D33_e = determ3(Y4_assoc4,3);

	D00_e = D00_e/DDD;
	D11_e = D11_e/DDD;
	D22_e = D22_e/DDD;
	D33_e = D33_e/DDD;

	double final_error_4 = sqrt(D00_e*D00_e*Y_e[0]*Y_e[0] + D11_e*D11_e*Y_e[1]*Y_e[1] + D22_e*D22_e*Y_e[2]*Y_e[2] +D33_e*D33_e*Y_e[3]*Y_e[3]);

	tmp.push_back(final_error_4);

	return tmp;

}









