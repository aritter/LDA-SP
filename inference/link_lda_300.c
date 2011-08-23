/*
This program was automatically generated using:
     __  ____   ____
    / / / /  \ / __/  HBC: The Hierarchical Bayes Compiler
   / /_/ / / // /     http://hal3.name/HBC/
  / __  / --</ /      
 / / / / /  / /___    Version 0.7 beta
 \/ /_/____/\____/    

HBC is a freely available compiler for statistical models.  This generated
code can be built using the following command:

  gcc -O3 -lm stats.c samplib.c link_lda_300.c -o link_lda_300.out

The hierarchical model that this code reflects is:

alpha		~ Gam(0.1,1)
eta		~ Gam(0.1,1)
beta1_{k}	~ DirSym(eta, V1)		, k \in [1,K]
beta2_{k}	~ DirSym(eta, V2)		, k \in [1,K]
theta_{d}	~ DirSym(alpha, K)		, d \in [1,D]
z1_{d,n}	~ Mult(theta_{d})		, d \in [1,D] , n \in [1,N_{d}]	
z2_{d,n}	~ Mult(theta_{d})		, d \in [1,D] , n \in [1,N_{d}]	
arg1_{d,n}	~ Mult(beta1_{z1_{d,n}})	, d \in [1,D] , n \in [1,N_{d}]
arg2_{d,n}	~ Mult(beta2_{z2_{d,n}})	, d \in [1,D] , n \in [1,N_{d}]

Generated using the command:

  hbc compile --define alpha 0.1 --define eta 0.1 --dump 50 z1 z2 ; --define K 300 --loadD /homes/gws/aritter/predType/hbc/data/arg1 arg1 V1 D N ; --loadD /homes/gws/aritter/predType/hbc/data/arg2 arg2 V2 D N ; --collapse theta --collapse beta1 --collapse beta2 --iter 1000 /homes/gws/aritter/predType/hbc/models/link-sp-lda.hier link_lda_300.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stats.h"


/**************************** SAMPLING ****************************/

void resample_post_theta(int D, int K, int* N, double** post_theta, int** z1, int** z2) {
  int d_13;
  double* tmpSP10;
  int n_4;
  double* tmpSP11;
  int n_22;
  double* vec_var_0;
  int dvv_loop_var_1;
  tmpSP10 = (double*) malloc(sizeof(double) * (1+((K) + (1))-(1)));
  tmpSP11 = (double*) malloc(sizeof(double) * (1+((K) + (1))-(1)));
  vec_var_0 = (double*) malloc(sizeof(double) * (1+((K) + (1))-(1)));
  for (d_13=1; d_13<=D; d_13++) {
    /* Implements direct sampling from the following distribution: */
    /*   Delta(post_theta_{d@13} | +(\sum_{n@4 \in [N_{d@13}]} IDR(z1_{d@13,n@4}, 1, K), \sum_{n@22 \in [N_{d@13}]} IDR(z2_{d@13,n@22}, 1, K)), K) */
    for (dvv_loop_var_1=1; dvv_loop_var_1<=K; dvv_loop_var_1++) {
      tmpSP10[dvv_loop_var_1-1] = 0.0;
    }
    tmpSP10[(K) + (1)-1] = (0.0) * (((1) + (K)) - (1));
    for (n_4=1; n_4<=N[d_13-1]; n_4++) {
      tmpSP10[(K) + (1)-1] += 1.0;
      tmpSP10[z1[d_13-1][n_4-1]-1] += 1.0;
    }
    for (dvv_loop_var_1=1; dvv_loop_var_1<=K; dvv_loop_var_1++) {
      tmpSP11[dvv_loop_var_1-1] = 0.0;
    }
    tmpSP11[(K) + (1)-1] = (0.0) * (((1) + (K)) - (1));
    for (n_22=1; n_22<=N[d_13-1]; n_22++) {
      tmpSP11[(K) + (1)-1] += 1.0;
      tmpSP11[z2[d_13-1][n_22-1]-1] += 1.0;
    }
    sample_Delta(post_theta[d_13-1], add_vec_r_1(vec_var_0, tmpSP10, tmpSP11, 1, K), K);
  }
  free(tmpSP10);
  free(tmpSP11);
  free(vec_var_0);
}

void resample_post_beta1(int D, int K, int* N, int V1, int** arg1, double** post_beta1, int** z1) {
  int k_11;
  double* tmpSP6;
  int d_2;
  int n_110;
  int dvv_loop_var_1;
  tmpSP6 = (double*) malloc(sizeof(double) * (1+((V1) + (1))-(1)));
  for (k_11=1; k_11<=K; k_11++) {
    /* Implements direct sampling from the following distribution: */
    /*   Delta(post_beta1_{k@11} | \sum_{d@2 \in [D]} \sum_{n@110 \in [N_{d@2}]} .*(=(k@11, z1_{d@2,n@110}), IDR(arg1_{d@2,n@110}, 1, V1)), V1) */
    for (dvv_loop_var_1=1; dvv_loop_var_1<=V1; dvv_loop_var_1++) {
      tmpSP6[dvv_loop_var_1-1] = 0.0;
    }
    tmpSP6[(V1) + (1)-1] = (0.0) * (((1) + (V1)) - (1));
    for (d_2=1; d_2<=D; d_2++) {
      for (n_110=1; n_110<=N[d_2-1]; n_110++) {
        tmpSP6[(V1) + (1)-1] += (1.0) * ((((k_11) == (z1[d_2-1][n_110-1])) ? 1 : 0));
        tmpSP6[arg1[d_2-1][n_110-1]-1] += (1.0) * ((((k_11) == (z1[d_2-1][n_110-1])) ? 1 : 0));
      }
    }
    sample_Delta(post_beta1[k_11-1], tmpSP6, V1);
  }
  free(tmpSP6);
}

void resample_post_beta2(int D, int K, int* N, int V2, int** arg2, double** post_beta2, int** z2) {
  int k_12;
  double* tmpSP8;
  int d_3;
  int n_111;
  int dvv_loop_var_1;
  tmpSP8 = (double*) malloc(sizeof(double) * (1+((V2) + (1))-(1)));
  for (k_12=1; k_12<=K; k_12++) {
    /* Implements direct sampling from the following distribution: */
    /*   Delta(post_beta2_{k@12} | \sum_{d@3 \in [D]} \sum_{n@111 \in [N_{d@3}]} .*(=(k@12, z2_{d@3,n@111}), IDR(arg2_{d@3,n@111}, 1, V2)), V2) */
    for (dvv_loop_var_1=1; dvv_loop_var_1<=V2; dvv_loop_var_1++) {
      tmpSP8[dvv_loop_var_1-1] = 0.0;
    }
    tmpSP8[(V2) + (1)-1] = (0.0) * (((1) + (V2)) - (1));
    for (d_3=1; d_3<=D; d_3++) {
      for (n_111=1; n_111<=N[d_3-1]; n_111++) {
        tmpSP8[(V2) + (1)-1] += (1.0) * ((((k_12) == (z2[d_3-1][n_111-1])) ? 1 : 0));
        tmpSP8[arg2[d_3-1][n_111-1]-1] += (1.0) * ((((k_12) == (z2[d_3-1][n_111-1])) ? 1 : 0));
      }
    }
    sample_Delta(post_beta2[k_12-1], tmpSP8, V2);
  }
  free(tmpSP8);
}

double resample_alpha(int D, int K, double alpha, double** post_theta) {
  double tmpSP0;
  int d_0;
  int cgds;
  /* Implements direct sampling from the following distribution: */
  /*   Gam(alpha | 0.1, /(1.0, -(1.0, /(1.0, \sum_{d@0 \in [D]} \sum_{cgds \in [K]} log(.*(/(1.0, sub(.+(alpha, post_theta_{d@0}), +(K, 1))), .+(alpha, post_theta_{d@0,cgds}))))))) */
  tmpSP0 = 0.0;
  for (d_0=1; d_0<=D; d_0++) {
    for (cgds=1; cgds<=K; cgds++) {
      tmpSP0 += log(((1.0) / ((alpha) + (post_theta[d_0-1][(K) + (1)-1]))) * ((alpha) + (post_theta[d_0-1][cgds-1])));
    }
  }
  alpha = sample_Gam(0.1, (1.0) / ((1.0) - ((1.0) / (tmpSP0))));
  return (alpha);
}

double resample_eta(int K, int V1, int V2, double eta, double** post_beta1, double** post_beta2) {
  double tmpSP2;
  int k_10;
  int cgds;
  double tmpSP4;
  int k_1;
  /* Implements direct sampling from the following distribution: */
  /*   Gam(eta | 0.1, /(1.0, -(1.0, /(1.0, +(\sum_{k@10 \in [K]} \sum_{cgds \in [V2]} log(.*(/(1.0, sub(.+(eta, post_beta2_{k@10}), +(V2, 1))), .+(eta, post_beta2_{k@10,cgds}))), \sum_{k@1 \in [K]} \sum_{cgds \in [V1]} log(.*(/(1.0, sub(.+(eta, post_beta1_{k@1}), +(V1, 1))), .+(eta, post_beta1_{k@1,cgds})))))))) */
  tmpSP2 = 0.0;
  for (k_10=1; k_10<=K; k_10++) {
    for (cgds=1; cgds<=V2; cgds++) {
      tmpSP2 += log(((1.0) / ((eta) + (post_beta2[k_10-1][(V2) + (1)-1]))) * ((eta) + (post_beta2[k_10-1][cgds-1])));
    }
  }
  tmpSP4 = 0.0;
  for (k_1=1; k_1<=K; k_1++) {
    for (cgds=1; cgds<=V1; cgds++) {
      tmpSP4 += log(((1.0) / ((eta) + (post_beta1[k_1-1][(V1) + (1)-1]))) * ((eta) + (post_beta1[k_1-1][cgds-1])));
    }
  }
  eta = sample_Gam(0.1, (1.0) / ((1.0) - ((1.0) / ((tmpSP2) + (tmpSP4)))));
  return (eta);
}

void resample_z1(int D, int* N, double alpha, int** arg1, double eta, double** post_beta1, double** post_theta, int** z1, int K, int V1) {
  int d_14;
  int n_23;
  double* tmp_post_z1_2;
  int tmp_idx_z1_2;
  int dvv_loop_var_1;
  tmp_post_z1_2 = (double*) malloc(sizeof(double) * (1+((K) + (1))-(1)));
  for (d_14=1; d_14<=D; d_14++) {
    for (n_23=1; n_23<=N[d_14-1]; n_23++) {
      post_beta1[z1[d_14-1][n_23-1]-1][(V1) + (1)-1] += (0.0) - ((1.0) * ((((z1[d_14-1][n_23-1]) == (z1[d_14-1][n_23-1])) ? 1 : 0)));
      post_beta1[z1[d_14-1][n_23-1]-1][arg1[d_14-1][n_23-1]-1] += (0.0) - ((1.0) * ((((z1[d_14-1][n_23-1]) == (z1[d_14-1][n_23-1])) ? 1 : 0)));
      post_theta[d_14-1][(K) + (1)-1] += (0.0) - (1.0);
      post_theta[d_14-1][z1[d_14-1][n_23-1]-1] += (0.0) - (1.0);
      /* Implements multinomial sampling from the following distribution: */
      /*   (Mult(arg1_{d@14,n@23} | .+(eta, sub(post_beta1, z1_{d@14,n@23}))))(Mult(z1_{d@14,n@23} | .+(alpha, post_theta_{d@14}))) */
      for (dvv_loop_var_1=1; dvv_loop_var_1<=K; dvv_loop_var_1++) {
        tmp_post_z1_2[dvv_loop_var_1-1] = 0.0;
      }
      tmp_post_z1_2[(K) + (1)-1] = (0.0) * (((1) + (K)) - (1));
      for (tmp_idx_z1_2=1; tmp_idx_z1_2<=K; tmp_idx_z1_2++) {
        tmp_post_z1_2[tmp_idx_z1_2-1] = (ldf_Mult_smooth(0, eta, arg1[d_14-1][n_23-1], post_beta1[tmp_idx_z1_2-1], 1, V1)) + (ldf_Mult_smooth(0, alpha, tmp_idx_z1_2, post_theta[d_14-1], 1, K));
      }
      normalizeLog(tmp_post_z1_2, 1, K);
      z1[d_14-1][n_23-1] = sample_Mult(tmp_post_z1_2, 1, K);
      post_theta[d_14-1][(K) + (1)-1] += 1.0;
      post_theta[d_14-1][z1[d_14-1][n_23-1]-1] += 1.0;
      post_beta1[z1[d_14-1][n_23-1]-1][(V1) + (1)-1] += (1.0) * ((((z1[d_14-1][n_23-1]) == (z1[d_14-1][n_23-1])) ? 1 : 0));
      post_beta1[z1[d_14-1][n_23-1]-1][arg1[d_14-1][n_23-1]-1] += (1.0) * ((((z1[d_14-1][n_23-1]) == (z1[d_14-1][n_23-1])) ? 1 : 0));
    }
  }
  free(tmp_post_z1_2);
}

void resample_z2(int D, int* N, double alpha, int** arg2, double eta, double** post_beta2, double** post_theta, int** z2, int K, int V2) {
  int d_15;
  int n_24;
  double* tmp_post_z2_2;
  int tmp_idx_z2_2;
  int dvv_loop_var_1;
  tmp_post_z2_2 = (double*) malloc(sizeof(double) * (1+((K) + (1))-(1)));
  for (d_15=1; d_15<=D; d_15++) {
    for (n_24=1; n_24<=N[d_15-1]; n_24++) {
      post_beta2[z2[d_15-1][n_24-1]-1][(V2) + (1)-1] += (0.0) - ((1.0) * ((((z2[d_15-1][n_24-1]) == (z2[d_15-1][n_24-1])) ? 1 : 0)));
      post_beta2[z2[d_15-1][n_24-1]-1][arg2[d_15-1][n_24-1]-1] += (0.0) - ((1.0) * ((((z2[d_15-1][n_24-1]) == (z2[d_15-1][n_24-1])) ? 1 : 0)));
      post_theta[d_15-1][(K) + (1)-1] += (0.0) - (1.0);
      post_theta[d_15-1][z2[d_15-1][n_24-1]-1] += (0.0) - (1.0);
      /* Implements multinomial sampling from the following distribution: */
      /*   (Mult(arg2_{d@15,n@24} | .+(eta, sub(post_beta2, z2_{d@15,n@24}))))(Mult(z2_{d@15,n@24} | .+(alpha, post_theta_{d@15}))) */
      for (dvv_loop_var_1=1; dvv_loop_var_1<=K; dvv_loop_var_1++) {
        tmp_post_z2_2[dvv_loop_var_1-1] = 0.0;
      }
      tmp_post_z2_2[(K) + (1)-1] = (0.0) * (((1) + (K)) - (1));
      for (tmp_idx_z2_2=1; tmp_idx_z2_2<=K; tmp_idx_z2_2++) {
        tmp_post_z2_2[tmp_idx_z2_2-1] = (ldf_Mult_smooth(0, eta, arg2[d_15-1][n_24-1], post_beta2[tmp_idx_z2_2-1], 1, V2)) + (ldf_Mult_smooth(0, alpha, tmp_idx_z2_2, post_theta[d_15-1], 1, K));
      }
      normalizeLog(tmp_post_z2_2, 1, K);
      z2[d_15-1][n_24-1] = sample_Mult(tmp_post_z2_2, 1, K);
      post_theta[d_15-1][(K) + (1)-1] += 1.0;
      post_theta[d_15-1][z2[d_15-1][n_24-1]-1] += 1.0;
      post_beta2[z2[d_15-1][n_24-1]-1][(V2) + (1)-1] += (1.0) * ((((z2[d_15-1][n_24-1]) == (z2[d_15-1][n_24-1])) ? 1 : 0));
      post_beta2[z2[d_15-1][n_24-1]-1][arg2[d_15-1][n_24-1]-1] += (1.0) * ((((z2[d_15-1][n_24-1]) == (z2[d_15-1][n_24-1])) ? 1 : 0));
    }
  }
  free(tmp_post_z2_2);
}

void resample_arg1(int D, int* N, int** arg1, double eta, double** post_beta1, int** z1, int V1) {
  int d_16;
  int n_25;
  for (d_16=1; d_16<=D; d_16++) {
    for (n_25=1; n_25<=N[d_16-1]; n_25++) {
      post_beta1[z1[d_16-1][n_25-1]-1][(V1) + (1)-1] += (0.0) - ((1.0) * ((((z1[d_16-1][n_25-1]) == (z1[d_16-1][n_25-1])) ? 1 : 0)));
      post_beta1[z1[d_16-1][n_25-1]-1][arg1[d_16-1][n_25-1]-1] += (0.0) - ((1.0) * ((((z1[d_16-1][n_25-1]) == (z1[d_16-1][n_25-1])) ? 1 : 0)));
      /* Implements direct sampling from the following distribution: */
      /*   Mult(arg1_{d@16,n@25} | .+(eta, sub(post_beta1, z1_{d@16,n@25}))) */
      arg1[d_16-1][n_25-1] = sample_Mult_smooth(eta, post_beta1[z1[d_16-1][n_25-1]-1], 1, V1);
      post_beta1[z1[d_16-1][n_25-1]-1][(V1) + (1)-1] += (1.0) * ((((z1[d_16-1][n_25-1]) == (z1[d_16-1][n_25-1])) ? 1 : 0));
      post_beta1[z1[d_16-1][n_25-1]-1][arg1[d_16-1][n_25-1]-1] += (1.0) * ((((z1[d_16-1][n_25-1]) == (z1[d_16-1][n_25-1])) ? 1 : 0));
    }
  }
}

void resample_arg2(int D, int* N, int** arg2, double eta, double** post_beta2, int** z2, int V2) {
  int d_17;
  int n_26;
  for (d_17=1; d_17<=D; d_17++) {
    for (n_26=1; n_26<=N[d_17-1]; n_26++) {
      post_beta2[z2[d_17-1][n_26-1]-1][(V2) + (1)-1] += (0.0) - ((1.0) * ((((z2[d_17-1][n_26-1]) == (z2[d_17-1][n_26-1])) ? 1 : 0)));
      post_beta2[z2[d_17-1][n_26-1]-1][arg2[d_17-1][n_26-1]-1] += (0.0) - ((1.0) * ((((z2[d_17-1][n_26-1]) == (z2[d_17-1][n_26-1])) ? 1 : 0)));
      /* Implements direct sampling from the following distribution: */
      /*   Mult(arg2_{d@17,n@26} | .+(eta, sub(post_beta2, z2_{d@17,n@26}))) */
      arg2[d_17-1][n_26-1] = sample_Mult_smooth(eta, post_beta2[z2[d_17-1][n_26-1]-1], 1, V2);
      post_beta2[z2[d_17-1][n_26-1]-1][(V2) + (1)-1] += (1.0) * ((((z2[d_17-1][n_26-1]) == (z2[d_17-1][n_26-1])) ? 1 : 0));
      post_beta2[z2[d_17-1][n_26-1]-1][arg2[d_17-1][n_26-1]-1] += (1.0) * ((((z2[d_17-1][n_26-1]) == (z2[d_17-1][n_26-1])) ? 1 : 0));
    }
  }
}


/************************* INITIALIZATION *************************/

double initialize_alpha() {
  double alpha;
  alpha = sample_Gam(1.0, 1.0);
  return (alpha);
}

double initialize_eta() {
  double eta;
  eta = sample_Gam(1.0, 1.0);
  return (eta);
}

void initialize_z1(int** z1, int D, int* N, int K) {
  int d_14;
  int n_23;
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=N[dvv_loop_var_1-1]; dvv_loop_var_2++) {
      z1[dvv_loop_var_1-1][dvv_loop_var_2-1] = 0;
    }
    z1[dvv_loop_var_1-1][(N[dvv_loop_var_1-1]) + (1)-1] = (0) * (((1) + (N[dvv_loop_var_1-1])) - (1));
  }
  for (d_14=1; d_14<=D; d_14++) {
    for (n_23=1; n_23<=N[d_14-1]; n_23++) {
      z1[d_14-1][n_23-1] = sample_MultSym(1, K);
    }
  }
}

void initialize_z2(int** z2, int D, int* N, int K) {
  int d_15;
  int n_24;
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=N[dvv_loop_var_1-1]; dvv_loop_var_2++) {
      z2[dvv_loop_var_1-1][dvv_loop_var_2-1] = 0;
    }
    z2[dvv_loop_var_1-1][(N[dvv_loop_var_1-1]) + (1)-1] = (0) * (((1) + (N[dvv_loop_var_1-1])) - (1));
  }
  for (d_15=1; d_15<=D; d_15++) {
    for (n_24=1; n_24<=N[d_15-1]; n_24++) {
      z2[d_15-1][n_24-1] = sample_MultSym(1, K);
    }
  }
}

void initialize_arg1(int** arg1, int D, int* N, int V1) {
  int d_16;
  int n_25;
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=N[dvv_loop_var_1-1]; dvv_loop_var_2++) {
      arg1[dvv_loop_var_1-1][dvv_loop_var_2-1] = 0;
    }
    arg1[dvv_loop_var_1-1][(N[dvv_loop_var_1-1]) + (1)-1] = (0) * (((1) + (N[dvv_loop_var_1-1])) - (1));
  }
  for (d_16=1; d_16<=D; d_16++) {
    for (n_25=1; n_25<=N[d_16-1]; n_25++) {
      arg1[d_16-1][n_25-1] = sample_MultSym(1, V1);
    }
  }
}

void initialize_arg2(int** arg2, int D, int* N, int V2) {
  int d_17;
  int n_26;
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=N[dvv_loop_var_1-1]; dvv_loop_var_2++) {
      arg2[dvv_loop_var_1-1][dvv_loop_var_2-1] = 0;
    }
    arg2[dvv_loop_var_1-1][(N[dvv_loop_var_1-1]) + (1)-1] = (0) * (((1) + (N[dvv_loop_var_1-1])) - (1));
  }
  for (d_17=1; d_17<=D; d_17++) {
    for (n_26=1; n_26<=N[d_17-1]; n_26++) {
      arg2[d_17-1][n_26-1] = sample_MultSym(1, V2);
    }
  }
}

void initialize_post_theta(double** post_theta, int D, int K, int* N, int** z1, int** z2) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=K; dvv_loop_var_2++) {
      post_theta[dvv_loop_var_1-1][dvv_loop_var_2-1] = 0.0;
    }
    post_theta[dvv_loop_var_1-1][(K) + (1)-1] = (0.0) * (((1) + (K)) - (1));
  }
  resample_post_theta(D, K, N, post_theta, z1, z2);
}

void initialize_post_beta1(double** post_beta1, int D, int K, int* N, int V1, int** arg1, int** z1) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  for (dvv_loop_var_1=1; dvv_loop_var_1<=K; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=V1; dvv_loop_var_2++) {
      post_beta1[dvv_loop_var_1-1][dvv_loop_var_2-1] = 0.0;
    }
    post_beta1[dvv_loop_var_1-1][(V1) + (1)-1] = (0.0) * (((1) + (V1)) - (1));
  }
  resample_post_beta1(D, K, N, V1, arg1, post_beta1, z1);
}

void initialize_post_beta2(double** post_beta2, int D, int K, int* N, int V2, int** arg2, int** z2) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  for (dvv_loop_var_1=1; dvv_loop_var_1<=K; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=V2; dvv_loop_var_2++) {
      post_beta2[dvv_loop_var_1-1][dvv_loop_var_2-1] = 0.0;
    }
    post_beta2[dvv_loop_var_1-1][(V2) + (1)-1] = (0.0) * (((1) + (V2)) - (1));
  }
  resample_post_beta2(D, K, N, V2, arg2, post_beta2, z2);
}


/**************************** DUMPING *****************************/

void dump_alpha(double alpha) {
  printf("alpha = ");
  printf("%g", alpha);
  printf("\n"); fflush(stdout);
}

void dump_eta(double eta) {
  printf("eta = ");
  printf("%g", eta);
  printf("\n"); fflush(stdout);
}

void dump_beta1(int K, int V1, double** beta1) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  printf("beta1 = ");
  for (dvv_loop_var_1=1; dvv_loop_var_1<=K; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=V1; dvv_loop_var_2++) {
      printf("%g", beta1[dvv_loop_var_1-1][dvv_loop_var_2-1]);
      printf(" ");
    }
    printf(" ; ");
  }
  printf("\n"); fflush(stdout);
}

void dump_beta2(int K, int V2, double** beta2) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  printf("beta2 = ");
  for (dvv_loop_var_1=1; dvv_loop_var_1<=K; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=V2; dvv_loop_var_2++) {
      printf("%g", beta2[dvv_loop_var_1-1][dvv_loop_var_2-1]);
      printf(" ");
    }
    printf(" ; ");
  }
  printf("\n"); fflush(stdout);
}

void dump_theta(int D, int K, double** theta) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  printf("theta = ");
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=K; dvv_loop_var_2++) {
      printf("%g", theta[dvv_loop_var_1-1][dvv_loop_var_2-1]);
      printf(" ");
    }
    printf(" ; ");
  }
  printf("\n"); fflush(stdout);
}

void dump_z1(int D, int* N, int** z1) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  printf("z1 = ");
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=N[dvv_loop_var_1-1]; dvv_loop_var_2++) {
      printf("%d", z1[dvv_loop_var_1-1][dvv_loop_var_2-1]);
      printf(" ");
    }
    printf(" ; ");
  }
  printf("\n"); fflush(stdout);
}

void dump_z2(int D, int* N, int** z2) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  printf("z2 = ");
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=N[dvv_loop_var_1-1]; dvv_loop_var_2++) {
      printf("%d", z2[dvv_loop_var_1-1][dvv_loop_var_2-1]);
      printf(" ");
    }
    printf(" ; ");
  }
  printf("\n"); fflush(stdout);
}

void dump_arg1(int D, int* N, int** arg1) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  printf("arg1 = ");
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=N[dvv_loop_var_1-1]; dvv_loop_var_2++) {
      printf("%d", arg1[dvv_loop_var_1-1][dvv_loop_var_2-1]);
      printf(" ");
    }
    printf(" ; ");
  }
  printf("\n"); fflush(stdout);
}

void dump_arg2(int D, int* N, int** arg2) {
  int dvv_loop_var_1;
  int dvv_loop_var_2;
  printf("arg2 = ");
  for (dvv_loop_var_1=1; dvv_loop_var_1<=D; dvv_loop_var_1++) {
    for (dvv_loop_var_2=1; dvv_loop_var_2<=N[dvv_loop_var_1-1]; dvv_loop_var_2++) {
      printf("%d", arg2[dvv_loop_var_1-1][dvv_loop_var_2-1]);
      printf(" ");
    }
    printf(" ; ");
  }
  printf("\n"); fflush(stdout);
}


/*************************** LIKELIHOOD ***************************/

double compute_log_posterior(int D, int K, int* N, int V1, int V2, double alpha, int** arg1, int** arg2, double** beta1, double** beta2, double eta, double** theta, int** z1, int** z2) {
  double ldfP5_0;
  double ldfP5_1;
  int n_23;
  int d_14;
  double ldfP6_0;
  double ldfP6_1;
  int n_24;
  int d_15;
  double ldfP7_0;
  double ldfP7_1;
  int n_25;
  int d_16;
  double ldfP8_0;
  double ldfP8_1;
  int n_26;
  int d_17;
  ldfP5_0 = 0.0;
  for (d_14=1; d_14<=D; d_14++) {
    ldfP5_1 = 0.0;
    for (n_23=1; n_23<=N[d_14-1]; n_23++) {
      ldfP5_1 += ldf_Mult(1, z1[d_14-1][n_23-1], theta[d_14-1], 1, K);
    }
    ldfP5_0 += ldfP5_1;
  }
  ldfP6_0 = 0.0;
  for (d_15=1; d_15<=D; d_15++) {
    ldfP6_1 = 0.0;
    for (n_24=1; n_24<=N[d_15-1]; n_24++) {
      ldfP6_1 += ldf_Mult(1, z2[d_15-1][n_24-1], theta[d_15-1], 1, K);
    }
    ldfP6_0 += ldfP6_1;
  }
  ldfP7_0 = 0.0;
  for (d_16=1; d_16<=D; d_16++) {
    ldfP7_1 = 0.0;
    for (n_25=1; n_25<=N[d_16-1]; n_25++) {
      ldfP7_1 += ldf_Mult(1, arg1[d_16-1][n_25-1], beta1[z1[d_16-1][n_25-1]-1], 1, V1);
    }
    ldfP7_0 += ldfP7_1;
  }
  ldfP8_0 = 0.0;
  for (d_17=1; d_17<=D; d_17++) {
    ldfP8_1 = 0.0;
    for (n_26=1; n_26<=N[d_17-1]; n_26++) {
      ldfP8_1 += ldf_Mult(1, arg2[d_17-1][n_26-1], beta2[z2[d_17-1][n_26-1]-1], 1, V2);
    }
    ldfP8_0 += ldfP8_1;
  }
  return ((ldf_Gam(1, alpha, 0.1, 1)) + ((ldf_Gam(1, eta, 0.1, 1)) + ((0.0) + ((0.0) + ((0.0) + ((ldfP5_0) + ((ldfP6_0) + ((ldfP7_0) + (ldfP8_0)))))))));
}

/****************************** MAIN ******************************/

int main(int ARGC, char *ARGV[]) {
  double loglik,bestloglik;
  int iter;
  int D;
  int K;
  int* N;
  int V1;
  int V2;
  double alpha;
  int** arg1;
  int** arg2;
  double eta;
  double** post_beta1;
  double** post_beta2;
  double** post_theta;
  int** z1;
  int** z2;
  int malloc_dim_1;

  int nIter;
  int dumpInterval;

  if(ARGC != 5) {
    fprintf(stderr, "usage: %s <arg1> <arg2> <nIter> <dumpInterval>\n", ARGV[0]);
    exit(1);
  }

  nIter = atoi(ARGV[3]);
  dumpInterval = atoi(ARGV[4]);

  fprintf(stderr, "-- This program was automatically generated using HBC (v 0.7 beta) from link_lda_300.hier\n--     see http://hal3.name/HBC for more information\n");
  fflush(stderr);
  setall(time(0),time(0));   /* initialize random number generator */


  /* variables defined with --define */
  alpha = 0.1;
  eta = 0.1;
  K = 300;

  fprintf(stderr, "Loading data...\n");
  fflush(stderr);
  /* variables defined with --loadD */
  arg1 = load_discrete2(ARGV[1], &D, &N, &V1);
  arg2 = load_discrete2(ARGV[2], &D, &N, &V2);

  /* variables defined with --loadM or --loadMI */

  fprintf(stderr, "Allocating memory...\n");
  fflush(stderr);
  post_beta1 = (double**) malloc(sizeof(double*) * (1+(K)-(1)));
  for (malloc_dim_1=1; malloc_dim_1<=K; malloc_dim_1++) {
    post_beta1[malloc_dim_1-1] = (double*) malloc(sizeof(double) * (1+((V1) + (1))-(1)));
  }

  post_beta2 = (double**) malloc(sizeof(double*) * (1+(K)-(1)));
  for (malloc_dim_1=1; malloc_dim_1<=K; malloc_dim_1++) {
    post_beta2[malloc_dim_1-1] = (double*) malloc(sizeof(double) * (1+((V2) + (1))-(1)));
  }

  post_theta = (double**) malloc(sizeof(double*) * (1+(D)-(1)));
  for (malloc_dim_1=1; malloc_dim_1<=D; malloc_dim_1++) {
    post_theta[malloc_dim_1-1] = (double*) malloc(sizeof(double) * (1+((K) + (1))-(1)));
  }

  z1 = (int**) malloc(sizeof(int*) * (1+(D)-(1)));
  for (malloc_dim_1=1; malloc_dim_1<=D; malloc_dim_1++) {
    z1[malloc_dim_1-1] = (int*) malloc(sizeof(int) * (1+((N[malloc_dim_1-1]) + (1))-(1)));
  }

  z2 = (int**) malloc(sizeof(int*) * (1+(D)-(1)));
  for (malloc_dim_1=1; malloc_dim_1<=D; malloc_dim_1++) {
    z2[malloc_dim_1-1] = (int*) malloc(sizeof(int) * (1+((N[malloc_dim_1-1]) + (1))-(1)));
  }


  fprintf(stderr, "Initializing variables...\n");
  fflush(stderr);
  initialize_z1(z1, D, N, K);
  initialize_z2(z2, D, N, K);
  initialize_post_theta(post_theta, D, K, N, z1, z2);
  initialize_post_beta1(post_beta1, D, K, N, V1, arg1, z1);
  initialize_post_beta2(post_beta2, D, K, N, V2, arg2, z2);

  for (iter=1; iter<=nIter; iter++) {
    fprintf(stderr, "iter %d", iter);
    fflush(stderr);
    resample_z1(D, N, alpha, arg1, eta, post_beta1, post_theta, z1, K, V1);
    resample_z2(D, N, alpha, arg2, eta, post_beta2, post_theta, z2, K, V2);

    loglik = compute_log_posterior(D, K, N, V1, V2, alpha, arg1, arg2, post_beta1, post_beta2, eta, post_theta, z1, z2);
    fprintf(stderr, "\t%g", loglik);
    if ((iter==1)||(loglik>bestloglik)) {
      bestloglik = loglik;
      fprintf(stderr, " *");
    }
    fprintf(stderr, "\n");
    fflush(stderr);
    if ((iter % dumpInterval) == 0) {
      printf("ll = %g\n", loglik);
      dump_z1(D, N, z1);
      dump_z2(D, N, z2);
    }
  }

  for (malloc_dim_1=1; malloc_dim_1<=D; malloc_dim_1++) {
    free(z2[malloc_dim_1-1]);
  }
  free(z2);

  for (malloc_dim_1=1; malloc_dim_1<=D; malloc_dim_1++) {
    free(z1[malloc_dim_1-1]);
  }
  free(z1);

  for (malloc_dim_1=1; malloc_dim_1<=D; malloc_dim_1++) {
    free(post_theta[malloc_dim_1-1]);
  }
  free(post_theta);

  for (malloc_dim_1=1; malloc_dim_1<=K; malloc_dim_1++) {
    free(post_beta2[malloc_dim_1-1]);
  }
  free(post_beta2);

  for (malloc_dim_1=1; malloc_dim_1<=K; malloc_dim_1++) {
    free(post_beta1[malloc_dim_1-1]);
  }
  free(post_beta1);

  for (malloc_dim_1=1; malloc_dim_1<=D; malloc_dim_1++) {
    free(arg2[malloc_dim_1-1]);
  }
  free(arg2);

  for (malloc_dim_1=1; malloc_dim_1<=D; malloc_dim_1++) {
    free(arg1[malloc_dim_1-1]);
  }
  free(arg1);

  free(N);


  return 0;
}
