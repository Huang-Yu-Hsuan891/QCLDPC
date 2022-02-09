#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

unsigned long long SEED = 3122891; // have value
//unsigned long long SEED = 3881111387891;
unsigned long long RANV;
int RANI = 0;

double Ranq1();
void normal(double sig, double *n1, double *n2);
double CHK(double L1, double L2);
double table[8] = {0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, 0};

int main(){
    // for declaration
    int i, j, k, m;                  //for counting
    int n, rc;              // n is column and rc is row
    int dv,dc;              // dv: column have #1 and dc: row have #1
    int a;                  // no need
    int **L = NULL;                 // check node connect 6 variable nodes
    int Llenrow = rc;        // rc = 408
    int Llencolumn = dc;     // dc = 6
    int **M = NULL;
    int Mlenrow = n;       // n = 816
    int Mlencolumn = dv;   //dv = 3
    int *codarray;          //codeword, Dynamic memory allocation, length = codarraylen
    int codarraylen = n;
    int *codearray;         //0->1; 1->-1,  Dynamic memory allocation, length = codearraylen
    int codearraylen = n;
    double *outp;           // codeword + noise, Dynamic memory allocation, length = outparray
    int outparray = n;
    int *output;            // result of interative algorithm decoding, Dynamic memory allocation, length = outputarray
    int outputarray = n;
    double *Lj;             // LLR
    int Ljlen = n;
    double **qij = NULL;    // from down to top
    int qijrow = dv;        // qijrow = dv = 3
    int qijcolumn = n;      // qijcolumn = n = 816
    double **uij = NULL;    // from top to down
    int uijrow = rc;        // uijrow = rc = 408
    int uijcolumn = dc;     // uijcloumn = dc 3
    double tempqij[5]; 
    double tempuij;
    double temp1uij[3];
    double temp1qij; 
    int *comput;
    int computlen = n;      // computlen = 816
    int *comput1;
    int comput1len = rc;    // comput1len = 408
    int *comput2;
    int comput2len = rc;
    double *qj;
    int qjlen = n;
    int *checkbit;
    int checkbitlen = rc;
    int krc;
    int **G;
    int Glenrow = krc;
    int Glencolumn = n;
    int *u;
    int ulen = krc;

    int step;               // test 6 EbN0(dB) result
    int s = 0;              // receive 100 error block
    int num = 0;            // do compute block
    int error;
    int totalerror=0;
    int restart = 0;
    double sigma;
    double ebn0;
    double ebn0s[6];
    double bers[6];
    double berscompare[6];
    double x, y;
    int valL;
    int valL2;
    double app;
    double app1;
    int stp;
    // declaration end

    for (i = 0; i < 6; i++) bers[i] = 0;   
    ebn0s[0] = 1.5;
    ebn0s[1] = 1.75;
    ebn0s[2] = 2.0;
    berscompare[0] = 4 * pow(10,-3);
    berscompare[1] = 3.47 * pow(10,-4);
    berscompare[2] = 1.29 * pow(10,-5); 

    // open file
    FILE *fpr;
    fpr=fopen("paritycheckmatrix.txt","r");
    fscanf(fpr,"%d",&n);
    fscanf(fpr,"%d",&rc);
    printf("column = %d\n", n);
    printf("row = %d\n", rc);
    fscanf(fpr,"%d",&dv);
    fscanf(fpr,"%d",&dc);
    printf("dv = %d\n", dv);
    printf("dc = %d\n", dc); 

    Llenrow = rc;
    Llencolumn = dc;
    Mlenrow = n; 
    Mlencolumn = dv;

    M = (int **)malloc(Mlenrow * sizeof(int *));
    for (i = 0; i < Mlenrow; i++) M[i] = (int *)malloc(Mlencolumn * sizeof(int));
    L = (int **)malloc(Llenrow * sizeof(int *));
    for (i = 0; i < Llenrow; i++) L[i] = (int *)malloc(Llencolumn * sizeof(int));

    for (j = 0; j < Mlenrow; j++) {
        for (i = 0; i < Mlencolumn; i++) {
            fscanf(fpr,"%d",&M[j][i]);
        }
    }
    int temp;
    for (j = 0; j < n; j++) {
        for (i = 0; i < dv; i++) {
            for (m = i; m < dv; m++) {
                if (M[j][m] < M[j][i]) {
                    temp = M[j][m];
                    M[j][m] = M[j][i];
                    M[j][i] = temp;
                }
            }
        }
    }
    for (i = 0; i < Llenrow; i++) {
        for (j = 0; j < dc; j++) {
            fscanf(fpr,"%d",&L[i][j]);
        }
    }
    for (i = 0; i < rc; i++) {
        for (j = 0; j < dc; j++) {
            for (m = j; m < dc; m++) {
                if (L[i][m] < L[i][j]) {
                    temp = L[i][m];
                    L[i][m] =L[i][j];
                    L[i][j] = temp;
                }
            }
        }
    }
    fclose(fpr);
    // close file

    FILE *fpr1;
    fpr1=fopen("generator1.txt","r");
    fscanf(fpr1,"%d",&krc);
    fscanf(fpr1,"%d",&n);
    printf("column = %d\n", n);
    printf("row = %d\n", krc);
    Glenrow = krc;
    Glencolumn = n;
    G = (int **)malloc(Glenrow * sizeof(int *));
    for (i = 0; i < Glenrow; i++) G[i] = (int *)malloc(Glencolumn * sizeof(int));
    for (i = 0; i < Glenrow; i++) {
        for (j = 0; j < Glencolumn; j++) {
            fscanf(fpr1,"%d",&G[i][j]);
            if (G[i][j] == 1) printf("%d ", j);
        }
        printf("\n");
    }
    // for CODE part
    fclose(fpr1);
    //printf("yes\n");
    codarraylen = n;
    codarray = (int *)malloc(codarraylen * sizeof(int));
    if( codarray == NULL ) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    codearraylen = n;
    codearray = (int *)malloc( codearraylen * sizeof(int) );
    if( codearray == NULL ) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    outparray = n;
    outp = (double *)malloc( outparray * sizeof(double) );
    if (outp == NULL) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    outputarray = n;
    output = (int *)malloc( outputarray * sizeof(int) );
    if (output == NULL) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    Ljlen = n;
    Lj = (double *)malloc(Ljlen * sizeof(double));
    if (Lj == NULL) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    qijrow = dv;        // qijrow = dv = 3
    qijcolumn = n;      // qijcolumn = n = 816
    qij = (double **)malloc(qijrow * sizeof(double *));
    for (i = 0; i < qijrow; i++) qij[i] = (double *)malloc(qijcolumn * sizeof(double));
    
    uijrow = rc;        // uijrow = rc = 408
    uijcolumn = dc;     // uijcolumn = dc = 6 
    uij = (double **)malloc(uijrow * sizeof(double *));
    for (i = 0; i < uijrow; i++) uij[i] = (double *)malloc(uijcolumn * sizeof(double));

    computlen = n;
    comput = (int *)malloc(computlen * sizeof(int));
    comput1len = rc;
    comput1 = (int *)malloc(comput1len * sizeof(int));
    comput2len = rc;
    comput2 = (int *)malloc(comput2len * sizeof(int));

    qjlen = n;
    qj = (double *)malloc(qjlen * sizeof(double));
    checkbitlen = rc;
    checkbit = (int *)malloc(checkbitlen * sizeof(int));

    ulen = krc;
    u = (int *)malloc(ulen * sizeof(int));
    
    //printf("yes\n");
    for (step = 0; step < 3; step++) {
        s = 0;
        num = 0;
        totalerror = 0;
        while (s < 100/*s < 100*/) {
            //num++;                                  // compute the number of transmit block 
             // initial
            for (i = 0; i < codarraylen; i++) {
                codarray[i] = 0;
            }
            //printf("yes\n");
            if (num == 0) {
                u[0] = 1;
                u[1] = 0;
                u[2] = 0;
                u[3] = 0;
                u[4] = 0;
                u[5] = 0;
                for (i = 0; i < krc - 6; i++) u[i + 6] = (u[i + 1] + u[i]) % 2;
            }
            else {
                u[0] = (u[krc-5] + u[krc-6]) % 2;
                u[1] = (u[krc-4] + u[krc-5]) % 2;
                u[2] = (u[krc-3] + u[krc-4]) % 2;
                u[3] = (u[krc-2] + u[krc-3]) % 2;
                u[4] = (u[krc-1] + u[krc-2]) % 2;
                u[5] = (u[krc-0] + u[krc-1]) % 2;
                for (i = 0; i < krc - 6; i++) u[i + 6] = (u[i + 1] + u[i]) % 2;
            } 
            //printf("yes\n");
            num++;
            for (i = 0; i < ulen; i++) {
                if (u[i] == 1) {
                    for (j = 0; j < n; j++) codarray[j] = (codarray[j] + G[i][j]) % 2;
                }
            }
            //printf("yes\n");
            //input to AWGN channel normalized to +-1
            for (i = 0; i < codearraylen; i++) {
                if (codarray[i] == 0) codearray[i] = 1;
                else codearray[i] = -1;
            }
            //printf("yes\n");
            ebn0 = ebn0s[step];
            sigma = sqrt(1.0 / (pow(10, ebn0/10)));
            for(i = 0; i < rc; i++) {
                normal(sigma, &x, &y);
                outp[2 * i] = codearray[2 * i] + x;
                outp[2 * i + 1] = codearray[2 * i + 1] + y;
            }
            //printf("yes\n");
            ebn0 = pow(10, ebn0/10);
            for(i = 0; i < Ljlen; i++) {
                Lj[i] =4 * 0.5 * ebn0 * outp[i];     //  0.5 * 1.2544 = Es/N0
            }

            // the interative decoding algotrithm
            for (j = 0; j < qijcolumn; j++) {                               // initialization
                for (i = 0; i < qijrow; i++) {
                        qij[i][j] = Lj[j];
                }  
            }
            //printf("yes\n");
            for (k = 0; k < 50/*k < 100*/ && restart != rc; k++) {         // message passing, for predetermined threshold = 100
                restart = 0;    
                for (i = 0; i < 5; i++) {                          // bottom-up
                    tempqij[i] = 0.0;
                }
                for (i = 0; i < computlen; i++) {
                    comput[i] = 0;
                }
                for (i = 0; i < rc; i++) {
                    for (j = 0; j < 6; j++) {
                        for (m = 0; m < 5; m++) {
                            if (m < j) {
                                valL = L[i][m]-1;
                                tempqij[m] = qij[comput[valL]][valL];
                            } 
                            else if (m >= j) {
                                valL = L[i][m+1]-1;
                                tempqij[m] = qij[comput[valL]][valL];
                            }
                        }
                        tempuij = tempqij[0];
                        for(m = 1; m < 5; m++) {
                            tempuij = CHK(tempuij, tempqij[m]);
                        }
                        uij[i][j] = tempuij;
                    }
                    for (m = 0; m < 6; m++) {
                            comput[L[i][m] - 1] += 1;
                    }
                }

                // top-down
                for(i = 0; i < 3; i++) {
                    temp1uij[i] = 0.0;
                }
                for (i = 0; i < comput1len; i++) comput1[i] = 0;
                for (j = 0; j < n; j++) {
                    for (i = 0; i < 3; i++) {
                        for (m = 0; m < 2; m++) {
                            if (m < i) { 
                                valL = M[j][m] - 1;
                                temp1uij[m] = uij[valL][comput1[valL]]; 
                            }
                            else if (m >= i) {
                                valL = M[j][m + 1] - 1;
                                temp1uij[m] = uij[valL][comput1[valL]];
                            }
                        }
                        temp1uij[2] = Lj[j];
                        qij[i][j] = temp1uij[0] + temp1uij[1] + temp1uij[2];
                    }
                    for (m = 0; m < 3; m++) {
                        comput1[M[j][m] - 1] += 1;
                    }
                }

                // decision
                for (i = 0; i < comput2len; i++) comput2[i] = 0;
                for (j = 0; j < n; j++) {
                    qj[j] = Lj[j];
                    for (i = 0; i < 3; i++) {
                        valL = M[j][i] - 1;
                        qj[j] += uij[valL][comput2[valL]];
                    }
                    if (qj[j] >= 0) output[j] = 0;
                    else if (qj[j] < 0) output[j] = 1;
                    for (i = 0; i < 3; i++) {
                        comput2[M[j][i] - 1] += 1;
                    }
                }

                // to check Hx=0     
                for (i = 0; i < rc; i++) {
                    checkbit[i] = 0;
                    for (j = 0; j < 6; j++) {
                        checkbit[i] += output[L[i][j] - 1];
                    }
                    checkbit[i] = checkbit[i] % 2;
                }

                for (i = 0; i < rc; i++) {
                    if (checkbit[i] == 0) restart += 1; // restart = 408 is success
                }
                stp = 0;
                if (k == 99 && restart != rc) {
                    stp = 1;
                    s++;
                }
            }
            /*if(k == 100) */printf("s = %d; k[%d] = %d\n", s, num, k);
            error = 0;
            for(i = 0; i < n; i++) {
                if (output[i] != codarray[i]) {
                    error += 1;
                }
            }
            if (error != 0 && stp == 0) s++;
            restart = 0;
            if(error != 0) printf("error = %d\n", error);       
            totalerror += error;
        }
        double ber;
        ber = (double)totalerror / (num * n);
        printf("totalerror = %d\n", totalerror);
        printf("BER = %g\n", ber);
        bers[step] = ber;
    }

    for(step = 0; step < 3; step++) {
        printf("enb0s[%d] = %g\n",step, ebn0s[step]);
        printf("bers[%d] = %g\n",step, bers[step]);
        printf("berscompare[%d] = %g\n\n", step, berscompare[step]);
    }
    // CODE  end

    // open write file
    FILE *outfp;
    
    outfp = fopen("qc_LDPC_result.txt","w");
    for (i = 0; i < 3; i++) {
         fprintf(outfp,"%g ",ebn0s[i]);
         fprintf(outfp,"%g ",berscompare[i]);
         fprintf(outfp,"%g ",bers[i]);
         fprintf(outfp,"\n");
    }
    fclose(outfp);
    // close write file

    // free dynametic
    free(codarray);
    free(codearray);
    free(outp);
    free(output);
    for (i = 0; i < Llenrow; i++) free(L[i]);
    free(L);
    for (i = 0; i < Mlenrow; i++) free(M[i]);
    free(M);
    for (i = 0; i < qijrow; i++) free(qij[i]);
    free(qij);
    for (i = 0; i < uijrow; i++) free(uij[i]);
    free(uij);
    free(Lj);
    free(comput);
    free(comput1);
    free(comput2);
    for (i = 0; i < krc; i++) free(G[i]);
    free(G);
    free(u);
    return 0;
}
void normal(double sig, double *n1, double *n2)
{   
    double x1,x2;
    double s;
    //printf("sigma = %g\n", sig);
    do{
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while (s >= 1.0);
    *n1 = sig * x1 * sqrt((-2.0 * log(s))/ s);
    *n2 = sig * x2 * sqrt((-2.0 * log(s))/ s);
    
}

double Ranq1() {
    if ( RANI == 0 ){
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV = RANV * 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;

    return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}

double CHK(double L1, double L2) 
{   
    double sgn1, sgn2, min;

    if(L1>0) sgn1 = 1;
    else if(L1 == 0) sgn1 = 0;
    else sgn1 = -1;
    if(L2>0) sgn2 = 1;
    else if(L2 == 0) sgn2 = 0;
    else sgn2 = -1;
    if(fabs(L1) >= fabs(L2)) min = fabs(L2);
    else min = fabs(L1);
    double temp1,temp2;
    double ope1,ope2;
    temp1 = L1+L2;
    temp2 = L1-L2;
    ope1 = fabs(temp1);
    ope2 = fabs(temp2);
    double lnup, lndown;
    if (ope1 >= 0 && ope1 < 0.196) lnup = table[0];
    else if (ope1 >= 0.196 && ope1 < 0.433) lnup = table[1];
    else if (ope1 >= 0.433 && ope1 < 0.71) lnup = table[2];
    else if (ope1 >= 0.71 && ope1 < 1.05) lnup = table[3];
    else if (ope1 >= 1.05 && ope1 < 1.508) lnup = table[4];
    else if (ope1 >= 1.508 && ope1 < 2.252) lnup = table[5];
    else if (ope1 >= 2.252 && ope1 < 4.5) lnup = table[6];
    else if (ope1 >= 4.5) lnup = table[7];

    if (ope2 >= 0 && ope2 < 0.196) lndown = table[0];
    else if (ope2 >= 0.196 && ope2 < 0.433) lndown = table[1];
    else if (ope2 >= 0.433 && ope2 < 0.71) lndown = table[2];
    else if (ope2 >= 0.71 && ope2 < 1.05) lndown = table[3];
    else if (ope2 >= 1.05 && ope2 < 1.508) lndown = table[4];
    else if (ope2 >= 1.508 && ope2 < 2.252) lndown = table[5];
    else if (ope2 >= 2.252 && ope2 < 4.5) lndown = table[6];
    else if (ope2 >= 4.5) lndown = table[7];

    return sgn1 * sgn2 * min +lnup-lndown /*answer*/;
}