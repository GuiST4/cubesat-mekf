#ifndef MEKF_WB_H
#define MEKF_WB_H

typedef struct
{
    float q[4];
    float b[3];
    float P[36];
    float V[36];
    float W[81];
    float dt;
} MEKF_State;

void expq(const float *phi, float *output);

void Rq(const float *q, float *output);
void Lq(const float *q, float *output);
void Gq(const float *q, float *output);

float norm(const float *v);

void mekf_wb(const MEKF_State *input, MEKF_State *output, const float *omega, const float *Br, const float *Nr);

void Mprod(const float *A, const float *B, const int rowsA, const int colsA, const int colsB, const float s, float *C);

void Madd(const float *A, const float *B, const int rows, const int cols, float *C);

void Msub(const float *A, const float *B, const int rows, const int cols, float *C);

void Mtranspose(const float *A, const int rows, const int cols, float *B);

void AkTransform(const float *Ak11, const float *Ak12, float *output);

void Qk1Transform(const float *Q, float *Qfull);

void Ck1Transform(const float *Ck_11, const float *Ck_21, const float *Ck_31, float *Ck1);

void LUSolve9x9(const float *A, const float *B, int nrhs, float *X);

#endif