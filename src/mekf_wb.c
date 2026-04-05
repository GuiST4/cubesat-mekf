#include "mekf_wb.h"
#include <math.h>
#include <string.h>

float norm(const float *v)
{
    return sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void expq(const float *phi, float *output)
{
    float norm_phi = 0;
    norm_phi = norm(phi);

    if (norm_phi < 1.0e-6f)
    {
        output[0] = 1.0f;
        output[1] = 0.0f;
        output[2] = 0.0f;
        output[3] = 0.0f;
    }
    else
    {
        output[0] = cosf(norm_phi);

        float k = sinf(norm_phi) / norm_phi;
        output[1] = phi[0] * k;
        output[2] = phi[1] * k;
        output[3] = phi[2] * k;
    }
}

void Rq(const float *q, float *output)
{
    // Row 1
    output[0] = q[0];
    output[1] = -q[1];
    output[2] = -q[2];
    output[3] = -q[3];

    // Row 2
    output[4] = q[1];
    output[5] = q[0];
    output[6] = q[3];
    output[7] = -q[2];

    // Row 3
    output[8] = q[2];
    output[9] = -q[3];
    output[10] = q[0];
    output[11] = q[1];

    // Row 4
    output[12] = q[3];
    output[13] = q[2];
    output[14] = -q[1];
    output[15] = q[0];
}

void Lq(const float *q, float *output)
{
    // Row 1
    output[0] = q[0];
    output[1] = -q[1];
    output[2] = -q[2];
    output[3] = -q[3];

    // Row 2
    output[4] = q[1];
    output[5] = q[0];
    output[6] = -q[3];
    output[7] = q[2];

    // Row 3
    output[8] = q[2];
    output[9] = q[3];
    output[10] = q[0];
    output[11] = -q[1];

    // Row 4
    output[12] = q[3];
    output[13] = -q[2];
    output[14] = q[1];
    output[15] = q[0];
}

void Gq(const float *q, float *output)
{
    // Row 1
    output[0] = -q[1];
    output[1] = -q[2];
    output[2] = -q[3];

    // Row 2
    output[3] = q[0];
    output[4] = -q[3];
    output[5] = q[2];

    // Row 3
    output[6] = q[3];
    output[7] = q[0];
    output[8] = -q[1];

    // Row 4
    output[9] = -q[2];
    output[10] = q[1];
    output[11] = q[0];
}

void Mprod(const float *A, const float *B, const int rowsA, const int colsA, const int colsB, const float s, float *C)
{
    for (int i = 0; i < rowsA; i++)
    {
        for (int j = 0; j < colsB; j++)
        {
            C[i * colsB + j] = 0.0f;
            for (int k = 0; k < colsA; k++)
            {
                C[i * colsB + j] += A[i * colsA + k] * B[k * colsB + j];
            }
            C[i * colsB + j] *= s;
        }
    }
}

void Madd(const float *A, const float *B, const int rows, const int cols, float *C)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            C[i * cols + j] = A[i * cols + j] + B[i * cols + j];
        }
    }
}

void Mtranspose(const float *A, const int rows, const int cols, float *B)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            B[j * rows + i] = A[i * cols + j];
        }
    }
}

void Msub(const float *A, const float *B, const int rows, const int cols, float *C)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            C[i * cols + j] = A[i * cols + j] - B[i * cols + j];
}

void AkTransform(const float *Ak11, const float *Ak12, float *Ak)
{
    for (int i = 0; i < 36; i++)
        Ak[i] = 0.0f;

    // Ak11
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Ak[i * 6 + j] = Ak11[i * 3 + j];

    // Ak12
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Ak[i * 6 + (j + 3)] = Ak12[i * 3 + j];

    Ak[21] = 1.0f;
    Ak[28] = 1.0f;
    Ak[35] = 1.0f;
}

void Qk1Transform(const float *Q, float *Qfull)
{
    for (int i = 0; i < 81; i++)
        Qfull[i] = 0.0f;

    // Top-left
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Qfull[i * 9 + j] = Q[i * 3 + j];

    // Middle
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Qfull[(i + 3) * 9 + (j + 3)] = Q[i * 3 + j];

    // Bottom-right
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Qfull[(i + 6) * 9 + (j + 6)] = Q[i * 3 + j];
}

void Ck1Transform(const float *Ck_11, const float *Ck_21, const float *Ck_31, float *Ck1)
{
    for (int i = 0; i < 54; i++)
        Ck1[i] = 0.0f;

    // Top-left Ck_11
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Ck1[i * 6 + j] = Ck_11[i * 3 + j];

    // Mid-left Ck_21
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Ck1[(i + 3) * 6 + j] = Ck_21[i * 3 + j];

    // Bottom-left Ck_21
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Ck1[(i + 6) * 6 + j] = Ck_31[i * 3 + j];
}

void LUSolve9x9(const float *A, const float *B, int nrhs, float *X)
{
    float LU[9][9];
    int piv[9];

    for (int i = 0; i < 9; i++)
    {
        piv[i] = i;
        for (int j = 0; j < 9; j++)
            LU[i][j] = A[i * 9 + j];
    }

    for (int k = 0; k < 9; k++)
    {
        int max_row = k;
        float max_val = fabsf(LU[k][k]);
        for (int i = k + 1; i < 9; i++)
            if (fabsf(LU[i][k]) > max_val)
            {
                max_val = fabsf(LU[i][k]);
                max_row = i;
            }

        if (max_row != k)
        {
            int tmp = piv[k];
            piv[k] = piv[max_row];
            piv[max_row] = tmp;
            for (int j = 0; j < 9; j++)
            {
                float t = LU[k][j];
                LU[k][j] = LU[max_row][j];
                LU[max_row][j] = t;
            }
        }

        for (int i = k + 1; i < 9; i++)
        {
            LU[i][k] /= LU[k][k];
            for (int j = k + 1; j < 9; j++)
                LU[i][j] -= LU[i][k] * LU[k][j];
        }
    }

    for (int c = 0; c < nrhs; c++)
    {
        float y[9];
        for (int i = 0; i < 9; i++)
            y[i] = B[piv[i] * nrhs + c];

        for (int i = 1; i < 9; i++)
            for (int j = 0; j < i; j++)
                y[i] -= LU[i][j] * y[j];

        // Note: The starting index for backward substitution changes from 5 to 8
        for (int i = 8; i >= 0; i--)
        {
            for (int j = i + 1; j < 9; j++)
                y[i] -= LU[i][j] * y[j];
            y[i] /= LU[i][i];
        }

        for (int i = 0; i < 9; i++)
            X[i * nrhs + c] = y[i];
    }
}

void mekf_wb(const MEKF_State *input, MEKF_State *output, const float *omega, const float *Br, const float *Nr)
{
    const float *qkk = input->q;
    const float *bkk = input->b;
    const float *Pkk = input->P;
    const float *W = input->W;
    const float *V = input->V;
    const float dt = input->dt;

    float omega_corr[3] = {0};
    omega_corr[0] = omega[0] - bkk[0];
    omega_corr[1] = omega[1] - bkk[1];
    omega_corr[2] = omega[2] - bkk[2];

    // 1. Prediction
    float k0 = 0.5 * dt;
    float phi[3] = {0};
    phi[0] = k0 * omega_corr[0];
    phi[1] = k0 * omega_corr[1];
    phi[2] = k0 * omega_corr[2];

    float dqkk[4] = {0};
    expq(phi, dqkk);

    float Lqkk[16] = {0};
    Lq(qkk, Lqkk); // Lqkk = L(qkk)
    float qk1k[4] = {0};
    Mprod(Lqkk, dqkk, 4, 4, 1, 1, qk1k); // qk1k = L(qkk)*dqkk

    float Rdqkk[16] = {0};
    Rq(dqkk, Rdqkk); // Rdqkk = R(dqkk)

    float Gqkk[12] = {0};
    Gq(qkk, Gqkk); // Gqkk = G(qkk)

    float temp0[12] = {0};
    Mprod(Rdqkk, Gqkk, 4, 4, 3, 1, temp0); // temp0 = R(dqkk)*G(qkk)

    float Gqk1k[12] = {0};
    Gq(qk1k, Gqk1k); // Gqk1k = G(qk1k)

    float Gqk1kT[12] = {0};
    Mtranspose(Gqk1k, 4, 3, Gqk1kT); // Gqk1kT = G(qk1k)'

    float Ak11[9] = {0}, Ak12[9] = {0}, Ak[36] = {0};
    Mprod(Gqk1kT, temp0, 3, 4, 3, 1, Ak11);  // Ak11 = G(qk1k)'*temp0 = G(qk1k)'*R(dqkk)*G(qkk)
    Mprod(Gqk1kT, Gqkk, 3, 4, 3, -k0, Ak12); // Ak12 = -0.5*dt*G(qk1k)'*G(qkk), k0 = 0.5*dt
    AkTransform(Ak11, Ak12, Ak);             // Ak = [Ak11, Ak12; zeros(3,3), eye(3)]

    float AkT[36] = {0}, Pk1k[36] = {0};
    Mtranspose(Ak, 6, 6, AkT); // AkT = Ak'
    float temp2[36] = {0}, temp3[36] = {0};
    Mprod(Pkk, AkT, 6, 6, 6, 1, temp2);  // temp2 = Pkk*Ak'
    Mprod(Ak, temp2, 6, 6, 6, 1, temp3); // temp3 = Ak*temp2 = Ak*Pkk*Ak'
    Madd(temp3, V, 6, 6, Pk1k);          // Pk1k = temp3 + V = Ak*Pkk*Ak' + V

    // 2. Innovation
    float H[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1}, HT[12] = {0};
    Mtranspose(H, 4, 3, HT); // HT = H'

    float Rqk1k[16] = {0}, Lqk1k[16] = {0}, temp4[12] = {0}, temp5[12] = {0};
    Rq(qk1k, Rqk1k); // Rqk1k = R(qk1k)
    Lq(qk1k, Lqk1k); // Lqk1k = L(qk1k)

    Mprod(Rqk1k, H, 4, 4, 3, 1, temp4); // temp4 = R(qk1k)*H

    float Lqk1kT[16] = {0};
    Mtranspose(Lqk1k, 4, 4, Lqk1kT);

    Mprod(Lqk1kT, temp4, 4, 4, 3, 1, temp5); // temp5 = L(qk1k)'*temp4 = L(qk1k)'*R(qk1k)*H

    float Qk1[9] = {0}, Qfull[81] = {0};
    Mprod(HT, temp5, 3, 4, 3, 1, Qk1); // Qk1 = HT*temp5 = HT*L(qk1k)'*R(qk1k)*H
    Qk1Transform(Qk1, Qfull);          // Qfull = [Qk1, zeros(3,3), zeros(3,3), Qk1]

    float zk1[9] = {0}, temp6[9] = {0};
    Mprod(Qfull, Nr, 9, 9, 1, 1, temp6); // temp6 = Qfull*Nr
    Msub(Br, temp6, 9, 1, zk1);          // zk1 = Br - temp6 = Br - Qfull*Nr

    float Nr_sun[3] = {Nr[0], Nr[1], Nr[2]}, Nr_mag[3] = {Nr[3], Nr[4], Nr[5]}, Nr_ehs[3] = {Nr[6], Nr[7], Nr[8]};
    float T[16] = {1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1};

    // Sun Sensor Measurements (Unit Vector)
    float temp7[4] = {0}, temp8[16] = {0}, temp9[16] = {0}, temp10[16] = {0}, temp11[16] = {0}, temp12[12] = {0};

    float Lqt7[16] = {0}, Rqt7[16] = {0}, Ck1_11[9] = {0};

    Mprod(H, Nr_sun, 4, 3, 1, 1, temp7); // temp7 = H*Nr_sun

    Lq(temp7, Lqt7);
    Rq(temp7, Rqt7);

    Mprod(Lqk1kT, Lqt7, 4, 4, 4, 1, temp8);   // temp8 = L(qk1k)'*L(H*Nr_sun)
    Mprod(Rqt7, T, 4, 4, 4, 1, temp9);        // temp9 = R(H*Nr_sun)*T
    Mprod(Rqk1k, temp9, 4, 4, 4, 1, temp10);  // temp10 = R(qk1k)*R(H*Nr_sun)*T
    Madd(temp8, temp10, 4, 4, temp11);        // temp11 = L(qk1k)'*L(H*Nr_sun) + R(qk1k)*R(H*Nr_sun)*T
    Mprod(temp11, Gqk1k, 4, 4, 3, 1, temp12); // temp12 = (L(qk1k)'*L(H*Nr_sun) + R(qk1k)*R(H*Nr_sun)*T)*G(qk1k)
    Mprod(HT, temp12, 3, 4, 3, 1, Ck1_11);    // Ck1_11 = H'*(L(qk1k)'*L(H*Nr_sun) + R(qk1k)*R(H*Nr_sun)*T)*G(qk1k)

    // Magnetometer Measurements (Unit Vector)
    float temp13[4] = {0}, temp14[16] = {0}, temp15[16] = {0}, temp16[16] = {0}, temp17[16] = {0}, temp18[12] = {0};

    float Lqt13[16] = {0}, Rqt13[16] = {0}, Ck1_21[9] = {0};

    Mprod(H, Nr_mag, 4, 3, 1, 1, temp13); // temp13 = H*Nr_mag

    Lq(temp13, Lqt13);
    Rq(temp13, Rqt13);

    Mprod(Lqk1kT, Lqt13, 4, 4, 4, 1, temp14); // temp14 = L(qk1k)'*L(H*Nr_mag)
    Mprod(Rqt13, T, 4, 4, 4, 1, temp15);      // temp15 = R(H*Nr_mag)*T
    Mprod(Rqk1k, temp15, 4, 4, 4, 1, temp16); // temp16 = R(qk1k)*R(H*Nr_mag)*T
    Madd(temp14, temp16, 4, 4, temp17);       // temp17 = L(qk1k)'*L(H*Nr_mag) + R(qk1k)*R(H*Nr_mag)*T
    Mprod(temp17, Gqk1k, 4, 4, 3, 1, temp18); // temp18 = (L(qk1k)'*L(H*Nr_mag) + R(qk1k)*R(H*Nr_mag)*T)*G(qk1k)
    Mprod(HT, temp18, 3, 4, 3, 1, Ck1_21);    // Ck1_21 = H'*(L(qk1k)'*L(H*Nr_mag) + R(qk1k)*R(H*Nr_mag)*T)*G(qk1k)

    // Earth Horizon Sensor Measurements (Unit Vector)
    float temp29[4] = {0}, temp30[16] = {0}, temp31[16] = {0}, temp32[16] = {0}, temp33[16] = {0}, temp34[12] = {0};

    float Lqt29[16] = {0}, Rqt29[16] = {0}, Ck1_31[9] = {0};

    Mprod(H, Nr_ehs, 4, 3, 1, 1, temp29); // temp29 = H*Nr_ehs

    Lq(temp29, Lqt29);
    Rq(temp29, Rqt29);

    Mprod(Lqk1kT, Lqt29, 4, 4, 4, 1, temp30); // temp30 = L(qk1k)'*L(H*Nr_ehs)
    Mprod(Rqt29, T, 4, 4, 4, 1, temp31);      // temp31 = R(H*Nr_ehs)*T
    Mprod(Rqk1k, temp31, 4, 4, 4, 1, temp32); // temp32 = R(qk1k)*R(H*Nr_ehs)*T
    Madd(temp30, temp32, 4, 4, temp33);       // temp33 = L(qk1k)'*L(H*Nr_ehs) + R(qk1k)*R(H*Nr_ehs)*T
    Mprod(temp33, Gqk1k, 4, 4, 3, 1, temp34); // temp34 = (L(qk1k)'*L(H*Nr_ehs) + R(qk1k)*R(H*Nr_ehs)*T)*G(qk1k)
    Mprod(HT, temp34, 3, 4, 3, 1, Ck1_31);    // Ck1_31 = H'*(L(qk1k)'*L(H*Nr_ehs) + R(qk1k)*R(H*Nr_ehs)*T)*G(qk1k)

    float Ck1[54] = {0}, Ck1T[54] = {0};
    Ck1Transform(Ck1_11, Ck1_21, Ck1_31, Ck1); // Ck1 = [[Ck1_11; Ck1_21; Ck1_31], zeros(9,3)]
    Mtranspose(Ck1, 9, 6, Ck1T);               // Ck1T =  Ck1'

    float Sk1[81] = {0}, temp19[54] = {0}, temp20[81] = {0};
    Mprod(Pk1k, Ck1T, 6, 6, 9, 1, temp19);  // temp19 = Pk1k*Ck1'
    Mprod(Ck1, temp19, 9, 6, 9, 1, temp20); // temp20 = Ck1*Pk1k*Ck1'
    Madd(temp20, W, 9, 9, Sk1);             // Sk1 = Ck1*Pk1k*Ck1' + W

    // 3. Kalman Gain

    // float Kk1[36] = {0}, invSk1[36] = {0};
    // Inv6x6(Sk1, invSk1);                    // invSk1 = Sk1^-1
    // Mprod(temp19, invSk1, 6, 6, 6, 1, Kk1); // Kk1 = Pk1k*Ck1'*Sk1^-1`

    // Solve Sk1*Kk1^T = temp19^T, Sk1 is symmetric
    float Kk1[54] = {0}, temp19T[54] = {0}, Kk1T_solve[54] = {0};
    Mtranspose(temp19, 6, 9, temp19T);
    LUSolve9x9(Sk1, temp19T, 6, Kk1T_solve);
    Mtranspose(Kk1T_solve, 9, 6, Kk1);

    // 4. Update
    float dx[6] = {0};
    Mprod(Kk1, zk1, 6, 9, 1, 1, dx);

    float phi_corr[3] = {dx[0], dx[1], dx[2]}, db[3] = {dx[3], dx[4], dx[5]};

    float bk1k1[3] = {0};
    Madd(bkk, db, 3, 1, bk1k1); // bk1k1 = bkk + db

    float qk1k1temp[4] = {0}, temp21[4] = {0};
    expq(phi_corr, temp21);
    Mprod(Lqk1k, temp21, 4, 4, 1, 1, qk1k1temp); // qk1k1 = L(qk1k)*expq(phi_corr)

    float q_norm = sqrtf(qk1k1temp[0] * qk1k1temp[0] + qk1k1temp[1] * qk1k1temp[1] + qk1k1temp[2] * qk1k1temp[2] + qk1k1temp[3] * qk1k1temp[3]);
    float qk1k1[4] = {0};
    qk1k1[0] = qk1k1temp[0] / q_norm;
    qk1k1[1] = qk1k1temp[1] / q_norm;
    qk1k1[2] = qk1k1temp[2] / q_norm;
    qk1k1[3] = qk1k1temp[3] / q_norm;

    float Kk1T[54] = {0}, temp22[54] = {0}, temp23[36] = {0}, temp24[36] = {0}, temp25[36] = {0}, temp26[36] = {0}, temp27[36] = {0}, temp28[36] = {0};
    Mtranspose(Kk1, 6, 9, Kk1T);            // Kk1T = Kk1'
    Mprod(W, Kk1T, 9, 9, 6, 1, temp22);     // temp22 = W*Kk1'
    Mprod(Kk1, temp22, 6, 9, 6, 1, temp23); // temp23 = Kk1*W*Kk1'
    Mprod(Kk1, Ck1, 6, 9, 6, 1, temp24);    // temp24 = Kk1*Ck1
    float eye6x6[36] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
    Msub(eye6x6, temp24, 6, 6, temp25);        // temp25 = I - Kk1*Ck1
    Mtranspose(temp25, 6, 6, temp26);          // temp26 = (I - Kk1*Ck1)'
    Mprod(Pk1k, temp26, 6, 6, 6, 1, temp27);   // temp27 = Pk1k*(I - Kk1*Ck1)'
    Mprod(temp25, temp27, 6, 6, 6, 1, temp28); // temp28 = (I - Kk1*Ck1)*Pk1k*(I - Kk1*Ck1)'
    float Pk1k1[36] = {0};
    Madd(temp28, temp23, 6, 6, Pk1k1); // Pk1k1 = (I - Kk1*Ck1)*Pk1k*(I - Kk1*Ck1)' + Kk1*W*Kk1'

    memcpy(output->q, qk1k1, 4 * sizeof(float));
    memcpy(output->b, bk1k1, 3 * sizeof(float));
    memcpy(output->P, Pk1k1, 36 * sizeof(float));
}