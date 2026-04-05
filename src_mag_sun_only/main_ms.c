#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mekf_wb_ms.h"

#define LINE_BUF 512

int main(void)
{
    MEKF_State state;
    state.q[0] = 1.0f;
    state.q[1] = 0.0f;
    state.q[2] = 0.0f;
    state.q[3] = 0.0f;
    state.b[0] = 0.0f;
    state.b[1] = 0.0f;
    state.b[2] = 0.0f;
    state.dt = 0.01f;

    float sig_att = 10.0f * 3.14159265f / 180.0f;
    float sig_bias = 0.5f * 3.14159265f / 180.0f;
    memset(state.P, 0, 36 * sizeof(float));
    state.P[0] = sig_att * sig_att;
    state.P[7] = sig_att * sig_att;
    state.P[14] = sig_att * sig_att;
    state.P[21] = sig_bias * sig_bias;
    state.P[28] = sig_bias * sig_bias;
    state.P[35] = sig_bias * sig_bias;

    float sig_w = 0.1f * 3.14159265f / 180.0f;
    float sig_b = 0.01f * 3.14159265f / 180.0f;
    memset(state.V, 0, 36 * sizeof(float));
    state.V[0] = sig_w * sig_w;
    state.V[7] = sig_w * sig_w;
    state.V[14] = sig_w * sig_w;
    state.V[21] = sig_b * sig_b;
    state.V[28] = sig_b * sig_b;
    state.V[35] = sig_b * sig_b;

    memset(state.W, 0, 36 * sizeof(float));
    state.W[0] = 0.01f * 0.01f;
    state.W[7] = 0.01f * 0.01f;
    state.W[14] = 0.01f * 0.01f;
    state.W[21] = 0.005f * 0.005f;
    state.W[28] = 0.005f * 0.005f;
    state.W[35] = 0.005f * 0.005f;

    FILE *fin = fopen("mekf_input.csv", "r");
    if (!fin)
    {
        printf("Cannot open mekf_input.csv\n");
        return 1;
    }

    FILE *fout = fopen("mekf_output.csv", "w");
    if (!fout)
    {
        printf("Cannot open mekf_output.csv\n");
        fclose(fin);
        return 1;
    }

    fprintf(fout, "t,q0,q1,q2,q3,b1,b2,b3\n");

    char line[LINE_BUF];
    fgets(line, LINE_BUF, fin);

    MEKF_State next;
    memcpy(&next, &state, sizeof(MEKF_State));

    int count = 0;
    while (fgets(line, LINE_BUF, fin))
    {
        float t, Nr[6], Br[6], gyro[3];

        sscanf(line, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
               &t,
               &Nr[0], &Nr[1], &Nr[2], &Nr[3], &Nr[4], &Nr[5],
               &Br[0], &Br[1], &Br[2], &Br[3], &Br[4], &Br[5],
               &gyro[0], &gyro[1], &gyro[2]);

        mekf_wb(&state, &next, gyro, Br, Nr);

        fprintf(fout, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                t,
                next.q[0], next.q[1], next.q[2], next.q[3],
                next.b[0], next.b[1], next.b[2]);

        memcpy(state.q, next.q, 4 * sizeof(float));
        memcpy(state.b, next.b, 3 * sizeof(float));
        memcpy(state.P, next.P, 36 * sizeof(float));

        count++;
    }

    printf("Processed %d samples\n", count);

    fclose(fin);
    fclose(fout);

    return 0;
}