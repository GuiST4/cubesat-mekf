#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mekf_wb.h"

#define LINE_BUF 1024

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

    // V MATRIX (Process Noise)
    float sig_w = (0.1f / 100.0f) * 3.14159265f / 180.0f;
    float sig_b = (0.01f / 40.0f) * 3.14159265f / 180.0f;

    memset(state.V, 0, 36 * sizeof(float));

    // Gyro Noise
    state.V[0] = sig_w * sig_w;
    state.V[7] = sig_w * sig_w;
    state.V[14] = sig_w * sig_w;

    // Bias Noise
    state.V[21] = sig_b * sig_b;
    state.V[28] = sig_b * sig_b;
    state.V[35] = sig_b * sig_b;

    // W MATRIX (Measurement Noise)
    float sig_ehs = (0.25f * 10.0f) * 3.14159265f / 180.0f;

    memset(state.W, 0, 81 * sizeof(float));

    // Sensor 1: Sun Sensor
    state.W[0] = 0.1f * 0.1f;
    state.W[10] = 0.1f * 0.1f;
    state.W[20] = 0.1f * 0.1f;

    // Sensor 2: Magnetometer
    state.W[30] = 0.05f * 0.05f;
    state.W[40] = 0.05f * 0.05f;
    state.W[50] = 0.05f * 0.05f;

    // Sensor 3: Earth Horizon Sensor
    state.W[60] = sig_ehs * sig_ehs;
    state.W[70] = sig_ehs * sig_ehs;
    state.W[80] = sig_ehs * sig_ehs;

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
        float t, Nr[9], Br[9], gyro[3];

        sscanf(line, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
               &t,
               &Nr[0], &Nr[1], &Nr[2], &Nr[3], &Nr[4], &Nr[5], &Nr[6], &Nr[7], &Nr[8],
               &Br[0], &Br[1], &Br[2], &Br[3], &Br[4], &Br[5], &Br[6], &Br[7], &Br[8],
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