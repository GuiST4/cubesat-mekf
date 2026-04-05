# NanoSat MEKF Attitude Estimator 🛰️

A C implementation of a Multiplicative Extended Kalman Filter (MEKF) for NanoSat attitude determination.

Designed for modern ARM Cortex-M microcontrollers with Hardware Floating-Point Units (like the **STM32 F7 and H7 series**), this project adapts double-precision simulation algorithms for 32-bit single-precision embedded environments.

## 🚀 Key Features

*   **Sensor Fusion:** Couples gyroscope rate integration with vector measurements from a **Sun Sensor and Magnetometer**.
*   **Joseph Form Covariance:** Maintains a positive semi-definite $P$ matrix to help mitigate 32-bit float truncation errors over time.
*   **Memory Efficiency:** Uses a ~6 KB stack footprint with flattened 1D arrays, avoiding dynamic memory allocation.

---

## 🛰️ Integration for Flight Software 

The filter is encapsulated into a single, state-transition function. It takes the current filter state and physical sensor vectors, and outputs the propagated state for the next time step.

```c
#include "mekf_wb.h"
#include <string.h>

MEKF_State current_state;
MEKF_State next_state;

void init() {
    // Initialize quaternion (e.g. identity), bias, and timestep
    current_state.q[0] = 1.0f; current_state.q[1] = 0.0f;
    current_state.q[2] = 0.0f; current_state.q[3] = 0.0f;
    current_state.b[0] = 0.0f; current_state.b[1] = 0.0f; current_state.b[2] = 0.0f;
    current_state.dt = 0.01f; // 100 Hz

    // Note: Covariance (P), Process Noise (V), and Measurement Noise (W) 
}

void loop() {
    // 1. Gather sensor data and reference vectors
    float gyro[3];  // Angular rates from Gyroscope (rad/s)
    float Br[6];    // Body frame measurements [Sun_x, Sun_y, Sun_z, Mag_x, Mag_y, Mag_z] (Unit vectors)
    float Nr[6];    // Inertial reference vectors [Sun_x, Sun_y, Sun_z, Mag_x, Mag_y, Mag_z] (Unit vectors)
    
    ReadSensors(gyro, Br, Nr); // Example hardware fetch

    // 2. Run the MEKF
    mekf_wb(&current_state, &next_state, gyro, Br, Nr);

    // 3. Propagate the state forward for the next cycle
    memcpy(current_state.q, next_state.q, 4 * sizeof(float));
    memcpy(current_state.b, next_state.b, 3 * sizeof(float));
    memcpy(current_state.P, next_state.P, 36 * sizeof(float));
}
```

## 📊 Validation & Performance

This filter was tested against a MATLAB/Simulink NanoSat Simulator. The embedded C implementation closely matches the numerical performance of the double-precision Simulink reference model.

### 1. Attitude Tracking
Tracks the simulated true state across all four quaternion components.


### 2. Gyro Bias Estimation
Estimates and removes gyroscope biases using the vector measurements for drift correction.




### 3. Absolute Attitude Error
Shows the physical pointing error of the C implementation compared to the simulation environment.



### 4. Implementation Accuracy (C vs. Simulink)
Compares the C code output directly against the Simulink MEKF. The Principal Rotation Angle difference stays around **~0.0001 degrees**, showing strong numerical agreement.



---

## 💻 Hardware Integration Notes

*   **FPU Requirement:** Ensure your compiler flags have hardware floating-point math enabled. Emulated floating-point math will noticeably impact execution speed.
*   **Sample Rates:** The filter structure assumes a higher frequency prediction step (Gyro) and a lower frequency innovation step (Sun/Mag). Ensure the `dt` parameter is updated in your hardware timer interrupts to handle multi-rate sensors.

## 📁 File Structure

*   `mekf_wb.c` / `mekf_wb.h` - Core filter logic and math operations.
*   `main.c` - Example test wrapper for processing `.csv` telemetry.
