function [xk1k1, Pk1k1] = mekf_wb(dt, xkk, Pkk, wk, V, W, Byk1, Nyk1)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Attitude Estimator - Multiplicative Extended Kalman Filter %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dt - Filter Estimation Period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % xkk - [qkk; bkk] Previous Estimated State %%%%%%%%%%%%%%%%%%
    % Pkk - Previous Estimated Covariance Matrix %%%%%%%%%%%%%%%%%
    % wk - Gyroscope Measurement (rad/s) %%%%%%%%%%%%%%%%%%%%%%%%%
    % V - Process Noise Covariance Matrix (Gyro) %%%%%%%%%%%%%%%%%
    % W - Measurement Noise Covariance Matrix (Sun + Mag) %%%%%%%%
    % Byk1 - [Br_sun; Br_mag] Measurements in the Body Frame %%%%%
    % Nyk1 - [Nr_sun; Nr_mag] Measurements in the Inertial Frame %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Outputs: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % xk1k1 - [qk1k1; bk1k1] Current State Estimation %%%%%%%%%%%%
    % Pk1k1 - Current Covariance Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xk1k1 = zeros(7,1);
    Pk1k1 = zeros(6,6);
    
    H = [zeros(1,3); eye(3)];
    T = diag([1, -1, -1, -1]);

    qkk = xkk(1:4);
    bkk = xkk(5:7);

    % 1. Prediction
    phi = 0.5*dt*(wk - bkk);
    dqkk = expq(phi);

    qk1k = Lq(qkk)*dqkk;

    Ak11 = Gq(qk1k)'*Rq(dqkk)*Gq(qkk);
    Ak12 = -0.5*dt*Gq(qk1k)'*Gq(qkk);
    Ak = [Ak11, Ak12; zeros(3,3), eye(3)];

    Pk1k = Ak*Pkk*Ak' + V;

    % 2. Innovation
    Qk1k = H'*Lq(qk1k)'*Rq(qk1k)*H;
    zk1 = Byk1 - [Qk1k, zeros(3,3); zeros(3,3), Qk1k]*Nyk1;
    
    NrSun = Nyk1(1:3);
    NrMag = Nyk1(4:6);
    Ck1_11 = H'*(Lq(qk1k)'*Lq(H*NrSun) + Rq(qk1k)*Rq(H*NrSun)*T)*Gq(qk1k);
    Ck1_21 = H'*(Lq(qk1k)'*Lq(H*NrMag) + Rq(qk1k)*Rq(H*NrMag)*T)*Gq(qk1k);
    Ck1 = [[Ck1_11; Ck1_21], zeros(6,3)];

    Sk1 = Ck1*Pk1k*Ck1' + W;

    % 3. Kalman Gain
    Kk1 = Pk1k*Ck1'/Sk1;
       
    % 4. Update
    dx = Kk1*zk1;               
    phi_corr = dx(1:3);         
    dbeta = dx(4:6);
    bk1k1 = bkk + dbeta;
    qk1k1 = Lq(qk1k)*[sqrt(1 - phi_corr'*phi_corr); phi_corr];
    Pk1k1 = (eye(6) - Kk1*Ck1)*Pk1k*(eye(6) - Kk1*Ck1)' + Kk1*W*Kk1';

    xk1k1 = [qk1k1; bk1k1];

end

function res = expq(v)
    norm_v = norm(v);
    qs = cos(norm_v);
    qv = v*sinc(norm_v/pi);
    res = [qs; qv];
end

function rq = Rq(q)
    rq = [q(1), -q(2:4)'; q(2:4), q(1)*eye(3) - hat(q(2:4))];
end

function lq = Lq(q)
    lq = [q(1), -q(2:4)'; q(2:4), q(1)*eye(3) + hat(q(2:4))];
end

function gq = Gq(q)
    H = [zeros(1,3); eye(3)];
    gq = Lq(q)*H;
end

function S = hat(v)
    S = [  0,   -v(3),  v(2);
          v(3),   0,   -v(1);
         -v(2),  v(1),   0  ];
end