%% pmsm_foc_sim.m
% Field-Oriented Control (FOC) of a PMSM, simulated as a discretely
% sampled digital drive (mirrors a real DSP/MCU implementation):
%
%   phase currents --Clarke--> alpha-beta --Park--> d-q
%                                                   |
%                     speed PI --> iq_ref     id_ref = 0
%                                                   |
%                         d/q current PI (+ cross-coupling decoupling)
%                                                   |
%                              inverse Park --> alpha-beta
%                                                   |
%                                   Space Vector PWM (SVPWM)
%                                                   |
%                       average inverter output voltage --> PMSM plant
%
% The PMSM electrical/mechanical plant is integrated exactly in the
% rotor dq frame (RK4). To exercise the *actual* Clarke/Park code on
% real three-phase quantities, the dq plant currents are converted to
% synthesized phase currents (inverse Park + inverse Clarke) each
% sample, as if read from current sensors, and then transformed back
% through Clarke+Park for the controller feedback. This gives bit-true
% transform code while keeping the electrical plant exact.

clear; clc; close all;

%% ---- Motor parameters ----
mp.R      = 2.875;     % Stator resistance (Ohm)
mp.Ld     = 8.5e-3;    % d-axis inductance (H)
mp.Lq     = 8.5e-3;    % q-axis inductance (H)   (Ld=Lq: surface PMSM)
mp.lambda = 0.175;     % PM flux linkage (Wb)
mp.P      = 4;         % Number of poles
mp.J      = 0.0008;    % Rotor inertia (kg.m^2)
mp.B      = 0.0001;    % Viscous friction (N.m.s)
Kt        = 1.5*(mp.P/2)*mp.lambda;   % Torque constant (N.m/A), Ld=Lq case

%% ---- Inverter / control settings ----
Vdc      = 300;            % DC bus voltage (V)
Ts       = 50e-6;          % Current-loop / PWM sample time (s) -> 20 kHz
Nspeed   = 20;              % Speed loop runs every Nspeed current cycles (1 kHz)
Tend     = 0.5;             % Simulation length (s)
N        = round(Tend/Ts);
Imax     = 10;               % Current (torque) limit, A
Vmax_lin = Vdc/sqrt(3);      % Max voltage magnitude for linear SVPWM

%% ---- Controller gains (pole-cancellation / IMC tuning) ----
wc_i = 2*pi*800;   % current-loop bandwidth (rad/s)
Kp_i = wc_i*mp.Ld; Ki_i = wc_i*mp.R;          % same L for d/q since Ld=Lq

wc_s = 2*pi*60;    % speed-loop bandwidth (rad/s), well below wc_i
Kp_s = wc_s*mp.J;  Ki_s = wc_s*mp.B;          % speed PI acting on Te

%% ---- References / disturbance profile ----
wm_ref_fun = @(t) 150*ones(size(t));          % speed reference, rad/s (~1432 RPM)
TL_fun     = @(t) 2.0*(t>=0.3);               % load torque step 2 N.m at t=0.3s

%% ---- State init ----
id = 0; iq = 0; wm = 0; theta_e = 0;
id_integ = 0; iq_integ = 0; speed_integ = 0;
iq_ref = 0; Te_max = Kt*Imax;

% Logging
t_log  = zeros(1,N); id_log = t_log; iq_log = t_log; wm_log = t_log;
idref_log = t_log; iqref_log = t_log; Te_log = t_log;
ia_log = t_log; ib_log = t_log; ic_log = t_log;
da_log = t_log; db_log = t_log; dc_log = t_log;

%% ---- Main discrete control loop ----
for k = 1:N
    t = (k-1)*Ts;
    we = (mp.P/2)*wm;                          % electrical speed

    % ---- "Measure" phase currents from the true dq plant state ----
    [ialpha_m, ibeta_m] = inv_park_transform(id, iq, theta_e);
    [ia_m, ib_m, ic_m]  = inv_clarke_transform(ialpha_m, ibeta_m);

    % ---- Clarke + Park on the measured phase currents (feedback path) ----
    [ialpha_c, ibeta_c] = clarke_transform(ia_m, ib_m, ic_m);
    [id_meas, iq_meas]  = park_transform(ialpha_c, ibeta_c, theta_e);

    % ---- Outer speed PI loop (slower update rate) ----
    if mod(k-1, Nspeed) == 0
        err_w = wm_ref_fun(t) - wm;
        [Te_ref, speed_integ] = pi_control(err_w, speed_integ, Kp_s, Ki_s, ...
                                            Ts*Nspeed, -Te_max, Te_max);
        iq_ref = saturate(Te_ref/Kt, -Imax, Imax);
    end
    id_ref = 0;   % zero d-axis current control (surface PMSM, no field weakening)

    % ---- Inner d/q current PI loops ----
    err_id = id_ref - id_meas;
    err_iq = iq_ref - iq_meas;
    [vd_pi, id_integ] = pi_control(err_id, id_integ, Kp_i, Ki_i, Ts, -Vdc, Vdc);
    [vq_pi, iq_integ] = pi_control(err_iq, iq_integ, Kp_i, Ki_i, Ts, -Vdc, Vdc);

    % ---- Cross-coupling decoupling feedforward ----
    vd_cmd = vd_pi - we*mp.Lq*iq_meas;
    vq_cmd = vq_pi + we*mp.Ld*id_meas + we*mp.lambda;

    % ---- Limit reference vector to the linear SVPWM range ----
    Vmag = hypot(vd_cmd, vq_cmd);
    if Vmag > Vmax_lin
        scale  = Vmax_lin/Vmag;
        vd_cmd = vd_cmd*scale;
        vq_cmd = vq_cmd*scale;
    end

    % ---- Inverse Park: dq -> alpha-beta voltage ----
    [valpha, vbeta] = inv_park_transform(vd_cmd, vq_cmd, theta_e);

    % ---- Space Vector PWM: alpha-beta -> duty cycles ----
    [da, db, dc] = svpwm(valpha, vbeta, Vdc);

    % ---- Reconstruct the actual average voltage the inverter applies ----
    van = Vdc*(da-0.5); vbn = Vdc*(db-0.5); vcn = Vdc*(dc-0.5);
    vcm = (van+vbn+vcn)/3;                          % common-mode component
    va_n = van-vcm; vb_n = vbn-vcm; vc_n = vcn-vcm;  % isolated-neutral star voltages
    [valpha_act, vbeta_act] = clarke_transform(va_n, vb_n, vc_n);
    [vd_act, vq_act]        = park_transform(valpha_act, vbeta_act, theta_e);

    % ---- Integrate the PMSM plant one Ts step (RK4) ----
    TL = TL_fun(t);
    [id, iq, wm] = pmsm_step(id, iq, wm, vd_act, vq_act, TL, mp, Ts);
    theta_e = wrap_to_pi(theta_e + we*Ts);

    % ---- Log ----
    t_log(k)=t; id_log(k)=id; iq_log(k)=iq; wm_log(k)=wm;
    idref_log(k)=id_ref; iqref_log(k)=iq_ref;
    Te_log(k) = 1.5*(mp.P/2)*(mp.lambda*iq + (mp.Ld-mp.Lq)*id*iq);
    ia_log(k)=ia_m; ib_log(k)=ib_m; ic_log(k)=ic_m;
    da_log(k)=da; db_log(k)=db; dc_log(k)=dc;
end

rpm_log = wm_log*60/(2*pi);

%% ---- Plots: control performance ----
figure('Name','FOC PMSM drive performance');

subplot(2,2,1);
plot(t_log, wm_ref_fun(t_log), 'k--', t_log, wm_log, 'b', 'LineWidth',1.3);
legend('\omega_{ref}','\omega_m','Location','best');
ylabel('Speed (rad/s)'); title('Speed response'); grid on;

subplot(2,2,2);
plot(t_log, Te_log, 'k', 'LineWidth',1.3);
ylabel('Torque (N.m)'); title('Electromagnetic torque'); grid on;

subplot(2,2,3);
plot(t_log, idref_log, 'k--', t_log, id_log, 'b', ...
     t_log, iqref_log, 'r--', t_log, iq_log, 'm', 'LineWidth',1.1);
legend('i_d^*','i_d','i_q^*','i_q','Location','best');
xlabel('Time (s)'); ylabel('Current (A)'); title('d/q currents'); grid on;

subplot(2,2,4);
plot(t_log, ia_log,'b', t_log, ib_log,'r', t_log, ic_log,'g');
xlim([0.30 0.32]);   % zoom in so the sinusoids are visible
xlabel('Time (s)'); ylabel('Current (A)'); title('Phase currents (zoom)'); grid on;

%% ---- Plot: SVPWM duty cycles ----
figure('Name','SVPWM duty cycles');
plot(t_log, da_log,'b', t_log, db_log,'r', t_log, dc_log,'g');
xlim([0.30 0.32]);
xlabel('Time (s)'); ylabel('Duty ratio');
legend('d_a','d_b','d_c'); title('SVPWM duty cycles (zoom)'); grid on;

%% ---- Illustration: actual switched PWM waveform for phase A ----
% Zero-order-hold the duty ratio found in each Ts and compare it against
% a symmetric (center-aligned) triangular carrier to produce the real
% gate/switch waveform, exactly as a PWM timer peripheral would.
k0 = round(0.30/Ts); n_periods = 15;
sub = 40; % carrier resolution per period
tt = []; carrier = []; Sa = [];
for kk = k0:(k0+n_periods-1)
    tau = linspace(0, Ts, sub);
    c   = [2*tau(tau<=Ts/2)/Ts, 2*(1-tau(tau>Ts/2)/Ts)];   % 0->1->0 triangle
    tt      = [tt, (kk-1)*Ts + tau]; %#ok<AGROW>
    carrier = [carrier, c]; %#ok<AGROW>
    Sa      = [Sa, double(c <= da_log(kk))]; %#ok<AGROW>
end
figure('Name','PWM generation detail (phase A)');
subplot(2,1,1);
plot(tt, carrier, 'Color',[0.6 0.6 0.6]); hold on;
stairs(tt, interp1((k0:(k0+n_periods-1))'*Ts-Ts, da_log(k0:(k0+n_periods-1)), tt, 'previous','extrap'), 'b','LineWidth',1.3);
ylabel('Duty / carrier'); title('Carrier vs. commanded duty ratio d_a'); grid on;
subplot(2,1,2);
stairs(tt, Sa, 'k','LineWidth',1.3); ylim([-0.2 1.2]);
xlabel('Time (s)'); ylabel('Switch state S_a'); title('Resulting PWM switching waveform'); grid on;

%% ======================= Local functions =======================

function [ialpha, ibeta] = clarke_transform(ia, ib, ic)
% Amplitude-invariant Clarke transform (abc -> alpha-beta)
    ialpha = (2/3)*(ia - 0.5*ib - 0.5*ic);
    ibeta  = (2/3)*(sqrt(3)/2*ib - sqrt(3)/2*ic);
end

function [ia, ib, ic] = inv_clarke_transform(ialpha, ibeta)
% Inverse Clarke transform (alpha-beta -> abc)
    ia = ialpha;
    ib = -0.5*ialpha + (sqrt(3)/2)*ibeta;
    ic = -0.5*ialpha - (sqrt(3)/2)*ibeta;
end

function [id, iq] = park_transform(ialpha, ibeta, theta)
% Park transform (alpha-beta -> rotor dq), theta = electrical angle
    id =  ialpha*cos(theta) + ibeta*sin(theta);
    iq = -ialpha*sin(theta) + ibeta*cos(theta);
end

function [valpha, vbeta] = inv_park_transform(vd, vq, theta)
% Inverse Park transform (dq -> alpha-beta)
    valpha = vd*cos(theta) - vq*sin(theta);
    vbeta  = vd*sin(theta) + vq*cos(theta);
end

function [da, db, dc] = svpwm(valpha, vbeta, Vdc)
% Standard 7-segment symmetric SVPWM. Returns duty ratios in [0,1].
    Vref  = hypot(valpha, vbeta);
    theta = atan2(vbeta, valpha);
    if theta < 0, theta = theta + 2*pi; end

    sector = floor(theta/(pi/3)) + 1;
    th_s   = theta - (sector-1)*pi/3;   % angle within the sector, [0, pi/3)

    T1 = (sqrt(3)/Vdc)*Vref*sin(pi/3 - th_s);   % normalized to Ts_pwm = 1
    T2 = (sqrt(3)/Vdc)*Vref*sin(th_s);
    T0 = 1 - T1 - T2;
    if T0 < 0                                    % overmodulation guard
        scale = 1/(T1+T2);
        T1 = T1*scale; T2 = T2*scale; T0 = 0;
    end

    switch sector
        case 1, Ta=T1+T2+T0/2; Tb=T2+T0/2;    Tc=T0/2;
        case 2, Ta=T1+T0/2;    Tb=T1+T2+T0/2; Tc=T0/2;
        case 3, Ta=T0/2;       Tb=T1+T2+T0/2; Tc=T2+T0/2;
        case 4, Ta=T0/2;       Tb=T1+T0/2;    Tc=T1+T2+T0/2;
        case 5, Ta=T2+T0/2;    Tb=T0/2;       Tc=T1+T2+T0/2;
        otherwise, Ta=T1+T2+T0/2; Tb=T0/2;    Tc=T1+T0/2;
    end
    da = Ta; db = Tb; dc = Tc;
end

function [out, integ_new] = pi_control(err, integ, Kp, Ki, Ts, out_min, out_max)
% Discrete PI controller with simple integrator-clamping anti-windup
    integ_try = integ + err*Ts;
    out_unsat = Kp*err + Ki*integ_try;
    out = saturate(out_unsat, out_min, out_max);
    if out == out_unsat
        integ_new = integ_try;      % not saturated: keep integrating
    else
        integ_new = integ;          % saturated: freeze integrator
    end
end

function y = saturate(x, lo, hi)
    y = min(max(x, lo), hi);
end

function [id_n, iq_n, wm_n] = pmsm_step(id, iq, wm, vd, vq, TL, p, Ts)
% RK4 integration of the PMSM electrical+mechanical dq-frame dynamics
    x0 = [id; iq; wm];
    f  = @(x) pmsm_deriv(x, vd, vq, TL, p);
    k1 = f(x0);
    k2 = f(x0 + Ts/2*k1);
    k3 = f(x0 + Ts/2*k2);
    k4 = f(x0 + Ts*k3);
    x1 = x0 + Ts/6*(k1+2*k2+2*k3+k4);
    id_n=x1(1); iq_n=x1(2); wm_n=x1(3);
end

function dx = pmsm_deriv(x, vd, vq, TL, p)
    id=x(1); iq=x(2); wm=x(3);
    we = (p.P/2)*wm;
    did = (vd - p.R*id + we*p.Lq*iq)/p.Ld;
    diq = (vq - p.R*iq - we*p.Ld*id - we*p.lambda)/p.Lq;
    Te  = 1.5*(p.P/2)*(p.lambda*iq + (p.Ld-p.Lq)*id*iq);
    dwm = (Te - TL - p.B*wm)/p.J;
    dx  = [did; diq; dwm];
end

function th = wrap_to_pi(th)
    th = mod(th+pi, 2*pi) - pi;
end
