close all;
clear;
clc;

global g m Iyy ut


%% calculating K matrix for linear terminal region (dual mode control) - manually calculated jacobian

g = 9.91;
m = 200;
Iyy = 5000;
ut = 200;

A = [
0.1, 0, 0, 0, 0, -g;
0.1, 0, 0, 0, ut, 0;
1, 0, 0, 0, 0, 0;
0, 1, 0, 0, 0, -ut;
0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 1, 0;
];

B = [
1/m, 0, 0;
0, 1/m, 0;
0, 0, 0;
0, 0, 0;
0, 0, Iyy;
0, 0, 0;
];

C = eye(6);

D = zeros(6, 3);

missile_lin = ss(A, B, C, D);

figure;
step(missile_lin);

poles_missile_lin = pole(missile_lin);

p = [-10; -10; 0; 0; -10; 0];
K = place(A,B,p);
Acl = A-B*K;
missile_lin_placed = ss(Acl,B,C,D);
poles_missile_lin_placed = pole(missile_lin_placed);

figure;
step(missile_lin_placed);

figure;
bode(missile_lin_placed);