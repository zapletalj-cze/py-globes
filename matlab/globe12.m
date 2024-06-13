clc
clear
axis equal

%global parameters
M=100000000;
R = 6380 *1000;
R= R/M;
proj = @gnom;
Du = 10;
Dv = 10;
du = 1;
dv = 1;

ua = 26.5651;      
ug = 52.6226 ;       
ud = 10.8123 ; 

A = [-ug, 0];
B = [-ug, 72];
C = [-ug, 144];
D = [-ug, 216]; 
E = [-ug, 288]; 

F = [-ud,0];
G = [ud, 36];
H = [-ud,72];
I = [ud, 108];
J = [-ud,144];
K = [ud, 180];
L = [-ud,216]; 
M = [ud, 252]; 
N = [-ud,288]; 
O = [ud, 324]; 

P = [ug, 36];
Q = [ug, 108];
RR = [ug, 180];
S = [ug, 252]; 
T = [ug, 324]; 

%% Face 1 (FGPTO, k.pol K1)
uk = ua;
vk = 0;

ub = [F(1), G(1), P(1), T(1), O(1), F(1)];
vb = [F(2), G(2), P(2), T(2), O(2), F(2)];
umin = -20;
umax = 62;
vmin = -41;
vmax = 41;

%figure(1)
subplot(3,4,1)
title('Face 1')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 2 (GHIQP, k.pol K2)
uk = ua;
vk = 72;

ub = [G(1), H(1), I(1), Q(1), P(1), G(1)];
vb = [G(2), H(2), I(2), Q(2), P(2), G(2)];
%u1 = -15; u2 = 57; v1 = 31; v2 = 112;
umin = -20;
umax = 62;
vmin = 31;
vmax = 112;

%figure(2)
subplot(3,4,2)
title('Face 2')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 3 (IJKRQ, k.pol K3)
uk = ua;
vk = 144;
umin = -20;
umax = 62;
vmin = 105;
vmax = 185;
ub = [I(1), J(1), K(1), RR(1), Q(1), I(1)];
vb = [I(2), J(2), K(2), RR(2), Q(2), I(2)];

%figure(3)
subplot(3,4,3)
title('Face 3')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 4 (LKRSM, k.pol K4)
uk = ua;
vk = 216;
umin = -20;
umax = 62;
vmin = 175;
vmax = 260;
ub = [L(1), K(1), RR(1), S(1), M(1), L(1)];
vb = [L(2), K(2), RR(2), S(2), M(2), L(2)];

%figure(4)
subplot(3,4,4)
title('Face 4')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 5 (ONMST, k.pol K5)
uk = ua;
vk = 288;
umin = -20;
umax = 62;
vmin = 240;
vmax = 330;
ub = [O(1), N(1), M(1), S(1), T(1), O(1)];
vb = [O(2), N(2), M(2), S(2), T(2), O(2)];

%figure(5)
subplot(3,4,5)
title('Face 5')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 6 (ABHGF, k.pol K6)
uk = -ua;
vk = 36;
umin = -62;
umax = 20;
vmin = -5;
vmax = 77;
ub = [A(1), B(1), H(1), G(1), F(1), A(1)];
vb = [A(2), B(2), H(2), G(2), F(2), A(2)];

%figure(6)
subplot(3,4,6)
title('Face 6')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 7 (BCJIH, k.pol K7)
uk = -ua;
vk = 108;
umin = -62;
umax = 20;
vmin = 65;
vmax = 155;
ub = [B(1), C(1), J(1), I(1), H(1), B(1)];
vb = [B(2), C(2), J(2), I(2), H(2), B(2)];

%figure(7)
subplot(3,4,7)
title('Face 7')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 8 (CJKLD, k.pol K8)
uk = -ua;
vk = 180;
umin = -62;
umax = 20;
vmin = 140;
vmax = 220;
ub = [C(1), J(1), K(1), L(1), D(1), C(1)];
vb = [C(2), J(2), K(2), L(2), D(2), C(2)];

%figure(8)
subplot(3,4,8)
title('Face 8')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 9 (EDLMN, k.pol K9)
uk = -ua;
vk = 252;
umin = -62;
umax = 20;
vmin = 210;
vmax = 295;
ub = [E(1), D(1), L(1), M(1), N(1), E(1)];
vb = [E(2), D(2), L(2), M(2), N(2), E(2)];

%figure(9)
subplot(3,4,9)
title('Face 9')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 10 (AENOF, k.pol K10)
uk = -ua;
vk = 324;
umin = -62;
umax = 20;
vmin = -75;
vmax = 5;
ub = [A(1), E(1), N(1), O(1), F(1), A(1)];
vb = [A(2), E(2), N(2), O(2), F(2), A(2)];

%figure(10)
subplot(3,4,10)
title('Face 10')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 11 (vrchni ploska PQRST, severni pol K11)
uk = 90;
vk = 0;
umin = 38;
umax = 90;
vmin = -180;
vmax = 180;
ub = [P(1), Q(1), RR(1), S(1), T(1), P(1)];
vb = [P(2), Q(2), RR(2), S(2), T(2), P(2)];

%figure(11)
subplot(3,4,11)
title('Face 11')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%% Face 12 (spodni ploska ABCDE, jizni pol K12)
uk = -90;
vk = 0;
umin = -90;
umax = -38;
vmin = -180;
vmax = 180;
ub = [A(1), B(1), C(1), D(1), E(1), A(1)];        
vb = [A(2), B(2), C(2), D(2), E(2), A(2)]; 

%figure(12)
subplot(3,4,12)
title('Face 12')
myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

