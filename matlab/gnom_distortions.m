clc
clear

format long g
hold on 
axis equal

%symbolic variables
syms R u v

%gnom projection
y = -R * tan(pi/2-u) * cos(v);
x = R * tan(pi/2-u) * sin(v);

%Partial derivatives
fu = diff(x,u);
fu = simplify(fu, 'Steps',20);
fv = diff(x,v);
fv = simplify(fv, 'Steps',20); %pocať pro par. der. pro x
gu = diff(y,u);
gu = simplify(gu, 'Steps',20);
gv = diff(y,v);
gv = simplify(gv, 'Steps',20);

%local linear scales
mp2 = simplify((fu^2+gu^2)/R^2);
mr2 = simplify((fv^2+gv^2)/(R*cos(u))^2);
mp = simplify(sqrt(mp2));
mr = simplify(sqrt(mr2));

% area scale
p = 2*(fu*fv + gu*gv)/(R^2*cos(u));
P1 = simplify((gu*fv-gv*fu)/(R^2*cos(u)));
P2 = simplify(mp*mr);

%max angular distorsion
d_omega = 2*asin(abs(mp-mr)/(mp+mr)); % pro tissotova indikatrix mp = b, mr= a

%convergence
sigma_p = atan2(gu,fu);
c = abs(pi/2 - sigma_p);

% numerical parameters - pro severni plosku
Rn = 1 ; 
un = 52.6226 * pi/180;
vn = 324 * pi/180;

% coordinates
Xn = double(subs(x, {R, u, v}, {Rn, un, vn}));
Yn = double(subs(y, {R, u, v}, {Rn, un, vn}));

FUn = double(subs(fu, {R, u, v}, {Rn, un, vn}));
FVn = double(subs(fv, {R, u, v}, {Rn, un, vn}));

GUn = double(subs(gu, {R, u, v}, {Rn, un, vn}));
GVn = double(subs(gv, {R, u, v}, {Rn, un, vn}));

% local linear scales
MPn = double(subs(mp, {R, u, v}, {Rn, un, vn}));
MRn = double(subs(mr, {R, u, v}, {Rn, un, vn}));

% areal scale
Pn = double(subs(p, {R, u, v}, {Rn, un, vn}));
P1n = double(subs(P1, {R, u, v}, {Rn, un, vn}));
P2n = double(subs(P2, {R, u, v}, {Rn, un, vn}));

% maximum angular distortion
d_omegan = double(subs(d_omega, {R, u, v}, {Rn, un, vn})) * 180/pi;

% convergence
sigma_pn = double(subs(sigma_p, {R, u, v}, {Rn, un, vn}));
if sigma_pn < 0
    sigma_pn = sigma_pn + 2*pi;
end
c_n = double(subs(c, {R, u, v}, {Rn, un, vn})) * 180/pi;



% Tissot's indicatrix
t = 0: 5*pi/180 : 2*pi;
sc = 0.1; 
Xe = MPn * cos(t) * sc;
Ye = MRn * sin(t) * sc;


matrix_r = [cos(sigma_pn) -sin(sigma_pn); sin(sigma_pn) cos(sigma_pn)];

% rotate ellipse
XYe = [Xe; Ye];
XYe_rotated = matrix_r * XYe;

Xer = XYe_rotated(1,:); 
Yer = XYe_rotated(2,:);

% shift ellipse
Xer = Xer + Xn;
Yer = Yer + Yn;

% vytvoření sítě
plot (Xer,Yer);
D_u = 10;
D_v = 10;
d_u = 1;
d_v = 1;
proj = @gnom;
%R = 6380;
u0 = 0;
umin = 38;
umax = 90;
vmin = -180;
vmax = 180;
uk = 90;
vk = 0;
[XM,YM,XP,YP] = my_graticule(umin, umax,vmin, vmax, D_u, D_v, d_u, d_v, Rn, proj, uk, vk);
plot (XM', YM', 'k'); 
plot (XP', YP', 'k');
title('Tissotova elipsa na vrchní plošce dvanáctistěnu')
