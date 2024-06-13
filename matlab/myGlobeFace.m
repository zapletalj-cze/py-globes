function[]= myGlobeFace(ub, vb, umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk)

axis equal
hold on

%Create graticule
[XM, YM, XP, YP] = my_graticule(umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk);

%Draw graticule
plot(XM',YM','k', 'LineWidth', 0.5);
plot(XP',YP','k', 'LineWidth', 0.5);

%Continents
continents('eur.txt', R, uk, vk, proj);
continents('amer.txt', R, uk, vk, proj);
continents('austr.txt', R, uk, vk, proj);
continents('anta.txt', R, uk, vk, proj);

%Set limits for drawing
cmin = -2*R;
cmax = -cmin;
xlim([cmin, cmax]); ylim([cmin, cmax]); 

%Boundary
[sb, db] = uv_to_sd(ub, vb, uk, vk);
[xb, yb] = proj(R, sb, db);
plot(xb, yb, 'r');

end