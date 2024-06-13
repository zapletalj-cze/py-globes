function [] = continents(C, R, uk, vk, proj)

%Loading coordinates
con = load(C);
u_con = con(:,1);
v_con = con(:,2);

%uv to sd
[s_con, d_con] = uv_to_sd(u_con, v_con, uk, vk);

%finding and deleting singularities
s_min = 5;
idx = find(s_con<s_min);

s_con(idx) = [];
d_con(idx) = [];

%projection
[x_con,y_con] = proj(R, s_con, d_con);

%plotting continent
plot(x_con, y_con, 'b', 'LineWidth', 1);

end