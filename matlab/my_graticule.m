function [XM, YM, XP, YP] = my_graticule(umin, umax, vmin, vmax, Du, Dv, du, dv, R, proj, uk, vk)

%parallels
XP=[];YP=[];
for u=umin:Du:umax
    %Create parallel
    vp = vmin:dv:vmax;
    up = u*ones(1,length(vp));
    
    %Oblique aspect
    [sp,dp] = uv_to_sd(up,vp,uk,vk); 

    %Project parallel
    [xp,yp] = proj(R,sp,dp);
    
    %Append row
    XP = [XP;xp];
    YP = [YP;yp];
end

%meridians
XM=[];YM=[];
for v=vmin:Dv:vmax
    %Create meridian
    um = umin:du:umax;
    vm = v*ones(1,length(um));
    
    %Oblique aspect
    [sm,dm] = uv_to_sd(um,vm,uk,vk);
    
    %project meridian
    [xm,ym] = proj(R,sm,dm);
    
    %Append row
    XM=[XM;xm];
    YM=[YM;ym];
end
end