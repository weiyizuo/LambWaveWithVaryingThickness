% 位移和应变的符号计算
clear
syms y h x ... % geometry parameter, which is real
     lambda mu A B C omega real % material parameter, which is real
syms k alpha beta... % wavenumbers of modes, which are complex
    k_2 alpha_2 beta_2
%     x y kyt kyl...
%     omegad kd cl ct
% global fx_plus fy_plus fx_minus fy_minus Sxy_plus Syx_plus Sxx_plus Syy_plus Sxy_minus Syx_minus Sxx_minus Syy_minus 


%% displacement of modes for symm 理论声学书 P461
mode = 'sym';
switch mode
    case 'sym'
        % h = d./2
        ux = 1i.*beta.*  ( k.^2.*cos(beta.*h).*cos(alpha.*y)+(beta.^2-k.^2)./2.*cos(alpha.*h).*cos(beta.*y) );
        uy = k.*( -alpha.*beta.*cos(beta.*h).*sin(alpha.*y)+(beta.^2-k.^2)./2.*cos(alpha.*h).*sin(beta.*y) );
        ux_2 = 1i.*beta_2.*  ( k_2.^2.*cos(beta_2.*h).*cos(alpha_2.*y)+(beta_2.^2-k_2.^2)./2.*cos(alpha_2.*h).*cos(beta_2.*y) );
        uy_2 = k_2.*( -alpha_2.*beta_2.*cos(beta_2.*h).*sin(alpha_2.*y)+(beta_2.^2-k_2.^2)./2.*cos(alpha_2.*h).*sin(beta_2.*y) );
    case 'asym'
        % h = d./2
        ux = 1i.*beta.*  ( k.^2.*sin(beta.*h).*sin(alpha.*y)+(beta.^2-k.^2)./2.*sin(alpha.*h).*sin(beta.*y) );
        uy = k.*( alpha.*beta.*sin(beta.*h).*cos(alpha.*y)-(beta.^2-k.^2)./2.*sin(alpha.*h).*cos(beta.*y) );
        ux_2 = 1i.*beta_2.*  ( k_2.^2.*sin(beta_2.*h).*sin(alpha_2.*y)+(beta_2.^2-k_2.^2)./2.*sin(alpha_2.*h).*sin(beta_2.*y) );
        uy_2 = k_2.*( alpha_2.*beta_2.*sin(beta_2.*h).*cos(alpha_2.*y)-(beta_2.^2-k_2.^2)./2.*sin(alpha_2.*h).*cos(beta_2.*y) );
end
% 考虑位移随x的变化, 需要将每一项后面乘以 e.^ikx
ux = ux.*exp(1i.*k.*x); uy = uy.*exp(1i.*k.*x);
ux_2 = ux_2.*exp(1i.*k_2.*x); uy_2 = uy_2.*exp(1i.*k_2.*x);
% 速度
vx_2 = -1i*omega*2*ux_2;
vy_2 = -1i*omega*2*uy_2;
%% strain
strain11_n = diff(ux,x);
strain12_n = 1./2.*( diff(ux,y) + diff(uy,x) );
strain22_n = diff(uy,y);
%% stress
stress11_n = (lambda+2*mu)*strain11_n + lambda*strain22_n;
stress12_n = 2*mu*strain12_n;
stress22_n = (lambda+2*mu)*strain22_n + lambda*strain11_n;



%% 体积力f_i和表面应力S_ij
S_xy = B*( diff(ux,x)+diff(uy,y) )*diff(uy,x) + ...
    A/4*( diff(uy,x)*diff(ux,x) +diff(uy,y)*diff(uy,x) ) + ...
    (lambda+B)*( diff(ux,x)+diff(uy,y) )*diff(ux,y) + ...
    (mu+A/4)*( diff(ux,x)*diff(uy,x) + diff(ux,y)*diff(uy,y) ...
        + diff(ux,x)*diff(ux,y) + diff(uy,x)*diff(uy,y) ...
        + diff(ux,x)*diff(ux,y) + diff(ux,y)*diff(uy,y));

f_x = (mu+A/4)*( diff(diff(ux,y),y)*diff(ux,x) + diff(diff(uy,x),x)*diff(uy,x)...
        + diff(diff(ux,y),y)*diff(ux,x) + diff(diff(uy,x),x)*diff(ux,y) ...
        + 2*diff(diff(ux,x),y)*diff(ux,y) + 2*diff(diff(ux,y),x)*diff(uy,x)) + ...
	(lambda+mu+A/4+B)*( diff(diff(ux,x),y)*diff(ux,y) + diff(diff(uy,x),x)*diff(uy,x) ...
        + diff(diff(ux,y),x)*diff(ux,y) + diff(diff(uy,x),y)*diff(ux,x) ) + ...
    (lambda+B)*( diff(diff(ux,x),x)*diff(uy,y) + diff(diff(ux,y),y)*diff(ux,x) ) + ...
    (A/4+B)*( diff(diff(ux,y),x)*diff(uy,x) + diff(diff(uy,x),y)*diff(ux,x)...
        + diff(diff(ux,x),x)*diff(ux,y) + diff(diff(uy,x),y)*diff(uy,x) ) + ...
    (B+2*C)*( diff(diff(ux,x),x)*diff(uy,y)+ diff(diff(uy,x),y)*diff(ux,x) );

S_yy = B*( diff(ux,x)+diff(uy,y) )*diff(uy,y) + ...
    A/4*( diff(uy,x)*diff(ux,y) +diff(uy,y)*diff(uy,y) ) + ...
    (lambda+B)*( diff(ux,x)+diff(uy,y) )*diff(uy,y) + ...
    (mu+A/4)*( diff(uy,x)*diff(uy,x) + diff(uy,y)*diff(uy,y) ...
        + diff(ux,y)*diff(ux,y) + diff(uy,y)*diff(uy,y) ...
        + diff(uy,x)*diff(ux,y) + diff(uy,y)*diff(uy,y)) + ...
    (lambda/2+B/2)*( diff(ux,y)*diff(ux,y)+diff(uy,x)*diff(uy,x) ) +...
    C*2*diff(ux,x)*diff(uy,y) + ...
    B/2*( diff(ux,y)*diff(uy,x)+diff(uy,x)*diff(ux,y) );

f_y = (mu+A/4)*( diff(diff(ux,y),y)*diff(ux,y) + diff(diff(uy,x),x)*diff(uy,y)...
        + diff(diff(ux,y),y)*diff(uy,x) + diff(diff(uy,x),x)*diff(uy,y) ...
        + 2*diff(diff(uy,x),y)*diff(ux,y) + 2*diff(diff(uy,y),x)*diff(uy,x)) + ...
	(lambda+mu+A/4+B)*( diff(diff(ux,y),y)*diff(ux,y) + diff(diff(uy,y),x)*diff(uy,x) ...
        + diff(diff(ux,y),x)*diff(uy,y) + diff(diff(uy,x),y)*diff(uy,x) ) + ...
    (lambda+B)*( diff(diff(uy,x),x)*diff(uy,y) + diff(diff(uy,y),y)*diff(ux,x) ) + ...
    (A/4+B)*( diff(diff(ux,y),x)*diff(uy,y) + diff(diff(uy,x),y)*diff(ux,y)...
        + diff(diff(ux,y),x)*diff(ux,y) + diff(diff(uy,y),y)*diff(uy,x) ) + ...
    (B+2*C)*( diff(diff(ux,y),x)*diff(uy,y)+ diff(diff(uy,y),y)*diff(ux,x) );
%% 功率项
df_vol_dy = 1/2*( conj(vx_2)*f_x + conj(vy_2)*f_y );
df_vol_dy = subs(df_vol_dy,x,0);

%
fid = fopen('df_vol_dy.txt', 'wt');
char_formula = char(df_vol_dy);
% replace operator for array calculation.
char_formula = strrep(char_formula,'/','./');
char_formula = strrep(char_formula,'*','.*');
char_formula = strrep(char_formula,'^','.^');
fprintf(fid, '%s;%% auto calculated expression\n', char_formula);
fclose(fid);
%%
f_surf_y = -1/2*( conj(vx_2)*S_xy + conj(vy_2)*S_yy );
f_surf_y = subs(f_surf_y,x,0);
f_surf = subs(f_surf_y,y,h) - subs(f_surf_y,y,-h);

% f_vol_y = int(df_vol_dy,y);
% f_vol_y = subs(f_vol_y,x,0);
% f_vol = subs(f_vol_y,y,h) - subs(f_vol_y,y,-h);

%
fid = fopen('f_surf.txt', 'wt');
char_formula = char(f_surf);
% replace operator for array calculation.
char_formula = strrep(char_formula,'/','./');
char_formula = strrep(char_formula,'*','.*');
char_formula = strrep(char_formula,'^','.^');
fprintf(fid, '%s;%% auto calculated expression\n', char_formula);
fclose(fid);

% fid = fopen('f_vol.txt', 'wt');
% char_formula = char(f_vol);
% % replace operator for array calculation.
% char_formula = strrep(char_formula,'/','./');
% char_formula = strrep(char_formula,'*','.*');
% char_formula = strrep(char_formula,'^','.^');
% fprintf(fid, '%s;%% auto calculated expression\n', char_formula);
% fclose(fid);


