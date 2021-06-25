  % 能流的符号计算
clear
syms y h x ... % geometry parameter, which is real
     lambda mu A B C omega real % material parameter, which is real
syms k alpha beta % wavenumbers of modes, which are complex
%     x y kyt kyl...
%     omegad kd cl ct
% global fx_plus fy_plus fx_minus fy_minus Sxy_plus Syx_plus Sxx_plus Syy_plus Sxy_minus Syx_minus Sxx_minus Syy_minus 


%% displacement of modes for symm
mode = 'sym';
switch mode
    case 'sym'
        % h = d./2
        ux = 1i.*beta.*  ( k.^2.*cos(beta.*h).*cos(alpha.*y)+(beta.^2-k.^2)./2.*cos(alpha.*h).*cos(beta.*y) );
        uy = k.*( -alpha.*beta.*cos(beta.*h).*sin(alpha.*y)+(beta.^2-k.^2)./2.*cos(alpha.*h).*sin(beta.*y) );
    case 'asym'
        % h = d./2
        ux = 1i.*beta.*  ( k.^2.*sin(beta.*h).*sin(alpha.*y)+(beta.^2-k.^2)./2.*sin(alpha.*h).*sin(beta.*y) );
        uy = k.*( alpha.*beta.*sin(beta.*h).*cos(alpha.*y)-(beta.^2-k.^2)./2.*sin(alpha.*h).*cos(beta.*y) );
end
% 考虑位移随x的变化, 需要将每一项后面乘以 e.^ikx
ux = ux.*exp(1i.*k.*x); uy = uy.*exp(1i.*k.*x);
% 速度
vx = -1i*omega*ux;
vy = -1i*omega*uy;
%% strain
strain11_n = diff(ux,x);
strain12_n = 1./2.*( diff(ux,y) + diff(uy,x) );
strain22_n = diff(uy,y);
%% stress
stress11_n = (lambda+2*mu)*strain11_n + lambda*strain22_n;
stress12_n = 2*mu*strain12_n;
stress22_n = (lambda+2*mu)*strain22_n + lambda*strain11_n;

%% 
dP_dy = 1i*( ux.*conj(stress11_n) - conj(ux).*stress11_n + uy.*conj(stress12_n) - conj(uy).*stress12_n);

fid = fopen('dP_dy.txt', 'wt');
char_formula = char(subs(dP_dy,x,0));
% replace operator for array calculation.
char_formula = strrep(char_formula,'/','./');
char_formula = strrep(char_formula,'*','.*');
char_formula = strrep(char_formula,'^','.^');
fprintf(fid, '%s;%% auto calculated expression\n', char_formula);
fclose(fid);

%%
% 能流，坡印廷矢量，理论上有解析解，实际算不出来。
% P_with_y = int(dP_dy,y); % 基于能流的归一化系数
% P_with_y = subs(P_with_y,x,0);
% P = subs(P_with_y,y,h) - subs(P_with_y,y,-h);
% % 归一化系数
% % C = sqrt(real(P))
% 
% fid = fopen('P.txt', 'wt');
% char_formula = char(P);
% % replace operator for array calculation.
% char_formula = strrep(char_formula,'/','./');
% char_formula = strrep(char_formula,'*','.*');
% char_formula = strrep(char_formula,'^','.^');
% fprintf(fid, '%s;%% auto calculated expression\n', char_formula);
% fclose(fid);

%%














