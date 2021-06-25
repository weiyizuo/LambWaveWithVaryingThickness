clear
%% 定义参数
mode = 1;
discrete_number = 1000;
amplitude = 0.01; % 位移的幅值。
len = 1;

geometry = 'uniform';
switch geometry
    case 'uniform'
        x = linspace(0,len,discrete_number); % propagation length
        d = linspace(0.005,0.005,discrete_number); % theckness
        frequency = 90e3; % circular frequency 
    case 'uphill'
        x = linspace(0,len,discrete_number); % propagation length
        d = linspace(0.003,0.005,discrete_number); % theckness
        frequency = 90e3; % circular frequency     
    case 'downhill'
        x = linspace(0,len,discrete_number); % propagation length
        d = linspace(0.005,0.003,discrete_number); % theckness
        frequency = 90e3; % circular frequency 
    case 'gauss'
        x = linspace(0,len,discrete_number); % propagation length
        d = 0.005-0.0015*exp(-100*(x-len*0.3).^2); % theckness
        frequency = 90e3; % circular frequency 
end
% frequency = 90e3; % circular frequency 
% material
c_l_real = 2740;
c_t_real = 1400;
gamma_l = 0.0035; % Np/m/Hz
gamma_t = 0.0053; % Np/m/Hz
c_l = c_l_real/(1+1i*gamma_l);
c_t = c_t_real/(1+1i*gamma_t);
lambda = 4.25e9;
mu = 2.32e9;
A = -61e9;
B = -3.5e9; %-34e9 - 1/2*A;
C = -69.5e9; %-73e9 - B;

%% 考虑/不考虑衰减的结果
for if_damping = 0:1
    % 加载频散曲线
    if if_damping == 1
        load('PMMA_s0.mat')
        c_l = c_l_real/(1+1i*gamma_l);
        c_t = c_t_real/(1+1i*gamma_t);
    else
        load('PMMA_s0_real.mat')
        c_l = c_l_real;
        c_t = c_t_real;
    end
    %
    h = d/2;
    omega = frequency*2*pi; % angular frequency
    od = omega*d;

    for j = length(kd_r_modes):-1:1
        if (isnan(kd_r_modes(mode,j)) && isnan(kd_i_modes(mode,j)) && isnan(od_modes(mode,j)))
            kd_r_modes(:,j) = [];
            kd_i_modes(:,j) = [];
            od_modes(:,j) = [];
        end
    end
%     kd_i_modes_if = kd_i_modes*if_damping;
    kd_r = interp1(od_modes(mode,:),kd_r_modes(mode,:),od,'spline');
    kd_i = interp1(od_modes(mode,:),kd_i_modes(mode,:),od,'spline');
    k_r = kd_r./d;
    k_i = kd_i./d;
    k = k_r + 1i*k_i;

    % 二倍频
    kd_r_2 = interp1(od_modes(mode,:),kd_r_modes(mode,:),od*2,'spline');
    kd_i_2 = interp1(od_modes(mode,:),kd_i_modes(mode,:),od*2,'spline');
    k_r_2 = kd_r_2./d;
    k_i_2 = kd_i_2./d;
    k_2 = k_r_2 + 1i*k_i_2;

    P =   P_energy_flow(k,  h, omega,  lambda,mu,c_l,c_t); % 能流是位移的二次项。
    P_2 = P_energy_flow(k_2,h, omega*2,lambda,mu,c_l,c_t);
    [f_surf,f_vol] = f_power(k_2,k,h, omega,lambda,mu,c_l,c_t,A,B,C);
    %
    uniform_P_2 = P_2./real(P_2); % 对能流归一化。
    uniform_f_surf = amplitude^2*f_surf./real(P)./sqrt(real(P_2)); % 基频位移的二次项，二倍频位移的一次项，分别进行归一化。
    uniform_f_vol =  amplitude^2*f_vol ./real(P)./sqrt(real(P_2));

    % 
    Qx = (uniform_f_surf+uniform_f_vol)./uniform_P_2 .* exp(1i*cumtrapz(x,2*k)); % 数值积分
    Px = -1i*k_2;
    Fx = cumtrapz(x,Px);
    C_1 = 0;%-Qx(1);
    if if_damping == 1
        A1 = C_1*exp(-Fx) + exp(-Fx).*cumtrapz(x,Qx.*exp(Fx));
    else
        A2 = C_1*exp(-Fx) + exp(-Fx).*cumtrapz(x,Qx.*exp(Fx));
    end
end

%%
figure
hold on
plot(x,abs(A1),'LineWidth',2)
plot(x,abs(A2))%,'--')
grid on
legend('衰减下的谐波','不考虑衰减下的谐波')
xlabel('传播距离 (m)')
ylabel('谐波振幅(m)')
title(geometry)







