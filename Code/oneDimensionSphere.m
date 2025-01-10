%一维球体SN方法，基于杜书华流程图求解中子输运方程

%常量 变量
R = 145.5436; %球半径(cm)
N = 8;    %SN方法近似阶数
K = 100;  %将半径R等分成K个区间，暂定100，用变量 k 取
dr = R/K; %差分步长
sigma_t = 0.050;    %总截面
sigma_s = 0.030;    %散射截面
musigma_f = 0.0225; %mu*裂变截面
epsilon = 1e-8; %判断计算上下两次迭代得到的keff判断是否达到稳定
mu = zeros(2*N+1,1);%高斯求积组初始化
omega_m = zeros(2*N+1,1);%高斯求积组初始化
%获取高斯求积元组
[mu(2:2:2*N),omega_m(2:2:2*N)] = getGaussianQuadTuple(N);
mu(1:2*N:2*N+1) = [-1,1];
mu(3:2:2*N-1) = 0.5*(mu(2:2:2*N-2)+mu(4:2:2*N));
r = 0:dr/2:R; %小区间坐标点
A = zeros(2*K+1,1); %中间变量，r^2
A(1:2:2*K+1)=r(1:2:2*K+1).^2; % 计算r_2k+-1的平方，杜书华书式9.3.3
V = zeros(2*K+1,1); %中间变量，1/3*(r^3-r^3)
V(2:2:2*K)=(r(3:2:(2*K+1)).^3-r(1:2:(2*K-1)).^3)/3; % 杜书华书式9.3.3
E = zeros(2*K,2*N); %中间变量
G = zeros(2*K,2*N); %中间变量
alpha = zeros(2*N+1,1); %杜书华书式9.3.3
for m = 1:N
   alpha(2*m+1) = alpha(2*m-1) - mu(2*m)*omega_m(2*m); %杜书华书P332(1)
end
phi = zeros(2*K+1,2*N+1); %中子角通量初值
S = 1.0*ones(2*K,1);      %源项初值
keff = 0.1*ones(2*K+1,1); %有效增殖系数
Phi = zeros(2*K+1,1);     %中子通量
old_Phi = 100*ones(2*K+1,1);%存储上一次迭代的中子通量，初始化为100
m = 0;
k = K;
freq = 0;   %迭代次数统计
signal = 1; %判断误差的指示，小于此值，结束循环

%循环计算
while (signal > 0)
    %mu = -1
    for k = K:(-1):1
        phi(2*k,1) = (V(2*k)*S(2*k)+(A(2*k+1)+A(2*k-1))*phi(2*k+1,1))/...
            (A(2*k+1)+A(2*k-1)+sigma_t*V(2*k));
        phi(2*k-1,1) = 2*phi(2*k,1) - phi(2*k+1,1);
    end
    %mu < 0
    for m = 1:N/2
         for k = K:(-1):1
            %杜书华书式9.3.19
            E(2*k,2*m) = abs(mu(2*m))*(A(2*k+1)+A(2*k-1));
            G(2*k,2*m) = (A(2*k+1)-A(2*k-1))*(alpha(2*m+1)+alpha(2*m-1))/omega_m(2*m);
            phi(2*k,2*m) = (V(2*k)*S(2*k)+E(2*k,2*m)*phi(2*k+1,2*m)+...
                G(2*k,2*m)*phi(2*k,2*m-1))/(E(2*k,2*m)+G(2*k,2*m)+sigma_t*V(2*k));
            if(phi(2*k,2*m) < 0)
                phi(2*k,2*m) = 0;
            end
            %以下两个通量式子均为，杜书华式9.3.18
            phi(2*k-1,2*m) = 2*phi(2*k,2*m) - phi(2*k+1,2*m);
            if(phi(2*k-1,2*m) < 0)
                phi(2*k-1,2*m) = 0;
            end
            phi(2*k,2*m+1) = 2*phi(2*k,2*m) - phi(2*k,2*m-1);
            if(phi(2*k,2*m+1) < 0)
                phi(2*k,2*m+1) = 0;
            end
        end
    end
     % mu > 0
    for m = (N/2+1):N
        phi(1,2*m) = phi(1,2*N+2-2*m); %球心对称条件，杜书华书9.3.8
        for k = 1:K
            %杜书华书式9.3.19
            E(2*k,2*m) = abs(mu(2*m))*(A(2*k+1)+A(2*k-1));
            G(2*k,2*m) = (A(2*k+1)-A(2*k-1))*(alpha(2*m+1)+alpha(2*m-1))/omega_m(2*m);
            phi(2*k,2*m) = (V(2*k)*S(2*k)+E(2*k,2*m)*phi(2*k-1,2*m)+...
                G(2*k,2*m)*phi(2*k,2*m-1))/(E(2*k,2*m)+G(2*k,2*m)+sigma_t*V(2*k));
            if(phi(2*k,2*m) < 0)
                phi(2*k,2*m) = 0;
            end
            %以下两个通量式子均为，杜书华式9.3.18
            phi(2*k+1,2*m) = 2*phi(2*k,2*m) - phi(2*k-1,2*m);
            if(phi(2*k+1,2*m) < 0)
                phi(2*k+1,2*m) = 0;
            end
            phi(2*k,2*m+1) = 2*phi(2*k,2*m) - phi(2*k,2*m-1);
            if(phi(2*k,2*m+1) < 0)
                phi(2*k,2*m+1) = 0;
            end
        end
    end
    %更新源项
    S(2:2:2*K) = (sigma_s + musigma_f)/2*phi(2:2:2*K,2:2:2*N)*omega_m(2:2:2*N);
    %判断是否达到稳定
    signal = 0;
    for k = 1:K
        Phi(2*k) = phi(2*k,2:2:2*N)*omega_m(2:2:2*N); %用角通量计算中子通量
        %计算当前计算面上的有效增殖系数，keff
        keff_temp = Phi(2*k)/old_Phi(2*k);
        %用上下两次迭代的keff的差判断是否达到稳定
        if(abs((keff_temp-keff(2*k))/keff(2*k)) > epsilon)
           signal = signal + 1;%不满足小于阈值的话，指示+1，继续下一次迭代
        end
        keff(2*k) = keff_temp;
    end
    old_Phi = Phi;%将本次迭代计算得到的中子通量记录，以便下一次迭代算keff
    freq = freq + 1;%统计迭代次数
end

fprintf('有效增值系数（平均值）： %.8f\n', mean(keff));
fprintf('有效增值系数（最大值）： %.8f\n', max(keff));
fprintf('有效增值系数（最小值）： %.8f\n', min(keff));
fprintf('迭代次数： %f\n', freq);
%画出中子通量分布
plot(r(2:2:2*K),Phi(2:2:2*K)/max(Phi));
xlim([0,R]);
title('一维球体SN中子通量分布')
xlabel('球半径R(cm)');
ylabel('中子通量(归一化)');
legendText = ['N = ', num2str(N)];
legend(legendText);
grid on;

