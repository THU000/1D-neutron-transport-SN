%一维平板SN方法，基于杜书华流程图求解中子输运方程
%左边为对称边界，右边为真空边界

%常量 变量
X = 66.0053; %平板的半厚度，单位cm
K = 100;     %将半厚度X等分成K个区间，暂定100，用变量 k 取
N = 4;       %SN方法近似阶数，偶数，可以修改
dx = X/K;    %差分步长
x = (dx/2):dx:(X-dx/2); %从x_1到x_K，即从0到X的K个区间的整数中点坐标
sigma_t = 0.05;     %总截面
sigma_s = 0.03;     %散射截面
muSigma_f = 0.0235; %mu*裂变截面
%高斯求积组，关于0对称
[mu_m,omega_m] = getGaussianQuadTuple(N);
L = zeros(N,1);       %杜书华9.2.8(a)式中间变量，初始化
phi = zeros(2*K+1,N); %中子角通量，前两列为mu<0的，后两列为mu>0的
S = 1*ones(2*K+1,1);  %中子源，迭代初值暂定1
keff = 0.1*ones(2*K+1,1); %有效增殖系数keff
epsilon = 1e-8; %判断计算上下两次迭代得到的keff判断是否达到稳定
Phi = zeros(2*K+1,1); %迭代计算得到的中子通量，初始化为0
old_Phi = 100*ones(2*K+1,1); %存储上一次迭代的中子通量，初始化为100
freq = 0;   %迭代次数统计
signal = 1; %判断误差的指示，小于此值，结束循环

%循环计算
while(signal > 0)
    for m = 1:(N/2) %真空边界条件计算，mu < 0
        phi(2*K+1,m) = 0; %右边真空边界条件，循环之前常量赋值已经将 k 初始化为 K
        L(m) = sigma_t*dx/mu_m(m); %杜书华9.2.8(a)式中间变量
        for k = K:(-1):1
            %杜书华9.2.8(a)式，计算k-0.5的phi
            phi(2*k-1,m) = (2+L(m))*phi(2*k+1,m)/(2-L(m)) - 2*L(m)*S(2*k)...
                /sigma_t/(2-L(m));
            %整数中点位置的中子角通量，phi
            phi(2*k,m) = (phi(2*k+1,m) + phi(2*k-1,m))/2;
        end
    end
    for m = (N/2 + 1):N %对称边界条件计算，mu > 0
        phi(1,m) = phi(1,N+1-m); %左边对称边界条件
        L(m) = sigma_t*dx/mu_m(m); %杜书华9.2.8(a)式中间变量
        for k = 1:K
            %杜书华9.2.8(b)式，计算k+0.5的phi
            phi(2*k+1,m) = (2-L(m))*phi(2*k-1,m)/(2+L(m)) + 2*L(m)*S(2*k)...
                /sigma_t/(2+L(m));
            %整数中点位置的中子角通量，phi
            phi(2*k,m) = (phi(2*k+1,m) + phi(2*k-1,m))/2;
        end
    end
    %更新迭代源
    for k = 1:K
        %杜书华(9.2.9)式
        S(2*k) = (sigma_s + muSigma_f)/2*dot(omega_m, phi(2*k,:)');
        %S(2*k) = (sigma_s + muSigma_f)/2*dot(omega_m, phi(2*k,:)');
        %杜书华书上有问题，源迭代应该还要除以keff_old
    end
    signal = 0; %指示等于0，便于下面判断
    for k = 1:(2*K+1)
        %用每个计算面上的中子角通量，计算2K+1个计算面的中子通量
        Phi(k) = dot(omega_m, phi(k,:)'); 
        %计算当前计算面上的有效增殖系数，keff
        keff_temp = Phi(k)/old_Phi(k);
        %修改需要在keff_temp上乘上keff_old
        %判断是否达到稳定
        if(abs(keff_temp-keff(k))/keff(k) > epsilon)
            signal = signal + 1; %不满足小于阈值的话，指示+1，继续下一次迭代
        end
        keff(k) = keff_temp; %将本次迭代计算的keff记录
    end
    old_Phi = Phi;   %将本次迭代计算得到的中子通量记录，以便下一次迭代算keff
    freq = freq + 1; %统计迭代次数
end

fprintf('有效增值系数（平均值）： %.8f\n', mean(keff));
fprintf('有效增值系数（最大值）： %.8f\n', max(keff));
fprintf('有效增值系数（最小值）： %.8f\n', min(keff));
fprintf('迭代次数： %f\n', freq);
%画出中子通量分布
plot(0:dx/2:X,Phi/max(Phi));
% ylim([0,1.0]);
xlim([0,X]);
title('一维平板SN中子通量分布')
xlabel('半厚度板x(cm)');
ylabel('中子通量(归一化)');
legendText = ['N = ', num2str(N)];
legend(legendText);
grid on;
