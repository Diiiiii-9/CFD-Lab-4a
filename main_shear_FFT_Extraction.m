% main_shear_FFT_Extraction.m
% 频域提取法：利用 2D FFT 完美分离低频物理漩涡与高频棋盘噪声

close all; clc; clear;
set(groot,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial', ...
    'DefaultAxesFontSize',12,'DefaultAxesLineWidth',1.1,'DefaultLineLineWidth',1.6);

[parm,flow] = build_structs;
[parm,flow] = set_params(parm,flow,'infile_shear.mat');
[parm,flow] = initialize(parm,flow);

target_step = round(parm.ntst * 0.8); % 在 80% 处提取
outdir = 'shear_fft_extraction'; if ~exist(outdir,'dir'), mkdir(outdir); end
x = linspace(0,parm.xl,parm.m); y = linspace(0,parm.yl,parm.n);

fprintf('Running shear case up to step %d for FFT extraction...\n', target_step);
for itst = 1:target_step
    [flow.rhsu,flow.rhsv] = rhs_ns(parm,flow);
    flow = runge_kutta_2d_vec(parm,flow);
    flow = direct_press_corr(parm,flow);
    flow = project(parm,flow);
end

% =========================================================================
% FFT 魔法剥离开始
% =========================================================================
fprintf('\n>> Performing 2D Fast Fourier Transform...\n');

P_raw = flow.p;
P_raw = P_raw - mean(P_raw(:)); % 去均值

% 1. 转换到频域并将零频移到中心
F = fftshift(fft2(P_raw));

% 2. 构建高斯滤波器 (Gaussian Filter)
[M, N] = size(F);
[U, V] = meshgrid(1:N, 1:M);
center_U = ceil(N/2); center_V = ceil(M/2);

% D 是距离频率中心的距离，D0 是截止频率（过滤掉核心 10% 的低频物理结构）
D = sqrt((U - center_U).^2 + (V - center_V).^2);
D0 = min(M, N) * 0.1; 

% 低通掩膜 (Low-pass) 保留物理漩涡，高通掩膜 (High-pass) 截获棋盘噪声
low_pass_mask = exp(-(D.^2) / (2 * (D0^2)));
high_pass_mask = 1 - low_pass_mask;

% 3. 在频域分离信号
F_smooth = F .* low_pass_mask;
F_checker = F .* high_pass_mask;

% 4. 逆变换回空间域
P_smooth = real(ifft2(ifftshift(F_smooth)));
P_checker = real(ifft2(ifftshift(F_checker)));

% =========================================================================
% 绘图验证
% =========================================================================
cl_raw = max(abs(P_raw(:)));
cl_chk = max(abs(P_checker(:)));

f1 = figure('Color','w','Position',[100 100 1200 400]);

subplot(1,3,1);
imagesc(x,y,P_raw'); set(gca,'YDir','normal'); axis equal tight; 
colormap(gca, parula); caxis([-cl_raw cl_raw]); 
title('原始场 (Raw): 漩涡 + 噪声'); xlabel('x'); ylabel('y');

subplot(1,3,2);
imagesc(x,y,P_smooth'); set(gca,'YDir','normal'); axis equal tight; 
colormap(gca, parula); caxis([-cl_raw cl_raw]); 
title('低频 (FFT Smooth): 纯净的漩涡'); xlabel('x'); ylabel('y');

subplot(1,3,3);
imagesc(x,y,P_checker'); set(gca,'YDir','normal'); axis equal tight; 
colormap(gca, bluewhitered(256)); caxis([-cl_chk cl_chk]); 
title('高频 (FFT Checker): 完美 $2\Delta x$ 噪声'); xlabel('x'); ylabel('y'); colorbar;

exportgraphics(f1, fullfile(outdir, 'Fig_FFT_Separation_2D.png'), 'Resolution', 300);

% --- 1D Cut 证明 ---
jc = round(parm.n/2);
f2 = figure('Color','w','Position',[150 150 900 400]);

subplot(1,2,1);
plot(1:parm.m, P_raw(:,jc), 'ko-', 'DisplayName', 'Raw Pressure'); hold on;
plot(1:parm.m, P_smooth(:,jc), 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed (Physical)');
grid on; box on; legend('Location','best'); title('1D Cut: 噪声是如何附着在物理低压上的');
xlabel('Grid index i'); ylabel('Pressure');

subplot(1,2,2);
stem(1:parm.m, P_checker(:,jc), 'filled', 'MarkerSize', 4);
yline(0, 'k:'); grid on; box on; title('纯正的 \lambda = 2\Delta x 棋盘效应');
xlabel('Grid index i'); ylabel('High-Frequency Noise');

exportgraphics(f2, fullfile(outdir, 'Fig_FFT_Separation_1D.png'), 'Resolution', 300);
fprintf('>> Extraction complete! Figures saved in %s.\n\n', outdir);

function cmap=bluewhitered(n)
    x=linspace(0,1,n)'; cmap=[min(1,2*x), 1-abs(2*x-1), min(1,2-2*x)]; 
end