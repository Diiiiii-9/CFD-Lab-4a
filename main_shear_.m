% =========================================================================
% main_shear_.m -- Assignment 4a, presentation version (Shear-layer flow)
%
% PURPOSE
%   1) Solve the periodic 2D unsteady Navier-Stokes shear-layer case using
%      the same numerical solver as main_channel.m.
%   2) demonstrate the CDS/collocated-grid odd-even pressure modes without
%      confusing them with the large physical pressure fluctuations caused
%      by the rolling-up shear layer and vortices.
%
% WHY MODAL / SPECTRAL DIAGNOSTICS ARE USED
%   A local four-neighbour filter is useful for a smooth channel pressure
%   field, but it can leak real vortex-driven pressure curvature into the
%   apparent checkerboard field. Since this shear-layer problem is periodic
%   in x and y, the grid-scale odd-even patterns can instead be isolated by
%   discrete modal projection.
%
%   The three CDS odd-even pressure modes are
%       x-Nyquist mode:   (-1)^i
%       y-Nyquist mode:   (-1)^j
%       xy-Nyquist mode:  (-1)^(i+j)
%
%   For an even periodic grid these are discrete Fourier modes and are
%   orthogonal to all other resolved Fourier modes. The corresponding
%   pressure coefficients therefore separate grid-scale oscillations from
%   the strong low-wavenumber physical pressure field.
%
% IMPORTANT FFT CONVENTION
%   After fftshift, the wavenumber axes for an even grid are
%       -m/2, ..., m/2-1  and  -n/2, ..., n/2-1.
%   Therefore the Nyquist locations appear at -m/2 and/or -n/2, not at
%   +m/2 and +n/2. Figure 3 marks all three relevant Nyquist modes.
%
% PLOTTING
%   Field panels use independent robust colour limits so weak structures
%   are not hidden by the largest snapshot. Spectrum panels use a common
%   log-power range so spectral amplitudes remain directly comparable.
% =========================================================================

close all;
clc;
clear;

set(groot, ...
    'DefaultAxesFontName','Arial', ...
    'DefaultTextFontName','Arial', ...
    'DefaultAxesFontSize',12, ...
    'DefaultAxesLineWidth',1.1, ...
    'DefaultLineLineWidth',1.6);

%% Setup
[parm,flow] = build_structs;
[parm,flow] = set_params(parm,flow,'infile_shear.mat');
[parm,flow] = initialize(parm,flow);

outdir = 'shear_presentation';
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Periodic grids must not duplicate the endpoint. Using dx and dy also
% guarantees consistency with the solver's actual mesh.
x = (0:parm.m-1)*parm.dx;
y = (0:parm.n-1)*parm.dy;
[X,Y] = meshgrid(x,y);

if mod(parm.m,2) ~= 0 || mod(parm.n,2) ~= 0
    warning(['At least one grid dimension is odd. A true discrete Nyquist ' ...
             'mode then does not exist in that direction. Modal correlations ' ...
             'are still reported, but the exact Nyquist interpretation and ' ...
             'the red FFT markers require even m and n.']);
end

track_steps = unique(max(1,min(parm.ntst, ...
    round([0.10 0.30 0.60 1.00]*parm.ntst))));
nsnap = numel(track_steps);

p_snap    = cell(nsnap,1);
spec_snap = cell(nsnap,1);
snap_t    = zeros(nsnap,1);

% Time histories.
t          = zeros(parm.ntst,1);
ke         = zeros(parm.ntst,1);
p_rms      = zeros(parm.ntst,1);
cb_proj_x  = zeros(parm.ntst,1);
cb_proj_y  = zeros(parm.ntst,1);
cb_proj_xy = zeros(parm.ntst,1);
cb_total   = zeros(parm.ntst,1);

fprintf('Running shear-layer case ...\n');

%% Time integration
for itst = 1:parm.ntst
    [flow.rhsu,flow.rhsv] = rhs_ns(parm,flow);
    flow = runge_kutta_2d_vec(parm,flow);
    flow = direct_press_corr(parm,flow);
    flow = project(parm,flow);

    % Exact odd-even modal diagnostics on an even doubly-periodic grid.
    [cb_proj_x(itst),cb_proj_y(itst),cb_proj_xy(itst)] = ...
        checkerboard_projection(flow.p);

    cb_total(itst) = sqrt(cb_proj_x(itst)^2 + ...
                          cb_proj_y(itst)^2 + ...
                          cb_proj_xy(itst)^2);

    t(itst) = itst*parm.dt;
    ke(itst) = 0.5*sum(flow.u(:).^2 + flow.v(:).^2)*parm.dx*parm.dy;
    p_rms(itst) = rms_field(remove_mean(flow.p));

    k = find(track_steps == itst,1);
    if ~isempty(k)
        p0 = remove_mean(flow.p);
        p_snap{k} = flow.p;

        % Normalised Fourier power: |P_hat/N|^2. This makes the scale
        % independent of the total number of points.
        spec_snap{k} = fftshift(abs(fft2(p0)/numel(p0)).^2);
        snap_t(k) = t(itst);
    end
end

% Snapshot modal coefficients, used by Figures 3 and 5.
cb_snap_x  = zeros(nsnap,1);
cb_snap_y  = zeros(nsnap,1);
cb_snap_xy = zeros(nsnap,1);
for k = 1:nsnap
    [cb_snap_x(k),cb_snap_y(k),cb_snap_xy(k)] = ...
        checkerboard_projection(p_snap{k});
end
cb_snap_total = sqrt(cb_snap_x.^2 + cb_snap_y.^2 + cb_snap_xy.^2);

%% Figure 1: physical shear-layer solution
% curl expects arrays arranged as rows in y and columns in x. The solver
% stores fields as (i,j), so transpose u and v for plotting/diagnostics.
[omega,~] = curl(X,Y,flow.u',flow.v');
speed = hypot(flow.u',flow.v');

f1 = figure('Color','w','Position',[80 90 1220 430]);
tl1 = tiledlayout(f1,1,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl1);
ke0 = max(abs(ke(1)),eps);
plot(ax,t,ke/ke0,'k-');
grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Time [s]');
ylabel(ax,'E_k/E_k(0)');
title(ax,'Kinetic-energy evolution');

ax = nexttile(tl1);
contourf(ax,X,Y,omega,40,'LineStyle','none');
axis(ax,'tight');
box(ax,'on');
colorbar(ax);
colormap(ax,bluewhitered(256));
c = robust_symmetric_limit(omega,99.5);
set_symmetric_clim(ax,c);
xlabel(ax,'x [m]');
ylabel(ax,'y [m]');
title(ax,sprintf('Vorticity, range = +/-%.2e',c));

ax = nexttile(tl1);
contourf(ax,X,Y,speed,40,'LineStyle','none');
hold(ax,'on');
streamslice(ax,X,Y,flow.u',flow.v');
axis(ax,'tight');
box(ax,'on');
colorbar(ax);
colormap(ax,parula(256));
xlabel(ax,'x [m]');
ylabel(ax,'y [m]');
title(ax,'Velocity magnitude and streamlines');

title(tl1,'Periodic shear-layer flow');
exportgraphics(f1,fullfile(outdir,'Fig1_Shear_Flow.png'),'Resolution',300);

%% Figure 2: raw pressure evolution
f2 = figure('Color','w','Position',[120 70 1120 750]);
tl2 = tiledlayout(f2,2,2,'TileSpacing','compact','Padding','compact');

for k = 1:nsnap
    ax = nexttile(tl2);
    pk = remove_mean(p_snap{k});

    imagesc(ax,x,y,pk');
    set(ax,'YDir','normal');
    axis(ax,'tight');
    box(ax,'on');
    colormap(ax,parula(256));
    colorbar(ax);

    % Independent colour range: best for seeing the spatial structure at
    % every time. The numerical range is printed to prevent miscomparison.
    c = robust_symmetric_limit(pk,99.5);
    set_symmetric_clim(ax,c);

    xlabel(ax,'x [m]');
    ylabel(ax,'y [m]');
    title(ax,sprintf('t = %.3g s, range = +/-%.2e',snap_t(k),c));
end

title(tl2,{ ...
    'Zero-mean raw pressure: physical vortex structures dominate', ...
    'Each panel uses an independent symmetric colour scale'});

exportgraphics(f2,fullfile(outdir,'Fig2_RawPressure_Evolution.png'), ...
    'Resolution',300);

%% Figure 3: 2D pressure power spectrum
% fftshift ordering. This definition works for both even and odd sizes.
kx = -floor(parm.m/2):ceil(parm.m/2)-1;
ky = -floor(parm.n/2):ceil(parm.n/2)-1;

% Common log-power scale across snapshots for honest time comparison.
all_log_power = [];
for k = 1:nsnap
    all_log_power = [all_log_power; log10(spec_snap{k}(:) + realmin)]; %#ok<AGROW>
end
spec_clim = robust_limits(all_log_power,0.5,99.8);

f3 = figure('Color','w','Position',[150 60 1160 780]);
tl3 = tiledlayout(f3,2,2,'TileSpacing','compact','Padding','compact');

for k = 1:nsnap
    ax = nexttile(tl3);
    log_power = log10(spec_snap{k}' + realmin);

    imagesc(ax,kx,ky,log_power);
    set(ax,'YDir','normal');
    axis(ax,'tight');
    box(ax,'on');
    colormap(ax,parula(256));
    colorbar(ax);
    caxis(ax,spec_clim);
    hold(ax,'on');

    % With fftshift, even-grid Nyquist bins are located on the negative
    % edge. Mark x-Nyquist, y-Nyquist and xy-Nyquist modes separately.
    if mod(parm.m,2) == 0
        plot(ax,-parm.m/2,0,'rx','MarkerSize',11,'LineWidth',2);
    end
    if mod(parm.n,2) == 0
        plot(ax,0,-parm.n/2,'r^','MarkerSize',8,'LineWidth',1.8);
    end
    if mod(parm.m,2) == 0 && mod(parm.n,2) == 0
        plot(ax,-parm.m/2,-parm.n/2,'r+', ...
            'MarkerSize',13,'LineWidth',2.2);
    end

    xlabel(ax,'k_x [DFT index]');
    ylabel(ax,'k_y [DFT index]');
    title(ax,sprintf(['t = %.3g s;  |A_x|=%.1e, |A_y|=%.1e, ' ...
                      '|A_{xy}|=%.1e'], ...
        snap_t(k),abs(cb_snap_x(k)),abs(cb_snap_y(k)),abs(cb_snap_xy(k))));
end

title(tl3,{ ...
    'Normalised pressure power spectrum: log_{10}|P_hat/N|^2', ...
    'red x: x-Nyquist; red triangle: y-Nyquist; red +: xy-Nyquist'});

exportgraphics(f3,fullfile(outdir,'Fig3_Pressure_PowerSpectrum.png'), ...
    'Resolution',300);

%% Figure 4: quantitative modal evolution
f4 = figure('Color','w','Position',[210 130 1000 510]);
ax = axes(f4);

semilogy(ax,t,max(abs(cb_proj_x),eps),'b-', ...
    'DisplayName','|x-Nyquist projection, (-1)^i|');
hold(ax,'on');
semilogy(ax,t,max(abs(cb_proj_y),eps),'-', ...
    'Color',[0 0.55 0], ...
    'DisplayName','|y-Nyquist projection, (-1)^j|');
semilogy(ax,t,max(abs(cb_proj_xy),eps),'m--', ...
    'DisplayName','|xy-Nyquist projection, (-1)^{i+j}|');
semilogy(ax,t,max(cb_total,eps),'r-', ...
    'LineWidth',2.2, ...
    'DisplayName','Combined odd-even modal amplitude');
semilogy(ax,t,max(p_rms,eps),'k:', ...
    'DisplayName','RMS zero-mean raw pressure');

grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Time [s]');
ylabel(ax,'Pressure scale');
title(ax,'Odd-even pressure modes versus physical pressure fluctuations');
legend(ax,'Location','best','NumColumns',2);

exportgraphics(f4,fullfile(outdir,'Fig4_Checkerboard_Amplitude.png'), ...
    'Resolution',300);

%% Figure 5: isolate the strongest odd-even mode in physical space
% A raw pressure cut is vortex dominated and may conceal the alternating
% signal. Reconstructing the measured odd-even modes shows the CDS mode
% directly, while the adjacent raw-pressure plots retain physical context.
[~,kpeak] = max(cb_snap_total);
ppeak = remove_mean(p_snap{kpeak});
pcb_peak = reconstruct_checkerboard(size(ppeak), ...
    cb_snap_x(kpeak),cb_snap_y(kpeak),cb_snap_xy(kpeak));

% Select cuts where the reconstructed checkerboard has the largest RMS.
row_strength = sqrt(mean(pcb_peak.^2,1));
col_strength = sqrt(mean(pcb_peak.^2,2));
[~,jc] = max(row_strength);
[~,ic] = max(col_strength);

f5 = figure('Color','w','Position',[180 60 1180 780]);
tl5 = tiledlayout(f5,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl5);
plot(ax,1:parm.m,ppeak(:,jc),'k.-','MarkerSize',11);
grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Grid index i');
ylabel(ax,'p - mean(p)');
title(ax,sprintf('Raw pressure along j = %d',jc));

ax = nexttile(tl5);
stem(ax,1:parm.m,pcb_peak(:,jc),'filled','MarkerSize',3);
yline(ax,0,'k:','HandleVisibility','off');
grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Grid index i');
ylabel(ax,'Reconstructed odd-even pressure');
title(ax,'Isolated alternating signal in x');

ax = nexttile(tl5);
plot(ax,1:parm.n,ppeak(ic,:),'k.-','MarkerSize',11);
grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Grid index j');
ylabel(ax,'p - mean(p)');
title(ax,sprintf('Raw pressure along i = %d',ic));

ax = nexttile(tl5);
stem(ax,1:parm.n,pcb_peak(ic,:),'filled','MarkerSize',3);
yline(ax,0,'k:','HandleVisibility','off');
grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Grid index j');
ylabel(ax,'Reconstructed odd-even pressure');
title(ax,'Isolated alternating signal in y');

title(tl5,sprintf( ...
    'Strongest tracked odd-even signal at t = %.3g s',snap_t(kpeak)));

exportgraphics(f5,fullfile(outdir,'Fig5_OddEven_PressureCuts.png'), ...
    'Resolution',300);

%% Final diagnostic output
fprintf('\nFinal shear-layer pressure diagnostics:\n');
fprintf('  x-Nyquist projection, (-1)^i:       %e\n',cb_proj_x(end));
fprintf('  y-Nyquist projection, (-1)^j:       %e\n',cb_proj_y(end));
fprintf('  xy-Nyquist projection:              %e\n',cb_proj_xy(end));
fprintf('  combined odd-even modal amplitude:  %e\n',cb_total(end));
fprintf('  RMS zero-mean raw pressure:         %e\n',p_rms(end));
fprintf('\nFigures saved in: %s\n',outdir);

% ============================= Helper Functions =========================

function [amp_x,amp_y,amp_xy] = checkerboard_projection(p)
    % Modal coefficients of the three CDS odd-even patterns. For an even,
    % doubly periodic grid these are exact discrete Fourier coefficients.

    p0 = remove_mean(p);
    [m,n] = size(p0);
    [I,J] = ndgrid(0:m-1,0:n-1);

    mask_x  = (-1).^I;
    mask_y  = (-1).^J;
    mask_xy = (-1).^(I+J);

    % Mean correction keeps the correlation well defined for odd grids.
    % On an even grid the masks already have exactly zero mean.
    mask_x  = mask_x  - mean(mask_x(:));
    mask_y  = mask_y  - mean(mask_y(:));
    mask_xy = mask_xy - mean(mask_xy(:));

    amp_x  = normalized_projection(p0,mask_x);
    amp_y  = normalized_projection(p0,mask_y);
    amp_xy = normalized_projection(p0,mask_xy);
end

function amp = normalized_projection(p,mask)
    denom = sum(mask(:).^2);
    if isfinite(denom) && denom > eps
        amp = sum(p(:).*mask(:))/denom;
    else
        amp = 0;
    end
end

function pcb = reconstruct_checkerboard(field_size, amp_x, amp_y, amp_xy)
    m = field_size(1); n = field_size(2);
    [I,J] = ndgrid(0:m-1, 0:n-1);

    mask_x  = (-1).^I;
    mask_y  = (-1).^J;
    mask_xy = (-1).^(I+J);

    % constant offset appears whenever m or n is odd
    mask_x  = mask_x  - mean(mask_x(:));
    mask_y  = mask_y  - mean(mask_y(:));
    mask_xy = mask_xy - mean(mask_xy(:));

    pcb = amp_x*mask_x + amp_y*mask_y + amp_xy*mask_xy;
end

function q = remove_mean(q)
    values = q(isfinite(q));
    if ~isempty(values)
        q = q - mean(values);
    end
end

function r = rms_field(q)
    q = q(isfinite(q));
    if isempty(q)
        r = 0;
    else
        r = sqrt(mean(q.^2));
    end
end

function lim = robust_symmetric_limit(q,percentile_value)
    q = abs(q(isfinite(q)));
    if isempty(q)
        lim = 1;
        return;
    end

    q = sort(q(:));
    percentile_value = min(max(percentile_value,0),100);
    idx = 1 + round((numel(q)-1)*percentile_value/100);
    idx = min(max(idx,1),numel(q));
    lim = q(idx);

    if ~isfinite(lim) || lim <= eps
        lim = max(q);
    end
    if ~isfinite(lim) || lim <= eps
        lim = 1;
    end
end

function limits = robust_limits(q,lower_percentile,upper_percentile)
    q = sort(q(isfinite(q)));
    if isempty(q)
        limits = [-16 0];
        return;
    end

    n = numel(q);
    i1 = 1 + round((n-1)*lower_percentile/100);
    i2 = 1 + round((n-1)*upper_percentile/100);
    i1 = min(max(i1,1),n);
    i2 = min(max(i2,1),n);
    limits = [q(i1) q(i2)];

    if limits(2) <= limits(1)
        limits = [limits(1)-1 limits(1)+1];
    end
end

function set_symmetric_clim(ax,c)
    if exist('clim','file') || exist('clim','builtin')
        clim(ax,[-c c]);
    else
        caxis(ax,[-c c]);
    end
end

function cmap = bluewhitered(n)
    if nargin < 1
        n = 256;
    end
    s = linspace(0,1,n)';
    cmap = [min(1,2*s), 1-abs(2*s-1), min(1,2-2*s)];
end
