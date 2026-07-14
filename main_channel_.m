% =========================================================================
% main_channel.m -- Assignment 4a, presentation version (Channel flow)
%
% GOAL OF THIS SCRIPT
%   1) Solve the 2D unsteady Navier-Stokes equations for the channel-flow
%      case using the projection method (RK time integration + direct
%      pressure-Poisson solve, central differences on a collocated grid).
%   2) Demonstrate and quantify the 2*dx / 2*dy checkerboard oscillations
%      that arise from CDS pressure-velocity coupling on a collocated grid.
%
% WHY CHECKERBOARD MODES EXIST
%   The CDS first derivative operator is
%       dp/dx|_i = (p(i+1,j) - p(i-1,j)) / (2*dx).
%   For a grid-scale alternating mode, e.g. p'(i,j) = (-1)^i,
%       p'(i+1,j) = p'(i-1,j),
%   so its CDS gradient is identically zero. The analogous y-alternating
%   and x-y alternating modes are also invisible to the corresponding CDS
%   pressure gradients. This is the odd-even decoupling problem.
%
% DIAGNOSTICS
%   (a) Local four-neighbour high-pass indicator
%       p_cb = p - mean(p_E,p_W,p_N,p_S)
%       This is used to show WHERE grid-scale oscillations occur. It is not
%       interpreted as an amplitude-normalised modal coefficient.
%
%   (b) Projections onto the three odd-even patterns
%       (-1)^i, (-1)^j, and (-1)^(i+j).
%       These quantify x-striped, y-striped, and full checkerboard modes.
%
% PLOTTING
%   Every field panel uses its own robust symmetric colour range. This
%   avoids hiding weak checkerboard patterns behind the largest snapshot.
%   The actual colour limit is printed in each panel title so amplitudes
%   are not confused when comparing panels.
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
[parm,flow] = set_params(parm,flow,'infile_channel.mat');
[parm,flow] = initialize(parm,flow);

outdir = 'channel_presentation';
if ~exist(outdir,'dir')
    mkdir(outdir);
end

x = linspace(0,parm.xl,parm.m);
y = linspace(0,parm.yl,parm.n);

if mod(parm.m,2) ~= 0
    warning(['parm.m (periodic x-direction) is odd. The (-1)^i mode is ' ...
        'then not exactly consistent with the periodic wrap-around, ' ...
        'so its projection should be read as an approximate indicator.']);
end

% Four approximately equally distributed presentation snapshots.
track_steps = unique(max(1,min(parm.ntst, ...
    round([0.10 0.30 0.60 1.00]*parm.ntst))));
nsnap = numel(track_steps);

p_snap  = cell(nsnap,1);
cb_snap = cell(nsnap,1);
snap_t  = zeros(nsnap,1);

% Time histories.
t          = zeros(parm.ntst,1);
ucentre    = zeros(parm.ntst,1);
cb_rms     = zeros(parm.ntst,1);  % RMS of local high-pass indicator
p_rms      = zeros(parm.ntst,1);  % RMS of zero-mean raw pressure
cb_proj_x  = zeros(parm.ntst,1);  % projection onto (-1)^i
cb_proj_y  = zeros(parm.ntst,1);  % projection onto (-1)^j
cb_proj_xy = zeros(parm.ntst,1);  % projection onto (-1)^(i+j)
cb_total   = zeros(parm.ntst,1);  % combined modal amplitude

fprintf('Running channel case ...\n');

%% Time integration
for itst = 1:parm.ntst
    [flow.rhsu,flow.rhsv] = rhs_ns(parm,flow);
    flow = runge_kutta_2d_vec(parm,flow);
    flow = direct_press_corr(parm,flow);
    flow = project(parm,flow);

    % Checkerboard diagnostics.
    pcb = checkerboard_indicator(parm,flow.p);

    [cb_proj_x(itst),cb_proj_y(itst),cb_proj_xy(itst)] = ...
        checkerboard_projection(flow.p);

    cb_total(itst) = sqrt(cb_proj_x(itst)^2 + ...
                          cb_proj_y(itst)^2 + ...
                          cb_proj_xy(itst)^2);

    t(itst)       = itst*parm.dt;
    ucentre(itst) = flow.u(round(parm.m/2),round(parm.n/2));
    cb_rms(itst)  = rms_field(pcb);
    p_rms(itst)   = rms_field(remove_mean(flow.p));

    k = find(track_steps == itst,1);
    if ~isempty(k)
        p_snap{k}  = flow.p;
        cb_snap{k} = pcb;
        snap_t(k)  = t(itst);
    end
end

%% Figure 1: physical channel solution / validation
f1 = figure('Color','w','Position',[100 100 1050 410]);
tl1 = tiledlayout(f1,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl1);
plot(ax,t,ucentre,'k-');
grid(ax,'on');
box(ax,'on');
xlabel(ax,'Time [s]');
ylabel(ax,'Centre velocity U [m/s]');
title(ax,'Channel-flow convergence to steady state');
axis(ax,'tight');

% Analytical Poiseuille profile.
u_exact = parm.grav/(2*parm.nu)*y.*(parm.yl-y);

ax = nexttile(tl1);
plot(ax,u_exact,y,'k-','DisplayName','Analytical Poiseuille');
hold(ax,'on');
plot(ax,flow.u(round(parm.m/2),:),y,'ro', ...
    'MarkerFaceColor','r','DisplayName','Simulation (4a)');
grid(ax,'on');
box(ax,'on');
xlabel(ax,'U [m/s]');
ylabel(ax,'y [m]');
title(ax,'Velocity-profile validation');
legend(ax,'Location','best');
axis(ax,'tight');

title(tl1,'Channel-flow solution and validation');
exportgraphics(f1,fullfile(outdir,'Fig1_Channel_Validation.png'), ...
    'Resolution',300);

%% Figure 2: raw pressure evolution
f2 = figure('Color','w','Position',[130 80 1120 750]);
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

    % Independent robust symmetric colour range for this panel.
    c = robust_symmetric_limit(pk,99.5);
    set_symmetric_clim(ax,c);

    xlabel(ax,'x [m]');
    ylabel(ax,'y [m]');
    title(ax,sprintf('t = %.3g s, colour range = +/-%.2e',snap_t(k),c));
end

title(tl2,{ ...
    'Zero-mean raw pressure field', ...
    'Independent symmetric colour scale in each panel'});

exportgraphics(f2,fullfile(outdir,'Fig2_RawPressure_Evolution.png'), ...
    'Resolution',300);

%% Figure 3: local high-pass checkerboard indicator
f3 = figure('Color','w','Position',[160 80 1120 750]);
tl3 = tiledlayout(f3,2,2,'TileSpacing','compact','Padding','compact');

for k = 1:nsnap
    ax = nexttile(tl3);
    q = cb_snap{k};

    h = imagesc(ax,x,y,q');
    set(ax,'YDir','normal');
    set(h,'AlphaData',isfinite(q'));
    axis(ax,'tight');
    box(ax,'on');
    colormap(ax,bluewhitered(256));
    colorbar(ax);

    % Independent robust symmetric colour range for this panel.
    c = robust_symmetric_limit(q,99.5);
    set_symmetric_clim(ax,c);

    xlabel(ax,'x [m]');
    ylabel(ax,'y [m]');
    title(ax,sprintf('t = %.3g s, colour range = +/-%.2e',snap_t(k),c));
end

title(tl3,{ ...
    'Local grid-scale odd-even pressure indicator', ...
    'Independent colour scale; non-periodic wall boundaries excluded'});

exportgraphics(f3,fullfile(outdir,'Fig3_Checkerboard_Indicator.png'), ...
    'Resolution',300);

%% Figure 4: odd-even signature along automatically selected cuts
% Select the snapshot with the largest RMS high-pass indicator.
snap_strength = cellfun(@rms_field,cb_snap);
[~,kpeak] = max(snap_strength);

pcut  = remove_mean(p_snap{kpeak});
cbcut = cb_snap{kpeak};

% Find the horizontal and vertical cuts with the strongest indicator.
qsearch = cbcut;
qsearch(~isfinite(qsearch)) = 0;

row_strength = sqrt(mean(qsearch.^2,1)); % one value for each j
col_strength = sqrt(mean(qsearch.^2,2)); % one value for each i

[~,jc] = max(row_strength);
[~,ic] = max(col_strength);

f4 = figure('Color','w','Position',[190 70 1120 760]);
tl4 = tiledlayout(f4,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl4);
plot(ax,1:parm.m,pcut(:,jc),'ko-','MarkerSize',4, ...
    'MarkerFaceColor','w');
grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Grid index i');
ylabel(ax,'p - mean(p)');
title(ax,sprintf('Raw pressure along j = %d',jc));

ax = nexttile(tl4);
stem(ax,1:parm.m,cbcut(:,jc),'filled','MarkerSize',3);
yline(ax,0,'k:','HandleVisibility','off');
grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Grid index i');
ylabel(ax,'High-pass indicator');
title(ax,'Alternating signature in x direction');

ax = nexttile(tl4);
plot(ax,1:parm.n,pcut(ic,:),'ko-','MarkerSize',4, ...
    'MarkerFaceColor','w');
grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Grid index j');
ylabel(ax,'p - mean(p)');
title(ax,sprintf('Raw pressure along i = %d',ic));

ax = nexttile(tl4);
stem(ax,1:parm.n,cbcut(ic,:),'filled','MarkerSize',3);
yline(ax,0,'k:','HandleVisibility','off');
grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Grid index j');
ylabel(ax,'High-pass indicator');
title(ax,'Alternating signature in y direction');

title(tl4,sprintf( ...
    'Grid-scale odd-even pressure signature at t = %.3g s',snap_t(kpeak)));

exportgraphics(f4,fullfile(outdir,'Fig4_OddEven_Signature.png'), ...
    'Resolution',300);

%% Figure 5: quantitative evolution of checkerboard modes
f5 = figure('Color','w','Position',[230 140 980 500]);
ax = axes(f5);

semilogy(ax,t,max(abs(cb_proj_x),eps),'b-', ...
    'DisplayName','|Projection onto (-1)^i|');
hold(ax,'on');

semilogy(ax,t,max(abs(cb_proj_y),eps),'-', ...
    'Color',[0 0.55 0], ...
    'DisplayName','|Projection onto (-1)^j|');

semilogy(ax,t,max(abs(cb_proj_xy),eps),'m--', ...
    'DisplayName','|Projection onto (-1)^{i+j}|');

semilogy(ax,t,max(cb_total,eps),'r-', ...
    'LineWidth',2.2, ...
    'DisplayName','Combined odd-even modal amplitude');

semilogy(ax,t,max(cb_rms,eps),'-.', ...
    'Color',[0.85 0.40 0], ...
    'DisplayName','RMS local high-pass indicator');

semilogy(ax,t,max(p_rms,eps),'k:', ...
    'DisplayName','RMS zero-mean raw pressure');

grid(ax,'on');
box(ax,'on');
axis(ax,'tight');
xlabel(ax,'Time [s]');
ylabel(ax,'Pressure scale');
title(ax,'Evolution of CDS odd-even pressure modes');
legend(ax,'Location','best','NumColumns',2);

exportgraphics(f5,fullfile(outdir,'Fig5_Checkerboard_Amplitude.png'), ...
    'Resolution',300);

%% Final diagnostic output
fprintf('\nFinal checkerboard diagnostics:\n');
fprintf('  x-odd-even projection, (-1)^i:       %e\n',cb_proj_x(end));
fprintf('  y-odd-even projection, (-1)^j:       %e\n',cb_proj_y(end));
fprintf('  xy-checkerboard projection:          %e\n',cb_proj_xy(end));
fprintf('  combined odd-even modal amplitude:   %e\n',cb_total(end));
fprintf('  RMS local high-pass indicator:       %e\n',cb_rms(end));
fprintf('  RMS zero-mean raw pressure:          %e\n',p_rms(end));
fprintf('\nFigures saved in: %s\n',outdir);

% ============================= Helper Functions =========================

function [pcb,psmooth] = checkerboard_indicator(parm,p)
    % Local four-neighbour high-pass indicator.
    %
    % For a periodic direction, neighbour indices from parm are retained.
    % For the wall-normal direction, only interior points are evaluated;
    % wall values are set to NaN to prevent boundary treatment from being
    % mistaken for checkerboard pressure noise.

    psmooth = nan(size(p));
    pcb     = nan(size(p));

    if parm.m < 1 || parm.n < 3
        return;
    end

    % x neighbours use the solver's indexing, preserving its periodicity.
    for i = 1:parm.m
        for j = 2:parm.n-1
            psmooth(i,j) = 0.25*( ...
                p(parm.ip(i),j) + p(parm.im(i),j) + ...
                p(i,j+1) + p(i,j-1));
        end
    end

    pcb(:,2:parm.n-1) = p(:,2:parm.n-1) - psmooth(:,2:parm.n-1);
end

function [amp_x,amp_y,amp_xy] = checkerboard_projection(p)
    % Projection/correlation with the three CDS odd-even pressure modes:
    %   mode_x  = (-1)^i
    %   mode_y  = (-1)^j
    %   mode_xy = (-1)^(i+j)
    %
    % Pressure and masks are mean-corrected so the diagnostic remains
    % robust if a grid dimension is odd.

    p0 = remove_mean(p);
    [m,n] = size(p0);
    [I,J] = ndgrid(1:m,1:n);

    mask_x  = (-1).^I;
    mask_y  = (-1).^J;
    mask_xy = (-1).^(I+J);

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

function q = remove_mean(q)
    finite_values = q(isfinite(q));

    if isempty(finite_values)
        return;
    end

    q = q - mean(finite_values);
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
    % Toolbox-independent percentile-based symmetric colour limit.
    % This suppresses isolated extrema without fixing all panels to the
    % global maximum of the entire simulation.

    if nargin < 2
        percentile_value = 99.5;
    end

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

function set_symmetric_clim(ax,c)
    % Compatibility helper for MATLAB releases with/without clim().
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

    x = linspace(0,1,n)';
    cmap = [min(1,2*x), ...
            1-abs(2*x-1), ...
            min(1,2-2*x)];
end
