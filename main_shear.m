% main_shear.m -- Assignment 4a presentation version
% Removes smooth pressure variations and highlights pressure that changes 
% sign from one collocated node to the next.

close all; clc; clear;
set(groot,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial', ...
    'DefaultAxesFontSize',12,'DefaultAxesLineWidth',1.1,'DefaultLineLineWidth',1.6);

[parm,flow] = build_structs;
[parm,flow] = set_params(parm,flow,'infile_shear.mat');
[parm,flow] = initialize(parm,flow);

outdir = 'shear_presentation'; if ~exist(outdir,'dir'), mkdir(outdir); end
x = linspace(0,parm.xl,parm.m); y = linspace(0,parm.yl,parm.n); [X,Y] = meshgrid(x,y);
track_steps = unique(max(1,min(parm.ntst,round([.10 .30 .60 1.00]*parm.ntst))));
p_snap = cell(numel(track_steps),1); cb_snap = cell(numel(track_steps),1);
snap_t = zeros(numel(track_steps),1); cb_rms = zeros(parm.ntst,1); p_rms = zeros(parm.ntst,1);
t = zeros(parm.ntst,1); ke = zeros(parm.ntst,1);

fprintf('Running shear-layer case ...\n');
for itst = 1:parm.ntst
    [flow.rhsu,flow.rhsv] = rhs_ns(parm,flow);
    flow = runge_kutta_2d_vec(parm,flow);
    flow = direct_press_corr(parm,flow);
    flow = project(parm,flow);
    
    % --- Filter extraction ---
    [pcb,~] = checkerboard_indicator(parm,flow.p);
    
    t(itst) = itst*parm.dt;
    ke(itst) = .5*sum(flow.u(:).^2+flow.v(:).^2)*parm.dx*parm.dy;
    cb_rms(itst) = rms_field(pcb); p_rms(itst) = rms_field(remove_mean(flow.p));
    
    k = find(track_steps==itst,1);
    if ~isempty(k), p_snap{k}=flow.p; cb_snap{k}=pcb; snap_t(k)=t(itst); end
end

%% Figure 1: physical shear-layer solution
[omega,~]=curl(X,Y,flow.u',flow.v'); speed=sqrt(flow.u'.^2+flow.v'.^2);
f1=figure('Color','w','Position',[100 90 1100 390]);
subplot(1,3,1); plot(t,ke/ke(1),'k-'); grid on; box on; xlabel('Time [s]'); ylabel('E_k/E_k(0)'); title('Kinetic-energy evolution');
subplot(1,3,2); contourf(X,Y,omega,40,'LineStyle','none'); axis equal tight; colorbar; xlabel('x'); ylabel('y'); title('Vorticity');
subplot(1,3,3); contourf(X,Y,speed,40,'LineStyle','none'); hold on; streamslice(X,Y,flow.u',flow.v'); axis equal tight; colorbar; xlabel('x'); ylabel('y'); title('Velocity magnitude and streamlines');
exportgraphics(f1,fullfile(outdir,'Fig1_Shear_Flow.png'),'Resolution',300);

%% Figure 2: raw pressure evolution
clp=max(cellfun(@(q)max(abs(remove_mean(q)),[],'all'),p_snap)); if clp==0,clp=1;end
f2=figure('Color','w','Position',[130 100 1020 680]);
for k=1:numel(track_steps)
    subplot(2,2,k); imagesc(x,y,remove_mean(p_snap{k})'); set(gca,'YDir','normal'); axis equal tight; colormap(gca,parula); caxis([-clp clp]); colorbar;
    title(sprintf('Raw pressure p, t = %.3g s',snap_t(k))); xlabel('x'); ylabel('y');
end
exportgraphics(f2,fullfile(outdir,'Fig2_RawPressure_Evolution.png'),'Resolution',300);

%% Figure 3: grid-scale checkerboard indicator
clcb=max(cellfun(@(q)max(abs(q),[],'all'),cb_snap)); if clcb==0,clcb=1;end
f3=figure('Color','w','Position',[160 120 1020 680]);
for k=1:numel(track_steps)
    subplot(2,2,k); imagesc(x,y,cb_snap{k}'); set(gca,'YDir','normal'); axis equal tight; colormap(gca,bluewhitered(256)); caxis([-clcb clcb]); colorbar;
    title(sprintf('2\\Delta x checkerboard indicator, t = %.3g s',snap_t(k))); xlabel('x'); ylabel('y');
end
exportgraphics(f3,fullfile(outdir,'Fig3_Checkerboard_Indicator.png'),'Resolution',300);

%% Figure 4: cut through the snapshot with the strongest indicator
[~,kpeak]=max(cellfun(@(q)rms_field(q),cb_snap)); jc=round(parm.n/2);
pcut=remove_mean(p_snap{kpeak}); cbcut=cb_snap{kpeak};
f4=figure('Color','w','Position',[200 150 900 400]);
subplot(1,2,1); plot(1:parm.m,pcut(:,jc),'ko-'); grid on; box on; xlabel('Grid index i'); ylabel('p - mean(p)'); title(sprintf('Raw-pressure cut, t = %.3g s',snap_t(kpeak)));
subplot(1,2,2); stem(1:parm.m,cbcut(:,jc),'filled','MarkerSize',4); yline(0,'k:','HandleVisibility','off'); grid on; box on; xlabel('Grid index i'); ylabel('checkerboard indicator'); title('Alternating sign at adjacent nodes');
exportgraphics(f4,fullfile(outdir,'Fig4_2dx_Signature.png'),'Resolution',300);

%% Figure 5: quantitative evolution
f5=figure('Color','w','Position',[230 170 720 400]);
semilogy(t,max(cb_rms,eps),'r-','DisplayName','RMS checkerboard indicator'); hold on; semilogy(t,max(p_rms,eps),'k--','DisplayName','RMS raw pressure');
grid on; box on; xlabel('Time [s]'); ylabel('RMS pressure scale'); title('Evolution of the grid-scale pressure component'); legend('Location','best');
exportgraphics(f5,fullfile(outdir,'Fig5_Checkerboard_Amplitude.png'),'Resolution',300);
fprintf('Figures saved in %s\n',outdir);

% --- Helper Functions ---
function [pcb,psmooth]=checkerboard_indicator(parm,p)
psmooth=zeros(size(p));
for i=1:parm.m
    for j=1:parm.n
        psmooth(i,j)=.25*(p(parm.ip(i),j)+p(parm.im(i),j)+p(i,parm.jp(j))+p(i,parm.jm(j)));
    end
end
pcb=p-psmooth;
end
function q=remove_mean(q), q=q-mean(q(:)); end
function r=rms_field(q), r=sqrt(mean(q(:).^2)); end
function cmap=bluewhitered(n), x=linspace(0,1,n)'; cmap=[min(1,2*x),1-abs(2*x-1),min(1,2-2*x)]; end