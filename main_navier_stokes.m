%**************************************************************************
% CFD Lab Summer Semester 2018, Assignment 8
%
% This code solves the 2D unsteady Navier-Stokes equation 
%
% For spatial discretisation a central difference scheme is used
% and for time integration a Runge-Kutta method is implemented.
%
% Two cases "Channel" and "Shear" are defined.
%
% authors: D.Quosdorf & Y.Sakai
% june, 2018
%**************************************************************************


% close figures, command window and clear memory
close all; clc; clear

% read infile 
infilename = 'infile_shear.mat';
fprintf('infilename is: %s\n', infilename)

% build structures 'parm' and 'flow'
[parm, flow] = build_structs;
fprintf('struct built\n')

% fill some fields of 'parm' and 'flow' with data from infile
[parm, flow] = set_params(parm, flow, infilename);
fprintf('parameters set\n')

% initialisation of flow field
[parm, flow] = initialize(parm, flow);
fprintf('flow field initialised\n')

% start time integration
[parm, flow] = timeint(parm, flow);
