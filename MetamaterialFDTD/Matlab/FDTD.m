%MATLAB FDTD Implementation
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%
%   Constants
%%%%%%%%%%%%%%%%%%%%%%%%
constants.cc    = 299792458;
constants.mu0   = 4e-7 * pi;
constants.eps0  = 1.0 / (constants.cc^2 * constants.mu0);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Source
%%%%%%%%%%%%%%%%%%%%%%%%%
source.freq = 1e9;
source.tw   = 0.5 / source.freq;
source.t0   = 4.0 * source.tw;


%%%%%%%%%%%%%%%%%%%%%%%%%
%   Grid
%%%%%%%%%%%%%%%%%%%%%%%%%
grid.CFL     = 1;
grid.NLambda = 20; 
grid.dx      = constants.cc / (source.freq * grid.NLambda);
grid.dt      = grid.CFL * grid.dx / constants.cc;

grid.sizeX          = round(grid.NLambda * 10);
grid.i_src          = grid.sizeX / 2;
grid.TFSF_start     = grid.NLambda;
grid.TFSF_end       = 5*grid.NLambda;
grid.maxTime        = 100;

%%%%%%%%%%%%%
%   Material
%%%%%%%%%%%%%
material(1).start = round(2*grid.NLambda);
material(1).end   = round(3*grid.NLambda);
material(1).epsr  = 1;  %Relative Permittivity
material(1).mur   = 1;  %Relative Permeability
material(1).electricLoss  = 0;  %Electric Loss [sig * dt /(2*eps)]
material(1).magneticLoss  = 0;  %Magnetic Loss [sigM * dt / (2*mu)]

material(2).start = round(3*grid.NLambda)+1;
material(2).end   = round(4*grid.NLambda);
material(2).epsr  = 1;  %Relative Permittivity
material(2).mur   = 1;  %Relative Permeability
material(2).electricLoss  = 0;  %Electric Loss [sig * dt /(2*eps)]
material(2).magneticLoss  = 0;  %Magnetic Loss [sigM * dt / (2*mu)]

%%%%%%%%%%%
%   Fields
%%%%%%%%%%%
field.Ez    = zeros(1,grid.sizeX);
field.Hy    = zeros(1,grid.sizeX);
field.Ceze  = 1;
field.Cezh  = grid.dt / (constants.eps0 * grid.dx);
field.Chyh  = 1;
field.Chye  = grid.dt / (constants.mu0 * grid.dx);

for t = 1:grid.maxTime
    
    %Magnetic Field Updates
    [field]    = applyBC_H(field,grid);
    [field]    = updateH(field,grid,material);
    [field]    = applySrc_H(field,grid,source,constants,t);
    
    %Electric Field Updates
    [field]    = applySrc_E(field,grid,source,constants,t);
    [field]    = applyBC_E(field,grid);
    [field]    = updateE(field,grid,material);
   
end

plot(field.Ez)