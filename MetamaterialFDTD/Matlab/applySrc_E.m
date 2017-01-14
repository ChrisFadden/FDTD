function [field] = applySrc_E(field,grid,src,constants,t)

%Hard/Soft Source Implementation
%field.Ez(grid.i_src) = sin(2*pi*src.freq * grid.dt * t);
%field.Ez(grid.i_src) = field.Ez(grid.i_src) + 2*exp(-((t*grid.dt - src.t0)/src.tw)^2); 

%TF/SF Source Implementation
offset = (grid.TFSF_end - grid.TFSF_start)*grid.dt;

%Gaussian Source
srcStart = exp(-((t*grid.dt - src.t0)/src.tw)^2);
srcEnd = exp(-((t*grid.dt - src.t0 - offset)/src.tw)^2);

%Harmonic Source (need no u(t - tfsf) since the "error" is just in
%                 scattered field, it won't reflect off a structure
%srcStart = sin(2*pi*t*grid.dt*src.freq);
%srcEnd = sin(2*pi*src.freq*(t*grid.dt - offset));

%Field Updates
field.Ez(grid.TFSF_start) = field.Ez(grid.TFSF_start) + srcStart;
field.Ez(grid.TFSF_end) = field.Ez(grid.TFSF_end) - ...
                          srcEnd;  
                           

% field.Ez(grid.TFSF_start) = field.Ez(grid.TFSF_start) + ...
%                             sin(2*pi*t*grid.dt*src.freq);
%  
%  offset = (grid.TFSF_end - grid.TFSF_start)*grid.dt;
%  
%  field.Ez(grid.TFSF_end) = field.Ez(grid.TFSF_end) - ...
%                              sin(2*pi*src.freq*(t*grid.dt - offset));


end

