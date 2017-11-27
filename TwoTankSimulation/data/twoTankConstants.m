%% Initial values

L1i=2;
L2i=3;
T1i=50;
T2i=50;

load('TwoTankOPTrim.mat')
Initials.F1=TwoTankOPTrim.Inputs(1).u;
Initials.F2=TwoTankOPTrim.Inputs(3).u;
Initials.F3=TwoTankOPTrim.Inputs(2).u;
Initials.F4=TwoTankOPTrim.Inputs(4).u;

%% Constants
Area=1;
b=0.5;
cp=1;
kl=0.128;
kv=0.01;
p=1E6;
aHeat=1.41E5;
percN = 0.1;
kv1=Initials.F1/50;
kv2=Initials.F2/50;
kv3=Initials.F3/50;
kv4=Initials.F4/50;

modelParams

%% Inputs
%the function createInputs parse input options for the model into the
%struct inputOpt. Calling without inputs gives default options. See the
%function file for more details on parameter options

inputOpt = createInputs();
%the following parameters have to be re-assigned since Simulink doesn't
%allow variant controls in struct form
T1Input = inputOpt.T1Input;
T2Input = inputOpt.T2Input;
T3Input = inputOpt.T3Input;
T4Input = inputOpt.T4Input;
F1Input = inputOpt.F1Input;
CL = inputOpt.CL;

% E.g. to change the simulation time and introduce a step change of 5
% degrees C in T1in:
% inputOpt = createInputs('SimTime',5000,'T1Input',inputType.STEP,'T1in',[200 25 stepSize 10]);
 
%% Set-points
SPL1=2;
SPL2=3;
SPT1=50;
SPT2=50;

%% Transport Delays
timeDelay.F1 = 1;
timeDelay.F1Out = 0.5;
timeDelay.F2 = 0.5;
timeDelay.F3 = 0.5;
timeDelay.F4 = 0.5;

%% Variable Units
units.F_1 = 'm^3/s';
units.F_2 = 'm^3/s';
units.F_3 = 'm^3/s';
units.F_4 = 'm^3/s';
units.T_1_i_n = '\circC';
units.T_1 = '\circC';
units.T_2_i_n = '\circC';
units.T_2 = '\circC';
units.T_3 = '\circC';
units.T_4 = '\circC';
units.L_1 = 'm';
units.L_2 = 'm';


