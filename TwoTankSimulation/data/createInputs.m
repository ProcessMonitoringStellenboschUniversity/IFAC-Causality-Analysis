function inputOpt = createInputs(varargin)
%Function that parses input options
% INPUT: paramater value pairs 

p = inputParser;
addParameter(p,'CL',false); %open or closed loop operation
addParameter(p,'T1in',[200 25 5 10]);
% inputType variant controls. Either CONSTANT, STEP, or OSCILLATING
addParameter(p,'T1Input',inputType.CONSTANT); 
addParameter(p,'T2in',[200 25 5 10]);
addParameter(p,'T2Input',inputType.CONSTANT);
addParameter(p,'T3',[200 100 5 10]);
addParameter(p,'T3Input',inputType.CONSTANT);
addParameter(p,'T4',[200 100 5 10]);
addParameter(p,'T4Input',inputType.CONSTANT);
addParameter(p,'F1',[200 0.1810 0.2 10]);
addParameter(p,'F1Input',inputType.CONSTANT);
addParameter(p,'addNoise',true);
addParameter(p,'simTime',500);
addParameter(p,'sampleTime',0.2);
addParameter(p,'piWriteEvent',false);
parse(p,varargin{:});
inputOpt = p.Results;

% T1in Input
inputOpt.T1in = [];
inputOpt.T1in.time = p.Results.T1in(1);% start time
inputOpt.T1in.start = p.Results.T1in(2);% initial value
inputOpt.T1in.end = p.Results.T1in(2)+ p.Results.T1in(3);%final value for step change (initial plus third entry)
inputOpt.T1in.amplitude = p.Results.T1in(3);%Oscillation amplitude
inputOpt.T1in.period = p.Results.T1in(4);%Oscillation period
% T2in Input
inputOpt.T2in = [];
inputOpt.T2in.time = p.Results.T2in(1);
inputOpt.T2in.start = p.Results.T2in(2);
inputOpt.T2in.end = p.Results.T2in(2) + p.Results.T2in(3);
inputOpt.T2in.amplitude = p.Results.T2in(3);
inputOpt.T2in.period = p.Results.T2in(4);
% T3 Input
inputOpt.T3= [];
inputOpt.T3.time = p.Results.T3(1);
inputOpt.T3.start = p.Results.T3(2);
inputOpt.T3.end = p.Results.T3(2) + p.Results.T3(3);
inputOpt.T3.amplitude = p.Results.T3(3);
inputOpt.T3.period = p.Results.T3(4);
% T4 Input
inputOpt.T4 = [];
inputOpt.T4.time = p.Results.T4(1);
inputOpt.T4.start = p.Results.T4(2);
inputOpt.T4.end = p.Results.T4(2) + p.Results.T4(3);
inputOpt.T4.amplitude = p.Results.T4(3);
inputOpt.T4.period = p.Results.T4(4);
% F1 Input
inputOpt.F1 = [];
inputOpt.F1.time = p.Results.F1(1);
inputOpt.F1.start = p.Results.F1(2);
inputOpt.F1.end = p.Results.F1(2) + p.Results.F1(3);
inputOpt.F1.amplitude = p.Results.F1(3);
inputOpt.F1.period = p.Results.F1(4);
end % function



