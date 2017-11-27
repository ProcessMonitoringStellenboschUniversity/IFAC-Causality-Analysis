% Model Parameters

%% Define Variant Controls
T1Oscillating = Simulink.Variant('T1Input == inputType.OSCILLATING');
T1Step = Simulink.Variant('T1Input == inputType.STEP');
T1Ramp = Simulink.Variant('T1Input == inputType.RAMP');
T1RW = Simulink.Variant('T1Input == inputType.RW');

T2Oscillating = Simulink.Variant('T2Input == inputType.OSCILLATING');
T2Step = Simulink.Variant('T2Input == inputType.STEP');

T3Oscillating = Simulink.Variant('T3Input == inputType.OSCILLATING');
T3Step = Simulink.Variant('T3Input == inputType.STEP');

T4Oscillating = Simulink.Variant('T4Input == inputType.OSCILLATING');
T4Step = Simulink.Variant('T4Input == inputType.STEP');

F1Oscillating = Simulink.Variant('F1Input == inputType.OSCILLATING');
F1Step = Simulink.Variant('F1Input == inputType.STEP');

closedLoop = Simulink.Variant('CL == true');
openLoop = Simulink.Variant('CL == false');
