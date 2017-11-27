function [ dataObj ] = parsLogSigs( logsout,units,source,faultSource,faultType,faultSize )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

signalNames=logsout.getElementNames;%logged variable names
Time=logsout.get(signalNames{1}).Values.Time;%time values (taken from first logged signal)

N=length(Time);%number of samples
M=logsout.numElements;%number of variables

%combined data matrix of all signals
dataMatrix=zeros(N,M);%initialisation
% figure
for m=1:M
    dataMatrix(:,m)=logsout.get(signalNames{m}).Values.Data;
%     subplot(2,1,m)
%     plot(dataMatrix(:,m))
%     ylabel(signalNames(m))
end
    
dataObj = FDDataObject(dataMatrix,units,Time,'min',signalNames,source,faultSource,faultType,faultSize);

end

