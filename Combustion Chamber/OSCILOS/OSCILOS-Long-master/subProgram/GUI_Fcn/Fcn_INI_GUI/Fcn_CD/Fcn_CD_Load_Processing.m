function Fcn_CD_Load_Processing
% This function is used to process the loaded data, 
% The objective is to:
% split the original shape, based on the TubeIndex
%
% first created: 2014-12-05
% last edited: 2014-12-05
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
%
global CI
switch CI.CD.isNoGC
    case 1      % No gradually area change section
        %
    case 0      % with gradually area change section
        GUI_INI_CD_Split
        NumSplit = CI.CD.GC.numSplitting;
        for ss = 1:length(CI.CD.indexGC)
            k = CI.CD.indexGC(ss);
            N = length(CI.CD.x_sample); % get the current length of x_sample
            CI.CD.x_sample((k+1:N) + NumSplit(ss)-2)        = CI.CD.x_sample(k+1:N); % move forward by NumSplit(ss)-1
            CI.CD.r_sample((k+1:N) + NumSplit(ss)-2)        = CI.CD.r_sample(k+1:N); % move forward by NumSplit(ss)-1
            CI.CD.SectionIndex((k+1:N) + NumSplit(ss)-2)    = CI.CD.SectionIndex(k+1:N);
            CI.CD.TubeIndex((k+1:N) + NumSplit(ss)-2)       = CI.CD.TubeIndex(k+1:N);
            %
            CI.CD.x_sample(k:k+NumSplit(ss)-1) = linspace(CI.CD.x_sample(k),CI.CD.x_sample(k+NumSplit(ss)-1), NumSplit(ss));
            CI.CD.r_sample(k:k+NumSplit(ss)-1) = linspace(CI.CD.r_sample(k),CI.CD.r_sample(k+NumSplit(ss)-1), NumSplit(ss));
            CI.CD.SectionIndex(k:k+NumSplit(ss)-2) = 0;
            CI.CD.TubeIndex(k:k+NumSplit(ss)-2) = 1;
            %
            if ss < length(CI.CD.indexGC)
                CI.CD.indexGC(ss+1:end) = CI.CD.indexGC(ss+1:end) + NumSplit(ss)-2;
            end            
        end
end
% %
assignin('base','CI',CI)
%
% -------------------------------end --------------------------------------