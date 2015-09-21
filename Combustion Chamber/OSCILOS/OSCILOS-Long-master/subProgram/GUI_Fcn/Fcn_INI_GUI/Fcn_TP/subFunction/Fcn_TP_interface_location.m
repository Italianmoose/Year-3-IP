function Fcn_TP_interface_location
global CI
% This function is used to locate the required interfaces
%
% first created: 2014-12-03
% last modified: 2014-12-08
% author: Jingxuan Li (jingxuan.li@imperial.ac.uk)
%
CI.indexHA     = find(CI.CD.SectionIndex == 10 | CI.CD.SectionIndex == 11);           % heat addition
CI.indexHP     = find(CI.CD.SectionIndex == 11);                                      % heat perturbations
CI.indexLiner  = find(CI.CD.SectionIndex == 30);                                      % left side of a liner
CI.indexDamper = find(CI.CD.SectionIndex == 2);                                       % damper
%
if ~isempty(CI.indexHA) == 1
    CI.isHA = 1;                                % with heat addition
else
    CI.isHA = 0;                                % without heat addition
end
if ~isempty(CI.indexHP) == 1
    CI.isHP = 1;                                % with heat perturbations
else
    CI.isHP = 0;                                % without heat addition
end
if ~isempty(CI.indexLiner) == 1
    CI.isLiner = 1;                             % with liner
else
    CI.isLiner = 0;                             % without liner
end
if ~isempty(CI.indexDamper) == 1
    CI.isDamper = 1;                            % with dampers
else
    CI.isDamper = 0;                            % without damper
end
%
assignin('base','CI',CI);
%
% -------------------------------end---------------------------------------