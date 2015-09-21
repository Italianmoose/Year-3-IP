function [Lf,tauf] = Fcn_flame_model_NL_JLi_AMorgans(uRatio,indexHP)
global CI
% This function is used to calculate the nonlinear gain ratio and time
% delay of flame describing function of our model
% indexHP is the index of heat perturbation interface
%
% first created: 2014-11-11
% last edit: 2014-12-15
%
%
HP      = CI.FM.HP{indexHP};  
Lf      = interp1(  CI.FM.HP{1,indexHP}.NL.Model3.uRatio,...
                    CI.FM.HP{1,indexHP}.NL.Model3.Lf,...
                    uRatio,'linear','extrap');
tauf    = HP.FTF.tauf + HP.NL.Model3.taufN.*(1-Lf);
%
% ---------------------------end-------------------------------------------
