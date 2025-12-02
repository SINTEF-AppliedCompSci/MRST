function [A, b] = setup_A_b(model)
%
% DESCRIPTION: defines inequalities required to force monotonic behaviour
%              for the point by point history matching, and some saturation
%              function parameterizations like the modified SKJAEVELAND
%
% SYNOPSIS:
%   [A, b] = setup_A_b(model)
%
% PARAMETERS:
%   model - struct with the history matching parameters setup
%
% RETURNS:
%   A, b - inequality matrixes in the format required by matlab optimizer
%   functions like fmincon
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
A = []; b = [];

if model.history_match.kr.status 
    if strcmp(model.history_match.kr.type,'POINT-BY-POINT') 
        if  model.history_match.pc.status
            if strcmp(model.history_match.pc.type,'POINT-BY-POINT')
                Sw_kr_hm = model.history_match.kr.Sw_hm;
                len_sw_kr = length(Sw_kr_hm);
                Sw_pc_hm = model.history_match.pc.Sw_hm;
                len_sw_pc = length(Sw_pc_hm);
                A = zeros(2*len_sw_kr-2+len_sw_pc-1, 2*len_sw_kr-2+len_sw_pc); 
                for i = 1:len_sw_kr-2
                    A(i,i+1) = -1; A(i,i) = 1;
                end
                for i = len_sw_kr:len_sw_kr*2-3
                    A(i,i) = -1; A(i,i+1) = 1;
                end
                for i = 1:len_sw_pc-1
                    A(i+len_sw_kr*2-2,i+len_sw_kr*2-2) = -1; 
                    A(i+len_sw_kr*2-2,i+len_sw_kr*2-1) = 1;
                end
                A(len_sw_kr-1,:) = []; A(2*len_sw_kr-3,:) = [];
                b = zeros(2*len_sw_kr-2+len_sw_pc-3,1);
            end
        else
            Sw_kr_hm = model.history_match.kr.Sw_hm;
            len_sw_kr = length(Sw_kr_hm);
            A = zeros(2*len_sw_kr-3, 2*len_sw_kr-2); 
            for i = 1:len_sw_kr-2
                A(i,i+1) = -1; A(i,i) = 1;
            end
            for i = len_sw_kr:len_sw_kr*2-3
                A(i,i) = -1; A(i,i+1) = 1;
            end
            A(len_sw_kr-1,:) = [];
            b = zeros(2*len_sw_kr-4,1);
        end
    elseif strcmp(model.history_match.kr.type,'MODIEFIED-COREY')
        if  model.history_match.pc.status
            if strcmp(model.history_match.pc.type,'POINT-BY-POINT')
                Sw_pc_hm = model.history_match.pc.Sw_hm;
                len_sw_pc = length(Sw_pc_hm);
                A = zeros(len_sw_pc-1,len_sw_pc+6); 
                for i = 1:len_sw_pc-1
                    A(i,i+6) = -1; A(i,i+7) = 1;
                end
                b = zeros(len_sw_pc-1,1);
            elseif strcmp(model.history_match.pc.type,'MODIFIED-SKJAEVELAND')
                A = zeros(2,12); A(1,11) = -1; A(1,1) = 1; A(2,2) = 1; A(2,12) = 1;
                b = [0.01 ; 0.99];
            end
        end
    elseif strcmp(model.history_match.kr.type,'LET')
        if  model.history_match.pc.status
            if strcmp(model.history_match.pc.type,'POINT-BY-POINT')
                Sw_pc_hm = model.history_match.pc.Sw_hm;
                len_sw_pc = length(Sw_pc_hm);
                A = zeros(len_sw_pc-1,len_sw_pc+10); 
                for i = 1:len_sw_pc-1
                    A(i,i+10) = -1; A(i,i+11) = 1;
                end
                b = zeros(len_sw_pc-1,1);
            elseif strcmp(model.history_match.pc.type,'MODIFIED-SKJAEVELAND')
                A = zeros(2,16); A(1,15) = -1; A(1,1) = 1; A(2,2) = 1; A(2,16) = 1;
                b = [0.01 ; 0.99];
            end
        end
    end
end
if model.history_match.pc.status && not(model.history_match.kr.status)
    if strcmp(model.history_match.pc.type,'POINT-BY-POINT')
        Sw_pc_hm = model.history_match.pc.Sw_hm;
        len_sw_pc = length(Sw_pc_hm);
        A = zeros(len_sw_pc-1,len_sw_pc); 
        for i = 1:len_sw_pc-1
            A(i,i) = -1; A(i,i+1) = 1;
        end
        b = zeros(len_sw_pc-1,1);
    end
end