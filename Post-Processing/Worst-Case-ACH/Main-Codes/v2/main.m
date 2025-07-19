% Comments: 
% Use eta_ef > 0.33 for now. I will fix eta_ef = 0 and let you know.
clear; clc; close all;
%% Input Parameters
% Enter manually

ACH_comp_list = [0.5 2.5 5 10];           % Air Changes per Hour 
H = 3.2; V = 320;           % Height and Volume of the room (m, m^3)
d = 8;                      % Source to Sink separation (m)
eta_ef = 1;                 % Filter efficiency
lam_v = 0.3/3600;           % Viral deactivation rate (From Bazant & Bush)
activity_list = {'Mouth_Breathing','Nose_Breathing','Singing','Soft_Singing',...
    'Speaking','Whispering'};
activity = 1;
dT = 300;        %Time between first appearance of source and sink.
T_si = 5*60;
T_so = 5*60 + dT;

% Prompts
% T_si = 60*input('Enter T_si (min): ');  %Time of occupation of Sink
% T_so = 60*input('Enter T_so (min): ');  %Time of occupation of Source
% d = input('Enter Source to Sink distance (m) (10 m > d > 0.4 m) : ');   %Source Sink separation
% H = input('Enter height of the room (m): ');
% V_room = input('Enter volume of the room (m): ');
    
for i=1:numel(ACH_comp_list)
    ACH = ACH_comp_list(i); ACS = ACH/3600;
%% Concentration from Well-Mixed model (n_wm)
    n_wm = get_nwm(ACH, H,V, eta_ef, lam_v, activity, max(T_si+dT,T_so));

%% Concentration for specified source-sink separation
    n_d = get_nd(n_wm, ACH, eta_ef,d,dT);
    
%     plot(n_d,'LineWidth',3); hold on;

%% Cumulative Exposure at Sink
    D_d(i) = get_exposure(n_d,ACS,V,T_si,T_so,dT);
end
legend({'ACH 0.5','ACH 2.5','ACH 5','ACH 10'});