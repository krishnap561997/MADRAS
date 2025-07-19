function [r, n_so,Qb] = choose_activity(activity,rc)

fpath = '/home/krishnap.kalivel/MATLAB_src/Higher-ACH-Better/Main-Codes/v2/Data';
% fpath = 'C:\Users\nzgheib\University of Florida\COVID-19 - NSF - First Paper\Krishna\Part 5\Evaluate n inf\Nadim Excel Files - n_so';
mouth_breathing = 'Bazant-Bush-mouth breathing (b-nm)'; %1
nose_breathing = 'Bazant-Bush-nose breathing (b-nn)';   %2
sing = 'Bazant-Bush-singing (aah-v-p)';                 %3
soft_sing = 'Bazant-Bush-singing softly (aah-w-p)';     %4
speak = 'Bazant-Bush-speaking (c-v-p)';                 %5
whisper = 'Bazant-Bush-whispering (c-w-p)';             %6

V = 4.5*180;

switch activity
    case 1
        text = mouth_breathing;
        Qb = 0.5/3600;
    case 2
        text = nose_breathing;
        Qb = 0.5/3600;
    case 3
        text = sing;
        Qb = 1/3600;
    case 4
        text = soft_sing;
        Qb = 1/3600;
    case 5
        text = speak;
        Qb = 0.75/3600;
    case 6
        text = whisper;
        Qb = 0.75/3600;
end

fname = strcat([fpath,'/',text]);
% disp(strcat(['Activity: ',text]));
M = readmatrix(fname);
r_digitized = M(:,1);
phi_s_digitized = M(:,2);

[r_digitized,ia,~] = unique(r_digitized);
phi_s_digitized     = phi_s_digitized(ia);

r_digitized     = r_digitized*1e-6;      % Nadim Krishna
phi_s_digitized = phi_s_digitized*1e-16; % Nadim Krishna

r_uniform = linspace(min(r_digitized),max(r_digitized),1000);
phi_s_uniform = interp1(r_digitized,phi_s_digitized,r_uniform);

r       = r_uniform(r_uniform<rc);
phi_s   = phi_s_uniform(r_uniform<rc);
% rc = max(r);        % Nadim
% phi_s = M(:,2);
%n_so is n_so*Vd
n_so = (0.65/3600).*(1+(r./(2.6e-06)).^2+ 0.3/0.65).*V.*phi_s./Qb;         % Nadim

% r = M(:,1);
% rc = max(r);
% phi_s = M(:,2);
% %n_so is n_so*Vd
% n_so = (0.65/3600).*(1+(r./rc).^2+0.3/0.65)...
%     .*V.*phi_s./Qb;
% 
