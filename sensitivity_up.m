close all;
clear all;
clc;
ys = zeros(1,46);
for i = 1:46
    ys(i) = 1;
end
yi = [100 80 70 60 50 0 0 40 0 30 0 20 0 10 0 50 0 0 20 0 40 0 10 0];
y0 = [yi,ys];
tspan = [0 40];
[t,y] = ode45(@sen_up,tspan,y0);
D = {'V', 'R', 'VR', 'VR2', 'NCK', 'NCKA', 'NCKA_VR2', 'PAK2', 'PAK2A', 'P38',...
    'P38A', 'MAPKAPK2', 'MAPKAPK2A', 'HSP27', 'HSP27A', 'SHB', 'SHBA', 'SHBA_VR2',...
    'PRAK', 'PRAKA', 'FAK', 'FAKA', 'PAXILLIN', 'PAXILLINA', 'dSh1/dt', 'dSh2/dt', ...
    'dSh3/dt', 'dSh4/dt', 'dSh5/dt', 'dSh6/dt', 'dSh7/dt', 'dSh8/dt', 'dSh9/dt', 'dSh10/dt', ...
    'dSh11/dt', 'dSh12/dt', 'dSh13/dt', 'dSh14/dt', 'dSh15/dt', 'dSh16/dt', 'dSh17/dt', ...
    'dSh18/dt', 'dSh19/dt', 'dSh20/dt', 'dSh21/dt','dSh22/dt','dSh23/dt'...
     'dSp1/dt', 'dSp2/dt', 'dSp3/dt', 'dSp4/dt', 'dSp5/dt', 'dSp6/dt', 'dSp7/dt', ...
     'dSp8/dt', 'dSp9/dt', 'dSp10/dt', 'dSp11/dt', 'dSp12/dt', 'dSp13/dt', 'dSp14/dt', ...
     'dSp15/dt', 'dSp16/dt', 'dSp17/dt', 'dSp18/dt', 'dSp19/dt', 'dSp20/dt', 'dSp21/dt','dSp22/dt','dSp23/dt'};
 
s = {'Sh1', 'Sh2', 'Sh3', 'Sh4', 'Sh5', 'Sh6', 'Sh7',...
    'Sh8', 'Sh9', 'Sh10', 'Sh11', 'Sh12', 'Sh13', 'Sh14',...
    'Sh15', 'Sh16', 'Sh17', 'Sh18', 'Sh19', 'Sh20', 'Sh21', 'Sh22', 'Sh23'...
    'Sp1', 'Sp2', 'Sp3', 'Sp4', 'Sp5', 'Sp6', 'Sp7', 'Sp8', 'Sp9', 'Sp10',...
    'Sp11', 'Sp12', 'Sp13', 'Sp14', 'Sp15', 'Sp16', 'Sp17', 'Sp18', 'Sp19',...
    'Sp20', 'Sp21', 'Sp22', 'Sp23'};

for i = [25,36,37,47,48,68,69,70]
    figure(i)   
    plot(t,y(:,i),'linewidth', 2)
    ylabel(s(i-24))
    xlabel('Time')
    str = D(i);
    title(str)
    file = s(i-24)+".jpg";
    saveas(figure(i),file ,'jpg');
end 

disp('Variation of Sensitivity of the final protein HSP27* w.r.t all rate constants')
vec = y(:,25:47);
A = vec'*vec;
[e,v] = eig(A);
disp("Values of the Eigen vector of the matrix formed by using sensitivity values wrt HSP27*")
disp(diag(v))
[val,ind] = max(diag(v));
disp("Values of the Eigen vector with Maximum value among vectors ")
% disp(e(:,ind))
disp('    Coefficient       Value')
disp([D(25:47)', cellstr(num2str(e(:, ind)))]);
% -------------------------------------------------------------------------------
% -------------------------------------------------------------------------------
disp('Variation of Sensitivity of the final protein Paxillin* w.r.t all rate constants')
vecp = y(:,48:70);
Ap = vecp'*vecp;
[ep,vp] = eig(Ap);
disp("Values of the Eigen vector of the matrix formed by using sensitivity values wrt Paxillin*")
disp(diag(vp))
[valp,indp] = max(diag(vp));
disp("Values of the Eigen vector with Maximum value among vectors ")
% disp(ep(:,indp))
disp('    Coefficient       Value')
disp([D(48:70)', cellstr(num2str(ep(:, indp)))]);