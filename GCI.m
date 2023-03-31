%% Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications
%% Fluids Engineering Division Of ASME (Celik et al.,2008)

clc
clear
close all

%% Initialization:

V = 1.0; % Total volume of the domain/mesh

numItr = 70; % Number of solver iterations
conv_crt = 0.0001; % Convergence criterion

nbr_elements = [18000; 8000; 4500];
phi = [6.063; 5.972; 5.863];

% If the values are in a CVS file

% datatable = readtable('data.csv');
% nbr_elements = datatable.nbr_elements;
% phi = datatable.phi;

table = zeros(17,1);

%% Solving:

[gci_21_p,gci_32_p,p,r21_p,r32_p,e_21ext_p,e_32ext_p,phi_21ext_p,phi_32ext_p,ea_21_p,ea_32_p] = gci(numItr,conv_crt,V,nbr_elements(1),phi(1),nbr_elements(2),phi(2),nbr_elements(3),phi(3));

%% Plot:

table(1) = nbr_elements(1);
table(2) = nbr_elements(2);
table(3) = nbr_elements(3);
table(4) = r21_p;
table(5) = r32_p;
table(6) = phi(1);
table(7) = phi(2);
table(8) = phi(3);
table(9) = p;
table(10) = phi_21ext_p;
table(11) = phi_32ext_p;
table(12) = ea_21_p;
table(13) = ea_32_p;
table(14) = e_21ext_p;
table(15) = e_32ext_p;
table(16) = gci_21_p;
table(17) = gci_32_p;

%% Functions:

function chart_error(vect)
    res = zeros(length(vect)-1,1);
    for i=2:length(vect)
        res(i-1) = abs((vect(i-1)-vect(i))/vect(i-1));
    end
    figure
    plot(res,'k-x','LineWidth',1.75)
    grid on
end 

function [gci_21,gci_32,p,r21,r32,e_21ext,e_32ext,phi_21ext,phi_32ext,ea_21,ea_32] = gci(numItr, conv_crt, V, N_1, phi_1, N_2, phi_2, N_3, phi_3)
    %% Initialisatiton:

    % 2D:
    h1=(V/N_1)^(1/2);
    h2=(V/N_2)^(1/2);
    h3=(V/N_3)^(1/2);

    % 3D:
%     h1=(V/N_1)^(1/3);
%     h2=(V/N_2)^(1/3);
%     h3=(V/N_3)^(1/3);
    
    r21=h2/h1;
    r32=h3/h2;
    
    epsiloN_32 = phi_3-phi_2;
    epsiloN_21 = phi_2-phi_1;
    
    s = sign(epsiloN_32/epsiloN_21);
    Q_func = @(p) log((r21^(p)-s)/(r32^(p)-s));
    P_func = @(q) (1/log(r21))*abs( log( abs(epsiloN_32/epsiloN_21) ) + q);
    
    %% Solving:
    
    p = fixed_point(P_func,Q_func,0,conv_crt,numItr);
    
    %% Results:
    
    phi_21ext = extrapolated_values(p,r21,phi_2,phi_1);
    phi_32ext = extrapolated_values(p,r32,phi_3,phi_2);
    
    ea_21 = approximate_relative_error(phi_2,phi_1);
    ea_32 = approximate_relative_error(phi_3,phi_2);
    
    e_21ext = extrapolated_relative_error(phi_21ext,phi_1);
    e_32ext = extrapolated_relative_error(phi_32ext,phi_2);
    
    gci_21 = gci_value(p,ea_21,r21);
    gci_32 = gci_value(p,ea_32,r32);
end

function res = gci_value(p,ea,r)
    res = (1.25*ea)/(r^p-1);
end

function e = extrapolated_relative_error(phi_ext,phi)
    e = abs((phi_ext-phi)/phi_ext);
end

function e = approximate_relative_error(phi_2,phi_1)
    e = abs((phi_1-phi_2)/phi_1);
end

function phi = extrapolated_values(p,r,phi_2,phi_1)
    phi = ((r^p)*phi_1-phi_2)/((r)^p-1);
end

function sol=fixed_point(P_func,Q_func,q_init,conv_crt,numItr)
    i=1;
    p1=feval(P_func,q_init);
    q1=feval(Q_func,p1);

    p2=feval(P_func,q1);
    q2=feval(Q_func,p2);

    if q2==q1 && p2==p1
        fprintf('The fixed point is %f', y)
    end

    while abs(q2-q1)>conv_crt && abs(p2-p1)>conv_crt && i+1<=numItr % Continue until the residual is smaler than the convergence criteria or if the number of iteration is exceeded
        
        i=i+1;
        
        p1=p2;
        q1=q2;

        p2=feval(P_func,q1);
        q2=feval(Q_func,p2);

    end

    if i == numItr
        fprintf('Algorithm stop! Too many iteration: %f', i)
    end

    sol = p2;

end


