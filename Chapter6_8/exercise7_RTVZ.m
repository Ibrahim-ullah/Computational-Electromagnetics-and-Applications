% Matlab code for 1D electromagnetic FDTD program


L = 5;  %domain length in meters
N = 505;    %spatial samples in domain
Niter = 800;    %number of iterations to perform
fs = 300e6;     %source frequency in Hz

ds = L/N;   %spatial steps in meters
dt = ds/300e6;  %magic time in step

eps0 = 8.854e-12;  %permittivity of free space
mu0 = pi*4e-7;      %permeability of free space
x = linspace(0,L,N);    %coordinate of spatial samples
showWKB=0;          %if=1 then show WKB approx at the end

%scale factors for E and H

ae = ones(N,1)*dt/(ds*eps0);
am = ones(N,1)*dt/(ds*mu0);
as = ones(N,1);
epsr = ones(N,1);
mur = ones(N,1);
sigma = zeros(N,1);



%Epsilon, sigma and mu profiles are specified. Try profile = 1,2,3,4,5,6
%in sequence. You can define epsr(i), mur(i), (relative permittivity and
%permeability and sigma(i) conductivity)

profile =2;

for i=1:N
    epsr(i) = 1;
    mur(i)  =1;
    w1 = 0.5;
    w2 = 1.5;
    if (profile==1)     %dielectric window
        if (abs(x(i)-L/2)<0.5) 
            epsr(i)=4;
        end
    end
    
    if (profile==2)     %dielectric window with smooth trase
       if (abs(x(i)-L/2)<1.5) 
            epsr(i)=1+3*(1+cos(pi*(abs(x(i)-L/2))));
       end
            if (abs(x(i)-L/2)<0.5) 
                epsr(i)=4;
            end
    end
       
    if (profile ==3)   %dielectric discontinuity
        if(x(i) > L/2) 
            epsr(i) = 9;
        end
    end
    
    if (profile==4)     %dielectric discontinuity with 1/4 wave matching layer
        if (x(i)>(L/2-0.1443))
            epsr(i)=3;
        end
        if (x(i)>L/2)
            epsr(i)=9;
        end
    end
    
    
    if (profile ==5)   %conducting half space
        if(x(i)>L/2)
            sigma(i) = 0.005;
        end
    end
    
    if (profile ==6)  %sinusoidal dielectric
        epsr(i) = 1+sin(2*pi*x(i)/L)^2;
    end
    
end

ae = ae./epsr;
am = am./mur;
ae = ae./(1+dt*(sigma./epsr)/(2*eps0));
as = (1-dt*(sigma./epsr)/(2*eps0))./(1*dt*(sigma./epsr)/(2*eps0));
    
    
    
%plot the permittivity , permeability & conductivity profiles

figure(1)
subplot(3,1,1);
plot(x,epsr);
grid on;
axis([3*ds L min(epsr)*0.9 max(epsr)*1.1]);
title('relative permittivity');
subplot(3,1,2);
plot(x,mur);
grid on;
axis([3*ds L min(mur)*0.9 max(mur)*1.1]);
title('relative permeability');
subplot(3,1,3);
plot(x,sigma);
grid on;
axis([3*ds L min(sigma)*0.9-0.001 max(sigma)*1.1+0.001]);
title('conductivity');


% The code you provided above (up to the definition of 'ae', 'am', 'as', etc.)

% Initialize fields
E = zeros(N, 1);
H = zeros(N, 1);

% Initialize plot
figure;
h = plot(x, E);
title('Traveling Wave');
xlabel('Position (m)');
ylabel('E-field');
xlim([0, L]);
ylim([-1, 1]);
drawnow;

% Main time-stepping loop
for t = 1:Niter
    % Update electric field
    E(2:N) = as(2:N) .* E(2:N) + ae(2:N) .* (H(2:N) - H(1:N-1));

    % Update magnetic field
    H(1:N-1) = H(1:N-1) + am(1:N-1) .* (E(2:N) - E(1:N-1));

    % Source excitation (Gaussian pulse)
    source = exp(-0.5 * ((t - 30) / 10)^2) * sin(2 * pi * fs * t * dt);
    E(1) = source;

    % Update plot
    set(h, 'YData', E);
    title(['Traveling Wave at Time Step ' num2str(t)]);
    drawnow;
end





