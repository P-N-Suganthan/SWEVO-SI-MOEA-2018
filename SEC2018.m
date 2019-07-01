%****************************************************************************************************
%  Ten test problems for the special issue in Swarm and Evolutionary Computation on many-objective optimization
%  Author: Dr. Hui Li, Xi'an Jiaotong University, China
%  History of Dates: 
%   ------ 2017.10.02 version 1
%   ------ 2017.11.05 version 2 
%   ------ 2018.5.8   version 3  (the scales of distance function g were changed)
%  If any bug is found in this version, please contact Dr. Hui Li at Email: lihui10@xjtu.edu.cn
%***************************************************************************************************

function [f] = SEC2018(inst, x, m, n)
% inst - name of the instance, e.g., 'MaOP1'-'MaOP10'
% x    - the vector of decision variables
% m, n - the number of objectives, variables
% f    - the vector of objectives
global nvar
global nobj
nvar = n;
nobj = m;
if ~strcmp(inst, 'MaOP1')==0
    f = MaOP1(x);
end
if ~strcmp(inst, 'MaOP2')==0
    f = MaOP2(x);
end
if ~strcmp(inst, 'MaOP3')==0
    f = MaOP3(x);
end
if ~strcmp(inst, 'MaOP4')==0
    f = MaOP4(x);
end
if ~strcmp(inst, 'MaOP5')==0
    f = MaOP5(x);
end

if ~strcmp(inst, 'MaOP6')==0
    f = MaOP6(x);
end

if ~strcmp(inst, 'MaOP7')==0
    f = MaOP7(x);
end

if ~strcmp(inst, 'MaOP8')==0
    f = MaOP8(x);
end

if ~strcmp(inst, 'MaOP9')==0
    f = MaOP9(x);
end

if ~strcmp(inst, 'MaOP10')==0
    f = MaOP10(x);
end

function [f] = MaOP1(x) 
global nvar;
global nobj;
f   = zeros(1, nobj);
g   = sum((x(nobj:nvar) - 0.5).^2 + (1 - cos(20*pi*(x(nobj:nvar)-0.5))))/nvar;  %    '+' before cos is changed to '-' 2018.5.8
tmp = 1;
for m = nobj:(-1):1
    id = nobj - m + 1;                            
    if m>1
        f(m) = (1 + g)*(1 - tmp*(1 - x(id)));      
        tmp  = tmp*x(id);                       
    else
        f(m) = (1 + g)*(1 - tmp);
    end    
    f(m) = (0.1 + 10*m)*f(m);
end


function [f] = MaOP2(x)         
global nvar;
global nobj;
f = zeros(1, nobj);
g = 0;
tmp1 = prod(sin(0.5*pi*x(1:(nobj-1))));
for n=nobj:1:nvar
    if mod(n,5)==0
        g = g + (x(n) - tmp1).^2;
    else
        g = g + (x(n) - 0.5).^2;    
    end;
end
g = g*200;    % The scale 200 is added to g.

tmp2 = 1;
for m = nobj:(-1):1
    p = 2^(mod(m,2)+1);
    if m==nobj
        f(m) = (1 + g)*sin(0.5*x(1)*pi)^p;
    elseif m<nobj&&m>=2        
        tmp2 = tmp2*cos(0.5*pi*x(nobj-m));
        f(m) = (1 + g)*(tmp2*sin(0.5*pi*x(nobj-m+1)))^p;        % wrong pow function
    else
        f(m) = (1 + g)*(tmp2*cos(0.5*pi*x(nobj-1)))^p;          % wrong pow function
    end
end


function [f] = MaOP3(x)
global nvar;
global nobj;
f = zeros(1, nobj);
g = 0;
tmp1 = prod(sin(0.5*pi*x(1:(nobj-1)))); 
for n=nobj:1:nvar
    if mod(n,5)==0
        g = g + n*abs(x(n) - tmp1).^(0.1);  
    else
        g = g + n*abs(x(n) - 0.5).^(0.1);  
    end;
end
tmp2 = 1;
for m = nobj:(-1):1
    if m==nobj
        f(m) = (1 + g)*sin(0.5*x(1)*pi);
    elseif m<nobj&&m>=2   
        tmp2 = tmp2*cos(0.5*pi*x(nobj-m));
        f(m) = (1 + g)*tmp2*sin(0.5*pi*x(nobj-m+1));        
    else
        f(m) = (1 + g)*tmp2*cos(0.5*pi*x(nobj-1));
    end
end

function [f] = MaOP4(x)  %% bias on variable x_1
global nvar;
global nobj;
f = zeros(1, nobj);
g = 0;
tmp1 = prod(sin(0.5*pi*x(1:(nobj-1))));
for n=nobj:1:nvar
    if mod(n,5)==0
        g = g + 2*sin(x(1)*pi)*(abs(-0.9*(x(n) - tmp1).^2) +  abs(x(n) - tmp1).^(0.6));
    else
        g = g + 2*sin(x(1)*pi)*(abs(-0.9*(x(n) - 0.5).^2) +  abs(x(n) - 0.5).^(0.6));
    end;
end
g = g*10;
tmp2 = 1;
for m = nobj:(-1):1
    if m==nobj
        f(m) = (1 + g)*sin(0.5*x(1)*pi);
    elseif m<nobj&&m>=2        
        tmp2 = tmp2*cos(0.5*pi*x(nobj-m));
        f(m) = (1 + g)*tmp2*sin(0.5*pi*x(nobj-m+1));        
    else
        f(m) = (1 + g)*tmp2*cos(0.5*pi*x(nobj-1));
    end
end

function [f] = MaOP5(x)
global nvar;
global nobj;
f = zeros(1, nobj);
g = zeros(1, nobj);
for m=1:nobj
    if m<=3
        g(m) = max(0, -1.4*cos(2*x(1)*pi)) + sum(abs(x(3:nvar) - x(1)*x(2)).^2);
    else
        g(m) = exp((x(m) - x(1)*x(2))^2) - 1;
    end;
    g(m) = 10*g(m);
end
alpha1 = cos(0.5*pi*x(1))*cos(0.5*pi*x(2));
alpha2 = cos(0.5*pi*x(1))*sin(0.5*pi*x(2));
alpha3 = sin(0.5*pi*x(1));
f(1) = (1 + g(1))*alpha1;
f(2) = 4*(1 + g(2))*alpha2;
f(3) = (1 + g(3))*alpha3;
for m = 4:nobj    
    f(m) = (1 + g(m))*(m*alpha1/nobj + (1-m/nobj)*alpha2 + sin(0.5*m*pi/nobj)*alpha3);  % n in g(n) -->g(m)
end

function [f] = MaOP6(x)
global nvar;
global nobj;
f = zeros(1, nobj);
g = zeros(1, nobj);
for m=1:nobj
    if m<=3
        g(m) =  max(0, 1.4*sin(4*x(1)*pi)) + sum(abs(x(3:nvar) - x(1)*x(2)).^2);
    else
        g(m) =  exp((x(m) - x(1)*x(2)).^2) - 1;
    end;
    g(m) = g(m)*10;
end
alpha1 = x(1)*x(2);
alpha2 = (1-x(2))*x(1);
alpha3 = (1-x(1));
f(1) = (1 + g(1))*alpha1;
f(2) = 2*(1 + g(2))*alpha2;
f(3) = 6*(1 + g(3))*alpha3;
for m = 4:nobj
    f(m) = (1 + g(m))*(m*alpha1/nobj + (1-m/nobj)*alpha2 + sin(0.5*m*pi/nobj)*alpha3);  % n in g(n) -->g(m)
end


function [f] = MaOP7(x)
global nvar;
global nobj;
f     = zeros(1, nobj);
alpha = zeros(1, nobj);
g = 0;
tmp = prod(sin(0.5*pi*x(1:(nobj-1))));  
for n=nobj:1:nvar
    if mod(n,5)==0
        g = g + (x(n) - tmp).^2;
    else
        g = g + (x(n) - 0.5).^2;
    end;
end
g = g*100;  % The scale 100 is multiplied by g. 2018.5.8

tau = sqrt(2)/2;
alpha(1) = -(2*x(1)-1).^3 + 1;
T = floor((nobj-1)/2);
for i=1:T
    alpha(2*i)    = x(1) + 2*x(i+1)*tau + tau*abs(2*x(i+1)-1).^(0.5+x(1));
    alpha(2*i+1)  = x(1) - (2*x(i+1)-2)*tau + tau*abs(2*x(i+1)-1).^(0.5+x(1));
end
if mod(nobj,2)==0
    alpha(nobj) = 1 - alpha(1);
end
f = (1 + g)*alpha;


function [f] = MaOP8(x)
global nvar;
global nobj;
f = zeros(1, nobj);
g = 0;
tmp = prod(sin(0.5*pi*x(1:(nobj-1))));  
for n=nobj:1:nvar
    if mod(n,5)==0
        g = g + (x(n) - tmp).^2;
    else
        g = g + (x(n) - 0.5).^2;
    end;
end
g = 100*g;   % The scale 100 is multiplied by g. 2018.5.8


alpha = zeros(1, nobj);
tau = sqrt(2)/2;

alpha(1) = -(2*x(1)-1).^3 + 1;
T = floor((nobj-1)/2);
for i=1:T
    alpha(2*i)    = x(1) + 2*x(i+1)*tau + tau*abs(2*x(i+1)-1).^(1-0.5*sin(4*pi*x(1)));
    alpha(2*i+1)  = x(1) - (2*x(i+1)-2)*tau + tau*abs(2*x(i+1)-1).^(1-0.5*sin(4*pi*x(1)));
end
if mod(nobj,2)==0
    alpha(nobj) = 1 - alpha(1);
end

f = (1 + g)*alpha;


function [f] = MaOP9(x)
global nvar;
global nobj;
f = zeros(1, nobj);
g = 0;
tmp = prod(sin(0.5*pi*x(1:(nobj-1))));  
for n=nobj:1:nvar
    if mod(n,5)==0
        g = g + (x(n) - tmp).^2;
    else
        g = g + (x(n) - 0.5).^2;
    end;
end
g = g*100;    % The scale 100 is multiplied by g. 2018.5.8


alpha = zeros(1, nobj);
tau = sqrt(2)/2;

alpha(1) = -(2*x(1)-1).^3 + 1;
T = floor((nobj-1)/2);
for i=1:T
    z = 2*(2*x(i+1) - floor(2*x(i+1))) - 1;
    alpha(2*i)    = x(1) + 2*x(i+1)*tau + tau*abs(z).^(0.5+x(1));
    alpha(2*i+1)  = x(1) - (2*x(i+1)-2)*tau + tau*abs(z).^(0.5+x(1));
end
if mod(nobj,2)==0
    alpha(nobj) = 1 - alpha(1);
end
f = (1 + g)*alpha;


function [f] = MaOP10(x)
global nvar;
global nobj;
f = zeros(1, nobj);
g = 0;
tmp = prod(sin(0.5*pi*x(1:(nobj-1))));  
for n=nobj:1:nvar
    if mod(n,5)==0
        g = g + (x(n) - tmp).^2;
    else
        g = g + (x(n) - 0.5).^2;
    end;
end

g = g*100;   % The scale 100 is multiplied by g. 2018.5.8

alpha = zeros(1, nobj);
tau = sqrt(2)/2;

alpha(1) = -(2*x(1)-1).^3 + 1;
T = floor((nobj-1)/2);

for i=1:T
    z = 2*(2*x(i+1) - floor(2*x(i+1))) - 1;
    if x(i+1)<0.5
        p = 0.5 + x(1);
    else
        p = 1.5 - x(1);
    end
    alpha(2*i)    = x(1) + 2*x(i+1)*tau + tau*abs(z).^p;
    alpha(2*i+1)  = x(1) - (2*x(i+1)-2)*tau + tau*abs(z).^p;    
end

if mod(nobj,2)==0
    alpha(nobj) = 1 - alpha(1);
end

f = (1 + g)*alpha;
