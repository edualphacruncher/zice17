% This script tests some matrix-vector things
%
% Stefan Wild. Argonne National Laboratory, January 2016.

% Set this to 1 or 2:
examplenum = 2;

rand('state',1);

if examplenum==1 % test backslash versus inverse matrix form
    for n=2.^[8:2:12];
        A = rand(n);
        b = rand(n,1);
        c = zeros(n,1);
        
        n
         
        tic;
        c = inv(A)*b;
        fullinv = toc
        
        tic
        c = A\b;
        backslash=toc
    end
elseif examplenum==2 % test Todd's dense example (I + u*v')*x
    for n=2.^[12:2:16];
        u = rand(n,1);
        v = rand(n,1);
        x = rand(n,1);
        
        c = zeros(n,1); % initialize
        
        n
        
        tic
        c = x+u*(v'*x);
        matvec=toc
        
        
        tic;
        c = (speye(n)+u*v')*x;
        form_and_mult = toc
    end
end

