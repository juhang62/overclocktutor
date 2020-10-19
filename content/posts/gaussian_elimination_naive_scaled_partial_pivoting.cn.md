---
title: "Matlab code for Gaussian elimination (naïve, partial pivoting, scaled partial pivoting)"
date: 2020-10-18T22:11:59-07:00
draft: false
---

Gaussian elimination is an algorithm to solve linear systems. Its naive version is usually taught as early as in your algebra class. When it is revisited later in a numerical analysis course, pivoting strategy is introduced to avoid amplified round-off error that plagues the Naive approach.

I frequently receive questions from students I tutor about how to code different versions of Gaussian elimination in Matlab. So here you go: all three versions of Gaussian elimination in one function:

```matlab
function x=ge(A,b,options)
% solve a linear system using Gaussian Elimination
% specify options to be one of the following: 
% ge: naive Gaussian Elimination
% pp: Gaussian Elimination with partial pivoting 
% spp: Gaussian Elimination with scaled partial pivoting 

% Step 1: Gaussian Elimination with scaled partial pivoting. 
%         Overwrite the matrix.
n=length(b);
p = (1:n)';	        % initialize the pivoting vector
if strcmp(options,'pp')
    s=ones(1,n);
elseif strcmp(options,'spp')
    s = max(abs(A'));% compute the scale of each row
end
for k = 1:(n-1)
    if strcmp(options,'pp') || strcmp(options,'spp')
        r = abs(A(p(k),k)/s(p(k)));
        kp = k;
        for i = (k+1):n
            t = abs(A(p(i),k)/s(p(i)));
            if t > r,  r = t;  kp = i;  end
        end
        l = p(kp);  p(kp) = p(k);  p(k) = l;   % interchange p(kp) and p(k)
    end
    for i = (k+1):n
        A(p(i),k) = A(p(i),k)/A(p(k),k);
        for j = (k+1):n
            A(p(i),j) = A(p(i),j)-A(p(i),k)*A(p(k),j);
        end
    end
end
%p                      % output the pivoting vector
%A                      % output the overwritten matrix

% Step 2: Forward substitution to solve  L*y = b, where
%         L(i,j) = 0 for j>i; L(i,i) = 1;  L(i,j) = A(p(i),j) for i>j.

y = zeros(n,1);        % initialize y to be a column vector
y(1) = b(p(1));
for i = 2:n
    y(i) = b(p(i));
    for j = 1:(i-1)
        y(i) = y(i)-A(p(i),j)*y(j);
    end
end
%y                      % output y

% Step 3: Back substitution to solve  U*x = y
%         U(i,j) = A(p(i),j) for j>=i; U(i,j) = 0 for i>j.

x = zeros(n,1);        % initialize x to be a column vector
x(n) = y(n)/A(p(n),n);
for i = (n-1): -1 : 1
    x(i) = y(i);
    for j = (i+1):n
        x(i) = x(i)-A(p(i),j)*x(j);
    end
    x(i) = x(i)/A(p(i),i);
end
end
```      
I have to admit that I did not write this code from scratch. I probably did it at some point of my student life and am not downplaying the importance to write the code yourself as a way to enhance learning. Nevertheless, be aware that there must be some code shared online for well-known algorithms (e.g. Gaussian elimination).  A simple Google search "scaled partial pivoting matlab" landed me to [this](http://www.math.ucsd.edu/~bli/teaching/math170Af06/gauss_sp.html).  After verifying it is a valid implementation of Gaussian elimination with scaled partial pivoting, I knew I just needed a few modifications to get the other two versions of Gaussian elimination. Pay attention to the if statements in the code that checks `option`.  That is, if we set scale factors to be ones, we get partial pivoting and if we totally skipped exchanging rows, we get naive Gaussian elimination. It is just that simple if you understand the differences between these three versions. 

That is for this post. If you have more questions, contact me. I am always here to help. 
