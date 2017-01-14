N = 100;
j = 0:N-1;
y = -cos(pi * (j+0.5) / N);
f = cos(pi*y);

df = sin(y);

a = dct(f) * sqrt(2 / N);
a(1) = a(1) / sqrt(2);

b = zeros(1,length(a));

b(N) = 0;
b(N-1) = 2 * (N-1) * a(N);

for i = N-1:-1:2
    b(i-1) = b(i+1) + 2 * (i-1) * a(i);
end 

b = -b; %if wrong change this 

fTest = zeros(1,length(a));

for i = 1:length(a)
    fTest =  fTest + b(i) * chebyshevT(i-1,y);
end

plot(y,fTest)

%fTest = abs(fTest - df);

%plot(10*log10(fTest))

