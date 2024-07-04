function [out,inds] = subamostra(x,nsamples)

s = sort(x);
if isinf(s(1))
    mx=s(2);
else
    mx=s(1);
end
subs = max(x)-mx;

frac = subs/nsamples;

intervalo = mx:frac:max(x);
out = zeros(size(intervalo));
inds= zeros(size(intervalo));

for i = 1:length(intervalo)
    test = abs(x - intervalo(i));
    ind = find(test==min(test));
    if length(ind) < 1
        disp(" ----------------  deu ruim ------------------");
    end
    out(i) = x(ind(1));
    inds(i)=ind(1);
end