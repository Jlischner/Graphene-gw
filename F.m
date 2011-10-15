function out = F(s1,s2,k1,k2);

absk1 = sqrt(sum(k1.^2,2));
absk2 = sqrt(sum(k2.^2,2));
cosk1k2 = transpose(k1*k2')./absk2./absk1;
out = (1+s1*s2*cosk1k2)/2;

endfunction;