function [h,ho]  = GCEHistogram_modification(f,b,w,alpha,g) 
%GCE Histogram Modification Algorithm
%
% Paper:
% T.  Arici,  S.  Dikbas  and  Y.  Altunbasak,  “A  histogram 
% modification  framework  and  its  application  for  image 
% contrast  enhancement”,  IEEE  Transactions  on  image 
% processing, Vol. 18, pp.1921 - 1935, 2009.
%
% SYNOPSIS: [h]  = GCEHistogram_modification(I) 
%
% INPUT f : input image
%		B&W strech parameters : b,w and 1/(alpha+1)
%		Level of enhancement: g                      
%
% gray-level range for B&W streching is [0,b] and [w,255]
%
% OUTPUT h : modified histogram
%
% REMARKS
%
% created with MATLAB ver.: 8.1.0.604 (R2013a) on Linux 4.4.8-1-MANJARO #1 SMP PREEMPT Thu Apr 21 00:21:08 UTC 2016 x86_64
%
% created by: Daniel Paredes
% DATE: 16-May-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thrs = 0;
umin = 10;

tmp = abs(f - [zeros(size(f,1),2) f(:,1:end-2) ]) ;
k = sum(tmp(:));
htmp = tmp<=thrs;
count = sum(htmp(:));

ho = zeros(256,1);
data =f(~htmp);

% Histogram
for n=1:length(data)
    ho(data(n)+1)=ho(data(n)+1)+1;
end

% normalization
kg =k*g;
k_prime= kg/pow2(2,round(log2(kg)));

u = min(count/256,umin);

h = zeros(256,1);
for n =1:256
    if (b < n) && (n < w)
        h(n) = (1-k_prime)*u+k_prime*ho(n);
    else
        h(n) = ( (1-k_prime)*u + k_prime*ho(n)  )/(1+alpha);
    end
end

sumb = 0;
lookup = zeros(256,1);
h = h/sum(h);
for k = 1:256
    sumb = sumb + h(k);
    lookup(k)  = sumb*255 ;
end

ho = lookup;
io = lookup(f+1);

figure; imshow(f,[]); title('Input image')
figure; imshow(histeq(f/255)); title('Using histeq ')
figure; imshow(io,[]); title('GCE hist modification algorithm')
