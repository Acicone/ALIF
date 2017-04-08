%  Ref. Example 3 page 17 
%  A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for 
%  Signal Decomposition and Instantaneous Frequency analysis'. Applied and 
%  Computational Harmonic Analysis, Volume 41, Issue 2, September 2016, 
%  Pages 384-411. doi:10.1016/j.acha.2016.03.001
%  ArXiv http://arxiv.org/abs/1411.6051


T=2*pi; % Period
N=5000; % Number of sample points
D=32;   % Delta in the frequencies
dt=T/N;
t=0:dt:T;

x=cos(-1/2*D/T*t.^2-4*t);


x1=[fliplr(x(2:end-1)) x ];
  
tt=[-fliplr(t(2:end-1)) t];  
figure; plot(tt,[x1])

y=cos(-1/2*D/T*t.^2-20*t);

ty1=0:dt:2*T-dt;
y1=[fliplr(y(2:end-1)) y ];
 
figure; plot(tt,[y1])

f=x1+y1+1;

figure; 
plot(tt,f)
set(gca,'fontsize', 20);
axis([-T T -1 3])


%% ALIF decomposition

opt = Settings_ALIF('ALIF.NIMFs',1,'plots',1,'saveplots',0,'ALIF.xi',2,'ALIF.delta',4*10^-6);

[IMF,mask_lengths]=ALIFv5_3(f,opt);

plot_imf_v8(IMF-[y1;x1+1],tt)

plot_imf_v7(IMF,[y1;x1+1],tt,'IMFs','Ground truth');


