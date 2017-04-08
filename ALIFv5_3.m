function [IMF,mask_lengths] = ALIFv5_3(f,options)

%
% Generate the decomposition of a signal f :
%
%  f = IMF(1,:) + IMF(2,:) + ... + IMF(size(IMF, 1), :)
%
% where the last component is the trend and other components are IMFs
%
%
%                    Input
%
%   f          Signal to be decomposed
%
%   options    Structure, generated using function decompSettings_v5, containing
%              all the parameters needed in the various algorithms
%
%                               Output
%
%   IMF           Matrices containg in row i the i-th IMF. The last row
%                  contains the remainder-trend.
%
%   mask_lengths  Mask length functions used for each IMF
%
%   See also SETTINGS_ALIF, IF_V5, SETTINGS_IF, GET_MASK_V1, MAXMINS_v3, PLOT_IMF_V8.
%
%  Ref: A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for 
%  Signal Decomposition and Instantaneous Frequency analysis'. Applied and 
%  Computational Harmonic Analysis, Volume 41, Issue 2, September 2016, 
%  Pages 384-411. doi:10.1016/j.acha.2016.03.001
%  ArXiv http://arxiv.org/abs/1411.6051



%% We deal with the input

if nargin == 0,  help ALIFv5_3; return; end
if nargin == 1, options = Settings_ALIF; end

extensionType = 'p'; % used in the calculations of mins and maxs

N = length(f);
if size(f,1)>size(f,2)
    f = f.';
end
if size(f,1)>1
    disp('Wrong dataset, the signal must be a single row vector')
end
IMF =[];

if options.plots>0
    if options.saveplots==1
        nameFile=input('Please enter the name of the file as a string using '' and '' <<  '); %'v081_Ex2';
    end
end

%% Main code


fprintf(['\n\n         ****************** WARNING ******************\n\n '...
    'We assume periodicity in the signal and\n\n  its instantaneous periods\n\n' ...
    '         *********************************************\n'])
pause(0.5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Adaptive Local Iterative Filtering           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

load('prefixed_double_filter','MM');
if length(f)>length(MM)
    fprintf(['\n\n      ********   Warning  *********\n\n'...
        ' The filter MM should contain more points\n'...
        ' to properly decompose the given signal\n\n'...
        ' We will use interpolation to generate a\n'...
        ' filter with the proper number of points\n\n'...
        '      *****************************\n\n'])
end

%%%%%%%%%%%%%%% Identify extreme points %%%%%%%%%%%%%%
maxmins_f=Maxmins_v3(f,extensionType);


while  length(maxmins_f) > (options.ALIF.ExtPoints) && size(IMF,1) < (options.ALIF.NIMFs) && toc < (options.maxTime)
    %% Outer loop
    
    if options.verbose>0
        fprintf('\n   ======================== IMF # %1.0d ========================\n',size(IMF,1)+1)
    end
    
    h = f;
    
    %%%%%%%%%%  Computing the mask length  %%%%%%%%%%%%%%%%
    
    if options.verbose>0
        fprintf('\n    --------------- Mask length computation ---------------\n')
    end
    
    T_f=[diff(maxmins_f) (maxmins_f(1)+N-maxmins_f(end))];
    temp_T_f=[T_f T_f T_f T_f T_f T_f T_f T_f T_f T_f T_f];
    temp_maxmins_f=[maxmins_f maxmins_f+N maxmins_f+2*N maxmins_f+3*N ...
        maxmins_f+4*N maxmins_f+5*N maxmins_f+6*N ...
        maxmins_f+7*N maxmins_f+8*N maxmins_f+9*N maxmins_f+10*N];
    temp_iT_f= interp1(temp_maxmins_f,temp_T_f,1:11*N,'cubic');
    iT_f = temp_iT_f(5*N+1:6*N);
    
    
    nTry=1;
    
    iT_f0=iT_f;
    
    if options.verbose>0
        fprintf('\n       ~~~~~~~~~~~~~~ Mask length using IF ~~~~~~~~~~~~~~\n')
    end
    
    OK=0;
    while OK==0
        opts=Settings_IF('IF.ExtPoints',3,'IF.NIMFs',nTry,'verbose',options.verbose,'IF.alpha',1,'IF.extensionType','p');
        IMF_iT_f = IF_v5(iT_f0,opts);  % We use IF algo for periodic signals to compute the mask length
        if 0>=min(IMF_iT_f(end,:)) && (size(IMF_iT_f,1)-1)==nTry
            nTry=nTry+1;
        elseif 0>=min(IMF_iT_f(end,:)) && not((size(IMF_iT_f,1)-1)==nTry)
            disp('Negative mask length')
            return
        else
            OK=1;
        end
    end
    
    if options.verbose>0
        fprintf('\n       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    end
    
    iT_f = IMF_iT_f(end,:);
    
    nn = 2*options.ALIF.xi;
    
    iT_f = nn*iT_f;
    
    if ceil(max(iT_f))>=floor(N/2)
        disp('The computation of the IMF requires a Mask length bigger than the signal itself')
        disp('From this IMF on you can try reducing the value of the paramter ALIF.Xi in Settings_ALIF')
        break
    else
    
    if options.plots>0
        figMask=figure;
        title(['Mask length IMF_{' num2str(size(IMF,1)+1) '}'])
        plot(iT_f,'b')
        hold on
        plot(maxmins_f,nn*T_f,'kx')
        plot(nn*sum(IMF_iT_f,1),'r')
        legend('Mask length','periods of the signal','Instantaneous periods')
        set(figMask,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        if options.saveplots==1
            saveas(figMask,[nameFile '_MaskLength_IMF_' num2str(size(IMF,1)+1)], 'fig')
            saveas(figMask,[nameFile '_MaskLength_IMF_' num2str(size(IMF,1)+1)], 'epsc')
            saveas(figMask,[nameFile '_MaskLength_IMF_' num2str(size(IMF,1)+1)], 'png')
        end
        pause(0.01)
    end 
    mask_lengths(size(IMF,1)+1,:)=iT_f;
     
    if options.verbose>0
        fprintf('\n    ------------------------------------------------------\n')
    end
    
    inStepN=0;
     
    W=zeros(N,N);
    for i=1:N        
        wn = get_mask_v1(MM, iT_f(i));
        wn=wn/norm(wn,1);
        len = (length(wn)-1)/2;
        wn = [reshape(wn,1,2*len+1) zeros(1,N-2*len-1)];
        W(i,:)=circshift(wn,[0 (i-len-1)]);
    end
    
    %%
    
    if options.verbose>0
        fprintf('\n IMF # %1.0d\n',size(IMF,1)+1)
        fprintf('\n  step #            SD\n\n')
    end
    
    SD=Inf;
    
    if options.plots>0
        figM=figure;
    end
    while SD > options.ALIF.delta && inStepN<500
        %% Inner loop
        
        inStepN=inStepN+1;
        
        if options.verbose>0
            
        end
        
        %%%%%%%% Compute the average %%%%%%%%%
        
        ave = W*h';
        
        %%%%%%%%%%%%%%%%%% generating h_n %%%%%%%%%%%%%%%%%%
        
        SD = norm(ave)^2/norm(h)^2;
        h = h - ave';
        
        if options.verbose>0
            fprintf('    %2.0d      %1.14f\n',inStepN,SD)
        end
        
        if options.plots>0 && rem(inStepN,1)==0
            
            textLeg{1}='f-h';
            textLeg{2}='h';
            titLeg=sprintf('IMF # %2.0d  Inner step =   %2.0d',size(IMF,1)+1,inStepN);
            
            title(sprintf('IMF # %2.0d  Inner step =   %2.0d',size(IMF,1)+1,inStepN))
            plot_imf_v8([h;f-h],1:length(h),titLeg,textLeg,figM);
            set(figM,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            
            if options.saveplots==1
                saveas(figM,[nameFile '_IMF' num2str(size(IMF,1)+1) '_step' num2str(inStepN)], 'fig')
                saveas(figM,[nameFile '_IMF' num2str(size(IMF,1)+1) '_step' num2str(inStepN)], 'epsc')
                saveas(figM,[nameFile '_IMF' num2str(size(IMF,1)+1) '_step' num2str(inStepN)], 'png')
            end
            pause(0.01)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%% Adding the new IMF %%%%%%%%%%%%%%%%%%%%%
    IMF=[IMF; h];
    
    %%%%%%%%%%%%%%%%%%%%% Updating the signal %%%%%%%%%%%%%%%%%%%%%%
    f = f-h;
    
    %%%%%%%%%%%%%%% Identify extreme points %%%%%%%%%%%%%%
    maxmins_f=Maxmins_v3(f,extensionType);
    
    if options.saveEnd == 1
        save('Decomp_ALIF_v5_3.mat')
    end
    
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%% Adding the trend %%%%%%%%%%%%%%%%%%%%%%%
IMF =[IMF; f];

if options.saveEnd == 1
    save('Decomp_ALIF_v5_3.mat')
end

end


%% Auxiliar functions

function a=get_mask_v1(y,k)
% get the mask with length 2*k+1
% k could be an integer or not an integer
% y is the area under the curve for each bar

n=length(y);
m=(n-1)/2;

if k<=m % The prefixed filter contains enough points
    
    if mod(k,1)==0     % if the mask_length is an integer
        
        a=zeros(1,2*k+1);
        
        for i=1:2*k+1
            s=(i-1)*(2*m+1)/(2*k+1)+1;
            t=i*(2*m+1)/(2*k+1);
            
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            %t2=ceil(t)-t;
            
            if floor(t)<1
                disp('Ops')
            end
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        
    else   % if the mask length is not an integer
        new_k=floor(k);
        extra = k-new_k;
        c=(2*m+1)/(2*new_k+1+2*extra);
        
        a=zeros(1,2*new_k+3);
        
        t=extra*c+1;
        t1=t-floor(t);
        %t2=ceil(t)-t;
        if k<0
            disp('Ops')
        end
        a(1)=sum(y(1:floor(t)))+t1*y(floor(t));
        
        for i=2:2*new_k+2
            s=extra*c+(i-2)*c+1;
            t=extra*c+(i-1)*c;
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        t2=ceil(t)-t;
        
        a(2*new_k+3)=sum(y(ceil(t):n))+t2*y(ceil(t));
    end
else % We need a filter with more points than MM, we use interpolation
    dx=0.01;
    % we assume that MM has a dx = 0.01, if m = 6200 it correspond to a
    % filter of length 62*2 in the physical space
    f=y/dx; % function we need to interpolate
    dy=m*dx/k;
    b=interp1(0:m,f(m+1:2*m+1),0:m/k:m);
    if size(b,1)>size(b,2)
        b=b.';
    end
    if size(b,1)>1
        fprintf('\n\nError!')
        disp('The provided mask is not a vector!!')
        a=[];
        return
    end
    a=[fliplr(b(2:end)) b]*dy;
    if abs(norm(a,1)-1)>10^-14
        %         fprintf('\n\nError\n\n')
        %         fprintf('Area under the mask = %2.20f\n',norm(a,1))
        %         fprintf('it should be equal to 1\nWe rescale it using its norm 1\n\n')
        a=a/norm(a,1);
    end
end

end


function varargout = Maxmins_v3(f,extensionType)

% Based on Version 2. 
% Modified the way zero-derivative regions are handled.
%
% Identify the maxima and minima of a signal f

if nargin == 1, extensionType = 'c'; end
N = length(f);
Maxs=zeros(1,N);
Mins=zeros(1,N);
df = diff(f);


h = 1;
cIn=0;
if strcmp(extensionType,'p') && df(1) == 0 && df(end) == 0
    while df(h)==0
        cIn=cIn+1;
        h=h+1;
    end
    if df(h) < 0
        Initial_df=-1;
    else
        Initial_df=+1;
    end
end

cmaxs=0;
cmins=0;

if strcmp(extensionType,'c') && df(1) == 0
    while df(h)==0
        h=h+1;
    end
    if df(h) < 0
        cmaxs=cmaxs+1;
        Maxs(cmaxs)=h;
    else
        cmins=cmins+1;
        Mins(cmins)=h;
    end
end

c = 0;

last_df=[];
for i=h:N-2
    if   df(i)*df(i+1) == 0
        if df(i) < 0
            last_df=-1;
            posc = i;
        elseif df(i) > 0
            last_df=+1;
            posc = i;
        end
        c = c + 1;
        if df(i+1) < 0
            if last_df==+1
                cmaxs=cmaxs+1;
                Maxs(cmaxs)=posc+floor((c-1)/2)+1;                
            end
            c=0;
        end
        if df(i+1) > 0
            if last_df==-1
                cmins=cmins+1;
                Mins(cmins)=posc+floor((c-1)/2)+1;                
            end
            c=0;
        end
        
    end
    if   df(i)*df(i+1) < 0
        if df(i) < df(i+1)
            cmins=cmins+1;
            Mins(cmins)=i+1;
            last_df=-1;
        else
            cmaxs=cmaxs+1;
            Maxs(cmaxs)=i+1;
            last_df=+1;
        end
    end
end
if c > 0
    if strcmp(extensionType,'p')
        disp('Code to be completed!')
        if Initial_df < 0
            if last_df==+1
                cmaxs=cmaxs+1;
                Maxs(cmaxs)=mod(posc+floor((c+cIn-1)/2)+1,N);                
            end            
        elseif Initial_df > 0
            if last_df==+1
                cmins=cmins+1;
                Mins(cmins)=mod(posc+floor((c+cIn-1)/2)+1,N);                 
            end   
        else
            disp('Code missing!')
        end        
    end        
    if strcmp(extensionType,'c')        
        if last_df > 0
            cmaxs=cmaxs+1;
            Maxs(cmaxs)=posc+1;
        else
            cmins=cmins+1;
            Mins(cmins)=posc+1;
        end
    end 
    if Mins(cmins)==0
        Mins(cmins)=N;
    end
    if Maxs(cmaxs)==0
        Maxs(cmaxs)=N;
    end
end

Maxs=Maxs(1:cmaxs);
Mins=Mins(1:cmins);
maxmins=sort([Maxs Mins]);

if strcmp(extensionType,'p') % we deal with a periodic signal
    disp('Code to be completed')
    if isempty(maxmins)
        maxmins = 1;
    else
        if maxmins(1)~=1 && maxmins(end)~=N
            if (f(maxmins(end)) > f(maxmins(end)+1) && f(maxmins(1)) > f(maxmins(1)-1)) || (f(maxmins(end)) < f(maxmins(end)+1) && f(maxmins(1)) < f(maxmins(1)-1))
                maxmins=[1 maxmins];
            end
        end
    end
elseif strcmp(extensionType,'c')
    if not(isempty(maxmins))
        if maxmins(1) ~= 1 && maxmins(end) ~= N && df(1)~=0 && df(end)~=0
            if Maxs(1) < Mins(1)
                Mins=[1 Mins];
            else
                Maxs=[1 Maxs];
            end
            if Maxs(end) < Mins(end)
                Maxs=[Maxs N];
            else
                Mins=[Mins N];
            end
            maxmins = [1, maxmins, N];
        elseif maxmins(1) ~= 1 && df(1)~=0
            maxmins = [1, maxmins];
            if Maxs(1) < Mins(1)
                Mins=[1 Mins];
            else
                Maxs=[1 Maxs];
            end
        elseif  maxmins(end) ~= N && df(end)~=0
            maxmins = [maxmins, N];
            if Maxs(end) < Mins(end)
                Maxs=[Maxs N];
            else
                Mins=[Mins N];
            end
        end
    end
elseif strcmp(extensionType,'r')
    disp('Code to be completed')
    if isempty(maxmins)
        maxmins = [1, N];
    else
        if maxmins(1) ~= f(1) && maxmins(end) ~= f(end)
            maxmins = [1, maxmins, N];
        elseif f(maxmins(1)) ~= f(1)
            maxmins = [1, maxmins];
        elseif  f(maxmins(end)) ~= f(end)
            maxmins = [maxmins, N];
        end
    end
end

if nargout<=1
    varargout{1}=maxmins;
elseif nargout==2
    varargout{1}=Maxs;
    varargout{2}=Mins;
end

end