function[h]=plot_imf_v8(imf,T,titLeg,legText,h,TextSize,date)

[m,n]=size(imf);
if nargin == 1, T=1:n; end

if isempty(T)
     T=1:n;
end

if nargin <6
    TextSize=30;
end

if isempty(TextSize)
    TextSize=30;
end

if nargin <7
    date=-1;
end

if nargin <5
    for j=1:ceil(m/5)
        h(j)=figure;
    end
end

if length(h)<ceil(m/5)
    for i= length(h)+1:ceil(m/5)
        h(i)=figure;
    end
end

for j=1:ceil(m/5)
    figure(h(j));
    if j<=floor(m/5)
        for i=1:5
            subplot(5,1,i);
            plot(T,imf(i+(j-1)*5,:),'k','LineWidth',2);
            set(gca,'fontsize', TextSize);
            if date>=0
                datetickzoom('x',date)
            end
            if nargin == 4 && not(isempty(legText))
                hh=legend(legText{i+(j-1)*5},'Location','NorthEastOutside');
                set(hh,'Interpreter','latex')
            end
            
            if CheckForIMF(imf(i+(j-1)*5,:))
                hold on;
                plot(T,zeros(1,n),'r','LineWidth',2)
                hold off
                
                if not(isnan(max(abs(imf(i+(j-1)*5,:))))) && not(max(abs(imf(i+(j-1)*5,:)))==0)
                    axis([T(1) T(end) -1.2*max(abs(imf(i+(j-1)*5,:))) 1.2*max(abs(imf(i+(j-1)*5,:)))])
                end
            else
                if (min(imf(i+(j-1)*5,:))==0 && max(imf(i+(j-1)*5,:))==0) || (isnan(max(imf(i+(j-1)*5,:))) || isnan(min(imf(i+(j-1)*5,:))))
                    axis([T(1) T(end) -1 1])
                else
                    axis([T(1) T(end) (min(imf(i+(j-1)*5,:))-0.1*abs(min(imf(i+(j-1)*5,:)))) (max(imf(i+(j-1)*5,:))+0.1*abs(max(imf(i+(j-1)*5,:))))])
                end
            end
            pause(0.1)
        end
    else
        for i=1:rem(m,5)
            subplot(rem(m,5),1,i);
            plot(T,imf(i+(j-1)*5,:),'k','LineWidth',2);
            set(gca,'fontsize', TextSize);
            if date>=0
                datetickzoom('x',date)
            end
            if nargin == 4 && not(isempty(legText))
                hh=legend(legText{i+(j-1)*5},'Location','NorthEastOutside');
                set(hh,'Interpreter','latex')
            end
            
            if CheckForIMF(imf(i+(j-1)*5,:))
                hold on;
                plot(T,zeros(1,n),'r','LineWidth',2)
                hold off
                
                if not(isnan(max(abs(imf(i+(j-1)*5,:))))) && not(max(abs(imf(i+(j-1)*5,:)))==0)
                    axis([T(1) T(end) -1.2*max(abs(imf(i+(j-1)*5,:))) 1.2*max(abs(imf(i+(j-1)*5,:)))])
                end
            else
                if (min(imf(i+(j-1)*5,:))==0 && max(imf(i+(j-1)*5,:))==0) || (isnan(max(imf(i+(j-1)*5,:))) || isnan(min(imf(i+(j-1)*5,:))))
                    axis([T(1) T(end) -1 1])
                else
                    axis([T(1) T(end) (min(imf(i+(j-1)*5,:))-0.1*abs(min(imf(i+(j-1)*5,:)))) (max(imf(i+(j-1)*5,:))+0.1*abs(max(imf(i+(j-1)*5,:))))])
                end
            end
            pause(0.1)
        end
    end
end

if nargin >= 3  && not(isempty(titLeg))
    set(gcf,'NextPlot','add');
    axes;
    tt= title(titLeg,'FontSize',TextSize);
    set(gca,'Visible','off');
    set(tt,'Visible','on');
end


end

function [YoN,checkAve]=CheckForIMF(f)

maxmins = Maxmins(f,'No');
if length(maxmins)<2
    YoN=1==0;
    checkAve=[];
else
    if f(maxmins(1))>f(maxmins(2))
        checkAve1=sum(f(maxmins(1:2:end))>0)/length(maxmins(1:2:end));
        checkAve2=sum(f(maxmins(2:2:end))<0)/length(maxmins(1:2:end));
    else
        checkAve1=sum(f(maxmins(2:2:end))>0)/length(maxmins(2:2:end));
        checkAve2=sum(f(maxmins(1:2:end))<0)/length(maxmins(1:2:end));
    end
    
    if checkAve1 >= 0.5 && checkAve2 >= 0.5
        YoN=1==1;
    else
        YoN=1==0;
    end
    checkAve=[checkAve1 checkAve2];
end

end

function maxmins = Maxmins(f,extensionType)

if nargin == 1, extensionType = 'p'; end

N = length(f);
maxmins=zeros(1,N);
df = diff(f);


h = 1;
cIn=0;
if strcmp(extensionType,'p') && df(1) == 0 && df(end) == 0
    while df(h)==0
        cIn=cIn+1;
        h=h+1;
    end
end

c = 0;
cmaxmins=0;
for i=h:N-2
    if   df(i)*df(i+1) <= 0
        if df(i+1) == 0
            if c == 0
                posc = i;
            end
            c = c + 1;
        else
            if c > 0
                cmaxmins=cmaxmins+1;
                maxmins(cmaxmins)=posc+floor((c-1)/2)+1;
                c = 0;
            else
                cmaxmins=cmaxmins+1;
                maxmins(cmaxmins)=i+1;
            end
        end
    end
end
if c > 0
    cmaxmins=cmaxmins+1;
    maxmins(cmaxmins)=mod(posc+floor((c+cIn-1)/2)+1,N);
    if maxmins(cmaxmins)==0
        maxmins(cmaxmins)=N;
    end
end

maxmins=maxmins(1:cmaxmins);

if strcmp(extensionType,'p') % we deal with a periodic signal
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
elseif strcmp(extensionType,'r')
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

end
