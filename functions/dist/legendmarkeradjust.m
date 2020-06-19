function legendmarkeradjust(varargin)

a1=ver('matlab');

if strcmp(a1.Release,'(R2016a)') | strcmp(a1.Release,'(R2016b)') | strcmp(a1.Release,'(R2017a)') | strcmp(a1.Release,'(R2015b)') | strcmp(a1.Release,'(R2015a)') | strcmp(a1.Release,'(R2014b)')

    
leg=get(legend);

legfontsize=leg.FontSize;
legstrings=leg.String;
legloc=leg.Location;

delete(legend)

[l1,l2,l3,l4]=legend(legstrings);


for n=1:length(l2)
    if sum(strcmp(properties(l2(n)),'MarkerSize'))
    l2(n).MarkerSize=varargin{1};
    end
end

for n=1:length(l2)
        if sum(strcmp(properties(l2(n)),'FontSize'))
        l2(n).FontSize=legfontsize;
        end
end

set(l1,'location',legloc)

else
    
s=get(legend);
s1=s.Children;
s2=[];

s2=findobj(s1,{'type','patch','-or','type','line'});

switch length(varargin)
    case 1
        marksize=varargin{1};
        for m=1:length(s2)
            set(s2(m),'markersize',marksize);
        end

    case 2
        marksize=varargin{1};
        lwidth=varargin{2};

        for m=1:length(s2)
            set(s2(m),'markersize',marksize,'linewidth',lwidth);
        end
end

end