function [fx] = funcionci(x)
 fx=zeros(length(x),1);
 for i=1:length(fx)
     if x(i)<0 
         fx(i)=-1;
     else if x(i)<1e-8
          fx(i)=0;
     else 
         fx(i)=1;
     end
end

end