funcprot(0)
function [diferencial]=diff(M,h)
    
    diferencial=M;
    diferencial(1)=((M(2))-M(1))/h;
for i=2:max(size(M))-1
    diferencial(i)=((M(i+1))-(M(i-1)))/(2*h);
end
    diferencial(max(size(M)))=(M(max(size(M)))-(M(max(size(M-1)))))/(h);
    
endfunction 
