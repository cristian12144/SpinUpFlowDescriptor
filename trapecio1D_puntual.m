function lzm=trapecio1D_puntual(nd,tf,ht,lz)

sum=0;   

       
    for j=2:nd-1
        
        sum=sum+lz(j); % Esto es por el t�rmino de sumatoria que hay en la f�rmula del c�lculo de la integral de forma num�rica.
        
    end
    
    lzm=1/tf*ht/2*(lz(1)+2*sum+lz(nd));
    
