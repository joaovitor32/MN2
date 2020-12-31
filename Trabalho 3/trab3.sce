/*n=20 //Nos no espaco
x=read("espaco.txt",-1,1)
y=read("concentracao.txt",-1,1)

plot(x, y)
title("Espaco X Concentração",'fontsize',4)
xlabel("L(cm)",'fontsize',4)
ylabel("Concentração mol/cm^3",'fontsize',4)*/

//n=100 //Nos no espaco
clf();
clear();

L =25;

A=200;
B=3.0;
C=6.5;
D=9.5;
sx = 1.0

u = 1;
c=0.8;

t=0;
tempo = 15;

nosEspaco = 500;

deltax = L / nosEspaco;

vetorEspaco(1) = 0

for i = 2:nosEspaco
    vetorEspaco(i) = vetorEspaco(i-1) + deltax;
end

function [oldVet, newVet]=inicializarVetores()
    for i = 1:nosEspaco
     
         value = %e^(-A*(vetorEspaco(i)-B)^2) 
         
         if C<=vetorEspaco(i) && vetorEspaco(i)<=D then
          
            value = value + sx;
          
         end
            
        oldVet(i)=value;
        newVet(i)=value;
        
    end
endfunction

function [psi]=psiUpwind(ponto)
  
    psi = 0;
     
endfunction

function [psi]=psiSuperbee(ponto)
  
    psi = max(0,min(1,2*ponto),min(2,ponto));
     
endfunction

function [psi]=psiValAlbada(ponto)
    
    psi = (((ponto^2) + ponto)/((ponto^2)+1));
     
endfunction

function [newVet]=Upwind()
    
    newVet=[];
    oldVet=[];
    deltat= c*deltax/u;
    tolerancia = (1 * 10) ** -6  
    
    [oldVet,newVet] = inicializarVetores()
    
    while t < tempo
        
        newVet(2) = oldVet(2) - c*(oldVet(2)-oldVet(1))
 
        oldVet(3) = oldVet(2)
        newVet(3) = oldVet(2)
        
        for j = 4:nosEspaco-1
            
             if(oldVet(j+1)-oldVet(j))==0 then
             
                aux=(10**(-12))
             
             else
             
                aux=(oldVet(j+1)-oldVet(j)) 
             
             end
            
            if (oldVet(j + 1) - oldVet(j)) <= tolerancia  then 
                 

                pontoPsiAvancado = (oldVet(j)-oldVet(j-1))/aux
            
            else
                  
                pontoPsiAvancado = (oldVet(j)-oldVet(j-1))/aux
            
            end     
        
            if (oldVet(j) - oldVet(j-1))==0 then
             
                aux=(10**(-12))
             
            else
             
                aux=(oldVet(j)-oldVet(j-1)) 
             
            end

            if  (oldVet(j) - oldVet(j-1)) <= tolerancia then 
                 
                pontoPsiRecuado = (oldVet(j-1)-oldVet(j-2))/aux
            
            else
                   
                pontoPsiRecuado = (oldVet(j-1)-oldVet(j-2))/(oldVet(j)-oldVet(j-1))
            
            end      
                 
            ponto1 = psiUpwind(pontoPsiAvancado)*(oldVet(j+1)-oldVet(j))
            ponto2 = psiUpwind(pontoPsiRecuado)*(oldVet(j)-oldVet(j-1))
        
            newVet(j) = oldVet(j) - c*(oldVet(j)-oldVet(j-1)) - (c/2)*(1-c)*(ponto1 - ponto2)
            
        end
       
        newVet(nosEspaco-1) = oldVet(nosEspaco-1)
        oldVet(nosEspaco-1) = oldVet(nosEspaco-2)
     
        oldVet=newVet;
        t = t + deltat;
    end

endfunction

function [newVet]=Superbee()
    
    newVet=[];
    oldVet=[];
    deltat= c*deltax/u;
    tolerancia = (1 * 10) ** -6  
    
    [oldVet,newVet] = inicializarVetores()
    
    while t < tempo
        
        newVet(2) = oldVet(2) - c*(oldVet(2)-oldVet(1))
 
        oldVet(3) = oldVet(2)
        newVet(3) = oldVet(2)
        
        for j = 4:nosEspaco-1
            
            if(oldVet(j+1)-oldVet(j))==0 then
             
                aux=(10**(-12))
             
             else
             
                aux=(oldVet(j+1)-oldVet(j)) 
             
             end
            
            if (oldVet(j + 1) - oldVet(j)) <= tolerancia  then 
                 

                pontoPsiAvancado = (oldVet(j)-oldVet(j-1))/aux
            
            else
                  
                pontoPsiAvancado = (oldVet(j)-oldVet(j-1))/aux
            
            end     
        
            if (oldVet(j) - oldVet(j-1))==0 then
             
                aux=(10**(-12))
             
            else
             
                aux=(oldVet(j)-oldVet(j-1)) 
             
            end

            if  (oldVet(j) - oldVet(j-1)) <= tolerancia then 
                 
                pontoPsiRecuado = (oldVet(j-1)-oldVet(j-2))/aux
            
            else
                   
                pontoPsiRecuado = (oldVet(j-1)-oldVet(j-2))/aux
            
            end   
                     
            ponto1 = psiSuperbee(pontoPsiAvancado)*(oldVet(j+1)-oldVet(j))
            ponto2 = psiSuperbee(pontoPsiRecuado)*(oldVet(j)-oldVet(j-1))
        
            newVet(j) = oldVet(j) - c*(oldVet(j)-oldVet(j-1)) - (c/2)*(1-c)*(ponto1 - ponto2)
            
        end
       
        newVet(nosEspaco-1) = oldVet(nosEspaco-1)
        oldVet(nosEspaco-1) = oldVet(nosEspaco-2)
     
        oldVet=newVet;
        t = t + deltat;
    end

endfunction

function [newVet]=vanAlbada()
    
    newVet=[];
    oldVet=[];
    deltat= c*deltax/u;
    tolerancia = (1 * 10) ** -6  
    
    [oldVet,newVet] = inicializarVetores()
    
    while t < tempo
        
        newVet(2) = oldVet(2) - c*(oldVet(2)-oldVet(1))
 
        oldVet(3) = oldVet(2)
        newVet(3) = oldVet(2)
        
        for j = 3:nosEspaco-1
                
             if(oldVet(j+1)-oldVet(j))==0 then
             
                aux=(10**(-12))
             
             else
             
                aux=(oldVet(j+1)-oldVet(j)) 
             
             end
            
            if (oldVet(j + 1) - oldVet(j)) <= tolerancia  then 
                 

                pontoPsiAvancado = (oldVet(j)-oldVet(j-1))/aux
            
            else
                  
                pontoPsiAvancado = (oldVet(j)-oldVet(j-1))/aux
            
            end     
        
            if (oldVet(j) - oldVet(j-1))==0 then
             
                aux=(10**(-12))
             
            else
             
                aux=(oldVet(j)-oldVet(j-1)) 
             
            end

            if  (oldVet(j) - oldVet(j-1)) <= tolerancia then 
                 
                pontoPsiRecuado = (oldVet(j-1)-oldVet(j-2))/aux
            
            else
                   
                pontoPsiRecuado = (oldVet(j-1)-oldVet(j-2))/aux
            
            end   
               
            ponto1 = psiValAlbada(pontoPsiAvancado)*(oldVet(j+1)-oldVet(j))
            ponto2 = psiValAlbada(pontoPsiRecuado)*(oldVet(j)-oldVet(j-1))
        
            newVet(j) = oldVet(j) - c*(oldVet(j)-oldVet(j-1)) - (c/2)*(1-c)*(ponto1 - ponto2)
            
        end
       
        newVet(nosEspaco-1) = oldVet(nosEspaco-1)
        oldVet(nosEspaco-1) = oldVet(nosEspaco-2)
     
        oldVet=newVet;
        t = t + deltat;
    end

endfunction

plot2d(vetorEspaco, Upwind(),2);
plot2d(vetorEspaco, Superbee(),3);
plot2d(vetorEspaco, vanAlbada(),19);

legends(['Upwind';'Superbee';'Van-Albada'],[2,3,19],opt="lr")
title("Gráfico Espaço X Concentração",'fontsize',3);
xlabel("M",'fontsize',3);
ylabel("C(mol/M^3)",'fontsize',3);
xgrid()

