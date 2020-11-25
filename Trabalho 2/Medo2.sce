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

L = 12;

A=200;
B=1.0;
C=4.5;
D=6.5;
sx = 1.0

u = 1;
c=0.9;

t=0;
tempo = 4;

nosEspaco = 400;

deltax = L / nosEspaco;

vetorEspaco(1) = 0

for i = 2:nosEspaco
    vetorEspaco(i) = vetorEspaco(i-1) + deltax;
end

function [oldVet,newVet] = inicializarVetores()
    for i = 1:nosEspaco
     
         value = %e^(-A*(vetorEspaco(i)-B)^2) 
         
         if C<=vetorEspaco(i) && vetorEspaco(i)<=D then
          
            value = value + sx;
          
         end
            
        oldVet(i)=value;
        newVet(i)=value;
        
    end
endfunction

function [newVet]=FTBS()
    
    newVet=[];
    oldVet=[];
    deltat= c*deltax/u;
    
    [oldVet,newVet] = inicializarVetores()
    
    while t < tempo
        
        newVet(2) = oldVet(1)
        oldVet(2) = oldVet(1)
        
        for j = 3:nosEspaco-1
        
            newVet(j) = (c)*oldVet(j-1) + (1-c)*oldVet(j);
            
        end
        
        newVet(nosEspaco-1) = oldVet(nosEspaco-1)
        oldVet(nosEspaco-1) = oldVet(nosEspaco-2)
     
    oldVet=newVet;
    t = t + deltat;
end

endfunction


function [newVet]=LaxFriedrichs()
    
   newVet=[];
   oldVet=[];
   deltat= (c*deltax)/u;

   [oldVet,newVet] = inicializarVetores()
   
    while t < tempo
        
        newVet(2) = oldVet(1)
        oldVet(2) = oldVet(1)
        
            
        for j = 3:nosEspaco-1
            
            //newVet(j) = c*oldVet(j-1) + (1-c)*oldVet(j);
            newVet(j) = (0.5-(c/2))*oldVet(j+1) + (0.5+(c/2))*oldVet(j-1);
            
        end
            
       newVet(nosEspaco-1) = oldVet(nosEspaco-1)
       oldVet(nosEspaco-1) = oldVet(nosEspaco-2)
         
       oldVet=newVet;
       t = t + deltat;
    end

endfunction


function [newVet]=LaxWendroff()
 
    newVet=[];
    oldVet=[];
    deltat= (c*deltax)/u;
    s = ((u^2)*(deltat^2))/(2*(deltax^2))

    [oldVet,newVet] = inicializarVetores()

    while t < tempo
        
        oldVet(2)=oldVet(1);
        newVet(2) = oldVet(1)
        
        for j = 3:nosEspaco-1
            
           newVet(j) = (1-2*s)*oldVet(j) + ( (c/2) + s )*oldVet(j-1) + (s - (c/2) )*oldVet(j+1);
            
        end
        
        newVet(nosEspaco-1) = oldVet(nosEspaco-1)
        oldVet(nosEspaco-1) = oldVet(nosEspaco-2)
     
        oldVet= newVet;
        t = t + deltat;
    end
    
endfunction

function [newVet]=BeamWarming()
    
    newVet=[];
    oldVet=[];
    deltat= (c*deltax)/u;
    s = ((u^2)*(deltat^2))/(2*(deltax^2))
    
    [oldVet,newVet] = inicializarVetores()

    while t < tempo
            
        
            newVet(2) = c*oldVet(1) + (1-c)*oldVet(2);
                 
            newVet(3) = newVet(2)
            oldVet(3) = oldVet(2)
            
            for j = 4:nosEspaco-1
            
               newVet(j) = (1+ s -(3*c/2))*oldVet(j) + (2*c - 2*s)*oldVet(j-1) + (s-(c/2))*oldVet(j-2);
                
            end
            
            newVet(nosEspaco-1) = oldVet(nosEspaco-1)
            oldVet(nosEspaco-1) = oldVet(nosEspaco-2)
         
        oldVet=newVet;
        t = t + deltat;
    end
   
endfunction

plot2d(vetorEspaco, FTBS(),2);
plot2d(vetorEspaco, LaxFriedrichs(),3);
plot2d(vetorEspaco, LaxWendroff(),1);
plot2d(vetorEspaco, BeamWarming(),5);
legends(['FTBS';'Lax-Friedrichs';'Lax-Wendroff';'Beam-Warming'],[2,3,1, 5],opt="lr")
title("Gráfico Espaço X Concentração",'fontsize',3);
xlabel("Cm",'fontsize',3);
ylabel("C(mol/cm^3)",'fontsize',3);
xgrid();
