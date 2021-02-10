
clf();
clear();

A=200;
B=3.0;
C=6.5;
D=9.5;
sx = 1.0
tam = 25;
u = 1;
c=0.95;
t=0;
tempo = 12;
nEsp = 30;
dx = tam / nEsp;
dt= c*dx/u;

vetorNovo=[];
vetorVelho=[];
vetEsp=[]; 

vetEsp(1) = 0
for i = 2:nEsp
    vetEsp(i) = vetEsp(i-1) + dx;
end

function [pontos]=Upwind(v1,v2)
        
    while t < tempo
      
        v1(2) = v2(1);
        v2(2) = v1(1);
        
        for it = 3:nEsp-1
            
            if (v1(it+1)- v1(it))== 0 then
                auxAvancado = 10**(-6)
            else
                auxAvancado = v1(it+1)- v1(it)
            end
            
            if  (v1(it) - v1(it-1)) ==0 then
                auxRecuado=10**(-6)
            else
                auxRecuado=v1(it) - v1(it-1)
            end 
            
            thetaRecuado = (v1(it-1)-v1(it-2))/auxRecuado
            thetaAvancado = (v1(it)-v1(it-1))/auxAvancado
            
            upwindAvancado = 0
            upwindRecuado = 0
            
            v2(it) = v1(it) - c*(v1(it)-v1(it-1)) - (c/2)*(1-c)*((upwindAvancado*(v1(it+1)-v1(it-1))-upwindRecuado*(v1(it)-v1(it-1))))
       
        end
    
        v2(nEsp-1) = v1(nEsp-1)
        v1(nEsp-1) = v1(nEsp-2)
        
        v1=v2
        t = t + dt
    
    end
    
    pontos=v2

endfunction

function [pontos]=VanAlbada(v1,v2)
        
    while t < tempo
  
        v1(2) = v2(1);
        v2(2) = v1(1);
        
        for it = 3:nEsp-1
            
            if (v1(it+1)- v1(it))== 0 then
                auxAvancado = 10**(-6)
            else
                auxAvancado = v1(it+1)- v1(it)
            end
            
            if  (v1(it) - v1(it-1)) ==0 then
                auxRecuado=10**(-6)
            else
                auxRecuado=v1(it) - v1(it-1)
            end 
            
            thetaRecuado = (v1(it-1)-v1(it-2))/auxRecuado
            thetaAvancado = (v1(it)-v1(it-1))/auxAvancado
            
            //Limitadores Van-Albada
            vanAlbadaAvancado = max(0,min(1,2*thetaAvancado),min(2,thetaAvancado))
            vanAlbadaRecuado = max(0,min(1,2*thetaRecuado),min(2,thetaRecuado))
            
            v2(it) = v1(it) - c*(v1(it)-v1(it-1)) - (c/2)*(1-c)*((vanAlbadaAvancado*(v1(it+1)-v1(it-1))-vanAlbadaRecuado*(v1(it)-v1(it-1))))
       
        end
    
        v2(nEsp-1) = v1(nEsp-1)
        v1(nEsp-1) = v1(nEsp-2)
        
        v1=v2
   
        t = t + dt
    
    end

    pontos=v2
    
endfunction


function [pontos]=Superbee(v1,v2)
        
    while t < tempo
      
        v1(2) = v2(1);
        v2(2) = v1(1);
        
        for it = 3:nEsp-1
            
            if (v1(it+1)- v1(it))== 0 then
                auxAvancado = 10**(-6)
            else
                auxAvancado = v1(it+1)- v1(it)
            end
            
            if  (v1(it) - v1(it-1)) ==0 then
                auxRecuado=10**(-6)
            else
                auxRecuado=v1(it) - v1(it-1)
            end 
            
            thetaRecuado = (v1(it-1)-v1(it-2))/auxRecuado
            thetaAvancado = (v1(it)-v1(it-1))/auxAvancado
            
            //Limitadores Superbee
            superbeeAvancado = (thetaAvancado*thetaAvancado + thetaAvancado)/(thetaAvancado*thetaAvancado + 1)
            superbeeRecuado = (thetaRecuado*thetaRecuado + thetaRecuado)/(thetaRecuado*thetaRecuado + 1)
            
            v2(it) = v1(it) - c*(v1(it)-v1(it-1)) - ((c-(c*c))/2)*((superbeeAvancado*(v1(it+1)-v1(it-1))-superbeeRecuado*(v1(it)-v1(it-1))))
       
        end
    
        v2(nEsp-1) = v1(nEsp-1)
        v1(nEsp-1) = v1(nEsp-2)
        v1=v2
        t = t + dt
    
    end
    
    pontos=v2
endfunction

function [pontos]= calcLimitedAdvec(typeLimiter)
    pontos=[]
    //Inserindo valores em vetores
    for j = 1:nEsp
     
         v = %e^(-A*(vetEsp(j)-B)^2) 
         
         if C<vetEsp(j) && vetEsp(j)<D then
          
            v = v + sx;
          
         end
            
        vetorVelho(j)=v;
        vetorNovo(j)=v;
        
   end 
     
   select typeLimiter,
   
    case 'Upwind' then pontos = Upwind(vetorVelho,vetorNovo);
    case 'Van-Albada' then pontos = VanAlbada(vetorVelho,vetorNovo);
    case 'Superbee' then pontos = Superbee(vetorVelho,vetorNovo);
           
   end
 

endfunction

plot2d(vetEsp, calcLimitedAdvec('Upwind'),5);
plot2d(vetEsp, calcLimitedAdvec('Van-Albada'),28);
plot2d(vetEsp, calcLimitedAdvec('Superbee'),25);

legends(['Upwind';'Van-Albada';'Superbee'],[5,28,25],opt="lr")
title("Espaço x Concentração",'fontsize',5);
xlabel("Metro",'fontsize',5);
ylabel("Concentração mol/M^3",'fontsize',5);
xgrid();
