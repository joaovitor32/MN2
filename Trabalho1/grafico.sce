/*n=20 //Nos no espaco
x=read("espaco.txt",-1,1)
y=read("concentracao.txt",-1,1)

plot(x, y)
title("Espaco X Concentração",'fontsize',4)
xlabel("L(cm)",'fontsize',4)
ylabel("Concentração mol/cm^3",'fontsize',4)*/

//n=100 //Nos no espaco

cini = 3;
cinj = 12;

L = 100;
alfa = 0.09;
u = 0.2;

tempo = 200;

nosEspaco = 40;
nosTempo = 20;

t=0;

deltax = L / nosEspaco;
deltat = tempo / nosTempo;
//deltat = 1 / (((2 * alfa) / (deltax * deltax)) + u / deltax);

vetorEspaco(1) = 0
for i = 2:nosEspaco
    vetorEspaco(i) = vetorEspaco(i-1) + deltax;
end

for i = 1:nosEspaco
    oldVet(i)=cini;
    newVet(i)=cini;
end


t=0;
erro = -1;

while t < tempo
        
        newVet(1) = oldVet(1) - ((deltat/deltax)*(u*(2*oldVet(1)-2*cinj)) - alfa*((oldVet(2)-3*oldVet(1)+2*cinj)/deltax));
        
        for j = 2:nosEspaco-1
           
            newVet(j) = oldVet(j) - (deltat/deltax) *(u)*(oldVet(j) -  oldVet(j - 1)) + (deltat/deltax) *(alfa)* ((oldVet(j + 1) - 2 * oldVet(j) + oldVet(j - 1))/deltax);
            
        end
        
        newVet(nosEspaco-1) = oldVet(nosEspaco-1)
        oldVet(nosEspaco-1) = oldVet(nosEspaco-2)
     
    oldVet=newVet;
    t = t + deltat;
end

plot2d(vetorEspaco, newVet,2);
title("Gráfico Espaço X Concentração",'fontsize',3);
xlabel("Cm",'fontsize',3);
ylabel("C(mol/cm^3)",'fontsize',3);
xgrid();
