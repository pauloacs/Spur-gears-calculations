# -*- coding: utf-8 -*-
"""
@author: Paulo Sousa
"""
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import fsolve

print('                      *Representação gráfica da pressão de contacto**')
print('                             -Indicar os seguintes parâmetros-')

    #Definição da função que condensa o código utilizado no programa 
    #Correção_de_dentado considerando o parâmetro K. Neste programa,
    #será necessária para obter o parâmetros geométricos após correção.
    
def corr(alfa,m,z2,z1)  : 
    
    
    a=(z1+z2)*m/2
    def invol(alfa):
        alfa_rad=alfa*math.pi/180
        return math.tan(alfa_rad)-alfa_rad
    
    
    def escorregamentos(a,x1,x2,alfa,alfa_linha,z1,z2,K):
        alfa_rad,alfa_linha_=alfa*math.pi/180 , alfa_linha*math.pi/180
        ra1,ra2=z1*m/2+(1+x1)*m-K*m, z2*m/2+(1+x2)*m-K*m
        rb1,rb2=z1*m/2*math.cos(alfa_rad), z2*m/2*math.cos(alfa_rad)
    
        raiz1,raiz2=(ra1**2-rb1**2)**0.5, (ra2**2-rb2**2)**0.5
            
        Gs1B_max=abs(1-raiz2/(a*math.sin(alfa_linha_)-raiz2)*(z1/z2))
        Gs2A_max=abs((raiz1/(a*math.sin(alfa_linha_)-raiz1))*(z2/z1)-1)
        return Gs1B_max,Gs2A_max
    
    
    if z1+z2<60: z2_=60-z1
    else: z2_=z2
    def eq(x):
        a_=(z1+z2_)*m/2
        s1B,s2A=escorregamentos(a_,x,-x,alfa,alfa,z1,z2_,0)
        return s1B-s2A
    x1=fsolve(eq,0.3)
    
    if z1+z2<60:
        def sist(z):
            alfa_rad=alfa*math.pi/180
            x2, a_linha, alfa_linha=z[0],z[1], z[2]
            a=(z1+z2)*m/2
            F0=a_linha*math.cos(alfa_linha*math.pi/180)-a*math.cos(alfa_rad)
            F1=invol(alfa)+2*math.tan(alfa_rad)*(x1+x2)/(z1+z2)-invol(alfa_linha)
            a1=(z1+z2)*m/2+(x1+x2)*m
            K=(a1-a_linha)/m
            Gs1B_linha,Gs2A_linha=escorregamentos(a_linha,x1,x2,alfa,alfa_linha,z1,z2,K)
            F2=Gs1B_linha-Gs2A_linha
            return [F0,F1,F2]
        x2,a_linha,alfa_linha=fsolve(sist,[0,a,alfa])
        a1=(z1+z2)*m/2+(x1+x2)*m
        K=(a1-a_linha)/m
        Gs1B,Gs2A=escorregamentos(a_linha,x1,x2,alfa,alfa_linha,z1,z2,K)
        
    else:
        x2=-x1
        alfa_linha=alfa
        a_linha=a
    return x1,x2,alfa_linha,a_linha

    #Input dos valores necessários:

alfa=int(input('Introduza o ângulo de pressão (em graus): '))
N=float(input('Introduza o valor da carga normal(N): '))
b=float(input('Largura do dente(mm): ') )
m=float(input('Módulo (mm): ') )
z1=float(input('Número de dentes do pinhão (Z1) : '))
z2=float(input('Número de dentes da roda (Z2) : '))

x1,x2,alfa,a=corr(alfa,m,z2,z1) #Obtenção dos parâmetros após correção de dentado 
                                #aplicando o método mais adequado para cada caso.

Fn=N/b    #Carga normal/largura do dente (N/mm)
alfa_rad=alfa*math.pi/180
r1,r2=z1*m/2,z2*m/2
T1T2=(a)*math.sin(alfa_rad)
rb1,rb2=r1*math.cos(alfa_rad),r2*math.cos(alfa_rad)
ra1,ra2=r1+m+x1*m,r2+m+x2*m
T2A=T1T2-(ra1**2-rb1**2)**0.5
T2B=(ra2**2-rb2**2)**0.5
pb=math.pi*m*math.cos(alfa_rad) #Passo de base
T2V=T2B-pb
T2W=T2A+pb

rho_R=[]
sigma_H1=[]
sigma_H2=[]
Rho_2=np.arange(0.01,T1T2+0.1,0.1)

    #Função para o cálculo da pressão de contacto Hertziana: 
    #(No caso específico de E=210000 Nmm-2, e sendo o coeficiente de Poisson  0,3)
    
def sigma(Fn,rho_2):
    T1T2=(a)*math.sin(alfa_rad)
    rho_1=T1T2-rho_2
    rho_r=1/(1/rho_1+1/rho_2)
    sigma=192*(Fn/rho_r)**0.5
    return rho_r,sigma

for rho_2 in Rho_2:
    rho_r,sigma_1=sigma(Fn,rho_2) #Cálculo da pressão de contacto nas zonas 
                                  #em que apenas um par de dentes está em contacto.
                                  
    rho_r,sigma_2=sigma(Fn/2,rho_2) #Pressão de contacto quando estão 2 pares de dentes em contacto.
                                    
    rho_R.append(rho_r)
    sigma_H1.append(sigma_1)                 #Vem em [N/mm^2]=[MPa]
    sigma_H2.append(sigma_2)                 #Vem em [N/mm^2]=[MPa]

plt.figure(dpi=200) #Definição da resolução da imagem que vai ser criada em seguida:

plt.plot(Rho_2,sigma_H1, lw = 1.5,label='1 par de dentes em contacto ') 
plt.plot(Rho_2,sigma_H2, lw = 1.5,label='2 pares de dentes em contacto ')
plt.legend(frameon = True,
           loc = 'upper center',
           fontsize = 8)
T2I=rb2*math.tan(alfa_rad) #Relação trigonométrica para obtenção da de T2I tendo em conta
                           #a geometria após correção.

    #Representação de linhas verticais auxiliares na interpretação do gráfico:
    
plt.axvline(x=T2I, ymin=0,ymax=sigma(Fn,T2I)[1]/sigma_H2[-5], linewidth=0.7,ls = '--',color='b')
plt.axvline(x=T2A, ymin=0, ymax=sigma(Fn/2,T2A)[1]/sigma_H2[-5],linewidth=1, color='orange')
plt.axvline(x=T2B, ymin=0, ymax=sigma(Fn/2,T2B)[1]/sigma_H2[-5],linewidth=1, color='orange')
plt.axvline(x=T2V,ymin=0, ymax=sigma(Fn/2,T2V)[1]/sigma_H2[-5],linewidth=0.7,ls = '--', color='b')
plt.axvline(x=T2W, ymin=0,ymax=sigma(Fn/2,T2W)[1]/sigma_H2[-5], linewidth=0.7,ls = '--',color='b')
plt.axvline(x=T2V,ymin=sigma(Fn/2,T2V)[1]/sigma_H2[-5], ymax=sigma(Fn,T2V)[1]/sigma_H2[-5],linewidth=1, color='b')
plt.axvline(x=T2W, ymin=sigma(Fn/2,T2W)[1]/sigma_H2[-5], ymax=sigma(Fn,T2W)[1]/sigma_H2[-5],linewidth=1, color='b')

    #Marcação dos pontos relevantes sobre o eixo das abcissas:
plt.xticks([0, T2A, T2B, T2V, T2I ,T2W, T1T2], 
           ['T2', 'A', 'B', 'V','I', 'W', 'T1'],
           fontsize = 8)   
    
    #Utilização de anotações para dar relevância aos pontos importantes:
    
plt.annotate('\u03C3[I] ='+str(np.round(sigma(Fn,T2I)[1],2)) + '[MPa]', xy =(T2I, sigma(Fn,T2I)[1]) , xycoords = 'data', xytext = (-90, +30),  textcoords = 'offset points',  fontsize = 5, color = 'r',arrowprops = dict(arrowstyle = '-|>', color = 'r',connectionstyle = 'arc3,rad=.2'))
plt.annotate('\u03C3[W] =' +str(np.round(sigma(Fn,T2W)[1],2)) + '[MPa]', xy =(T2W, sigma(Fn,T2W)[1]) , xycoords = 'data', xytext = (30, 30),  textcoords = 'offset points',  fontsize = 5, color = 'r',arrowprops = dict(arrowstyle = '-|>', color = 'r',connectionstyle = 'arc3,rad=.2'))

plt.xlabel('<--Reta de engrenamento-->', fontsize = 10, labelpad = 10,color="Blue")
plt.ylabel('\u03C3H '+ '[MPa]', fontsize = 10, labelpad = 20, color="Blue")
plt.axis([0,T1T2,0,sigma_H2[-5]])
plt.title('Tensão de contacto'.upper(), y = 1.1, fontsize = 14, color = 'purple', fontweight = 'bold')
plt.legend()
plt.show()
