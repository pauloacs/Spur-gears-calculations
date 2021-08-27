# -*- coding: utf-8 -*-
"""
@author: Paulo Sousa
"""

from math import sin, cos, tan, pi
from scipy.optimize import fsolve
#Bibliotecas necessárias importadas nas linhas anteriores.

print('                            ***Correção de dentado***')
print('                      --Indique os seguintes parâmetros--')

#Introdução dos parâmetros:

alfa=float(input('Ângulo de pressão: '))
m=float(input('Módulo(m): '))
z2=float(input('Nº de dentes da roda (Z2): '))
z1=float(input('Nº de dentes do pinhão (Z1): '))

    #Definição da expressão que calcula os valores de escorregamento específico máximo:

def escorregamentos(a,x1,x2,alfa,alfa_linha,z1,z2,K):
        alfa_rad,alfa_linha_=alfa*pi/180 , alfa_linha*pi/180
        ra1,ra2=z1*m/2+(1+x1)*m-K*m, z2*m/2+(1+x2)*m-K*m
        rb1,rb2=z1*m/2*cos(alfa_rad), z2*m/2*cos(alfa_rad)
    
        raiz1,raiz2=(ra1**2-rb1**2)**0.5, (ra2**2-rb2**2)**0.5
            
        Gs1B_max=abs(1-raiz2/(a*sin(alfa_linha_)-raiz2)*(z1/z2))
        Gs2A_max=abs((raiz1/(a*sin(alfa_linha_)-raiz1))*(z2/z1)-1)
        return Gs1B_max,Gs2A_max
Gs1b_max, Gs2A_max=escorregamentos((z1+z2)*m/2,0,0,alfa,alfa,z1,z2,0) 
    #Cálculo dos escorregamentos específicos

    #Por motivos de comparação, apresentam-se os valores ainda sem correção:
print('\n\nSem correção os escorregamentos específicos vêm:')
print('Gs1B_max=%.5f e Gs2A_max=%.5f' %(Gs1b_max, Gs2A_max))
print()

x=input('Considerar K?\n<y or n>')

if x=='n' or x=='N': # O programa contempla a opção do utilizador em aplicar ou não
                     #a correcção de altura 'K' que visa manter a distância entre 
                     #a cabeça do dente de uma roda, e o pé do dente da roda engrenada.

    a=(z1+z2)*m/2 #entre-eixo 
    
    def invol(alfa):  #Definição da função involuta
        alfa_rad=alfa*pi/180
        return tan(alfa_rad)-alfa_rad
    
    if z1+z2<60: z2_=60-z1 #No caso da situação de z1+z2 pequeno, calcula-se x1 considerando z2*=60-z1 
                           #Neste caso x1 e x2 não serão simétricos
                           
    else: z2_=z2 #Se z1+z2 suficientemente alto, procede-se da forma mais direta
                 #Neste caso x1=-x2
    
    def eq(x):      #Igualando os escorregamentos específicos máximos considerando correção simétrica
        
        a_=(z1+z2_)*m/2
        s1B,s2A=escorregamentos(a_,x,-x,alfa,alfa,z1,z2_,0)
        return s1B-s2A
    xSOL=fsolve(eq,0) #Utilizando a função fsolve que foi, no início, importada da biblioteca 
                      #scipy.optimize, é possível, por resolução numérica obter o valor da correção
                      #Note-se que, para o fazer, é necessário definir uma iteração inícial,
                      #e deve-se, portanto, garantir que é obtida convergência.
                      #Tendo isso em conta, apresenta-se como valor inicial: 0 , visto que
                      #Se verifica que os valores da correção andam em torno deste,
                      #e testanto o programa não foram verificadas impossibilidades.
    
    
    #No caso de z1+z2 suficientemente alto (pelo menos 60), está, neste ponto, resolvida a questão da correção 
    #na roda e no pinhão(sendo estas simétricas), nesta etapa já se poderia aprensentar 
    #x1=xSOL  e x2=-x1
    
    if z1+z2<60: #Para o caso de z1+z2 pequeno:
                 #Considera-se o procedimento de Henriot, que será o desenvolvido em seguida:
                 
        x1=xSOL
        
        #Note-se que neste procedimento, tanto o entre-eixo como o ângulo de pressão,
        #virão alterados, é necessário, portanto, conhecido o valor de x1, resolver 
        #um sistema de equações com x2, a'(novo entre-eixo), e alfa' (novo ângulo de pressão)
        #como incógnitas
        
        def sist(z): 
                    
            alfa_rad=alfa*pi/180
            x2, a_linha, alfa_linha=z[0],z[1], z[2]
            a=(z1+z2)*m/2
            F0=a_linha*cos(alfa_linha*pi/180)-a*cos(alfa_rad)
            F1=invol(alfa)+2*tan(alfa_rad)*(x1+x2)/(z1+z2)-invol(alfa_linha)       
            Gs1B_linha,Gs2A_linha=escorregamentos(a_linha,x1,x2,alfa,alfa_linha,z1,z2,0)
            F2=Gs1B_linha-Gs2A_linha
            return [F0,F1,F2]
        
        #De forma semelhante ao caso anterior, a resolução será feita numericamente,
        #e novamente se realça a importância da primeira iteração.
        #Neste caso considera-se x2=0 pelo motivo apresentado para o caso anterior, 
        #e a' e alfa' iguais aos valores de a e alfa sem correção, visto que estes variam relativamente pouco
        
        x2,a_linha,alfa_linha=fsolve(sist,[0,a,alfa])
        
        #Aqui torna-se importante perceber que a convergência do método não será possível
        #para valores muito pouco usuais(fora da gama considerada pelo representação gráfica de Henriot)
        #para o lado esquerdo.
        
        Gs1B,Gs2A=escorregamentos(a_linha,x1,x2,alfa,alfa_linha,z1,z2,0)
        print("\n\nSituação Z1+Z2<60:\nExige novo valor de entre-eixo(a') e ângulo de pressão(alfa');")
        print("x1:%.5f e x2:%.5f \na':%.3f \nalfa': %.2f \nOs escorregamentos específicos máximos:\nGs1B_max=%.5f e Gs2A_max=%.5f"%(x1,x2,a_linha,alfa_linha,Gs1B,Gs2A))
        
    else:
        Gs1B,Gs2A=escorregamentos(a,xSOL,-xSOL,alfa,alfa,z1,z2,0)
        print("\n\nSituação Z1+Z2>=60:")
        print("\nx1=%0.5f e x2=%.5f \n\nOs escorregamentos específicos máximos:\nGs1B_max=%.5f e Gs2A_max=%.5f" %(xSOL,-xSOL,Gs1B,Gs2A))
        
elif x=='Y' or x=='y':
    print('                               O K será considerado!')
    
    
    
    
    a=(z1+z2)*m/2
    def invol(alfa):
        alfa_rad=alfa*pi/180
        return tan(alfa_rad)-alfa_rad
    
    if z1+z2<60: z2_=60-z1
    else: z2_=z2
    def eq(x):
        a_=(z1+z2_)*m/2
        s1B,s2A=escorregamentos(a_,x,-x,alfa,alfa,z1,z2_,0)
        return s1B-s2A
    xSOL=fsolve(eq,0.3)
    
    if z1+z2<60:
        x1=xSOL
        def sist(z):
            alfa_rad=alfa*pi/180
            x2, a_linha, alfa_linha=z[0],z[1], z[2]
            a=(z1+z2)*m/2
            F0=a_linha*cos(alfa_linha*pi/180)-a*cos(alfa_rad)
            F1=invol(alfa)+2*tan(alfa_rad)*(x1+x2)/(z1+z2)-invol(alfa_linha)
            a1=(z1+z2)*m/2+(x1+x2)*m
            K=(a1-a_linha)/m
            Gs1B_linha,Gs2A_linha=escorregamentos(a_linha,x1,x2,alfa,alfa_linha,z1,z2,K)
            F2=Gs1B_linha-Gs2A_linha
            return [F0,F1,F2]
        x2,a_linha,alfa_linha=fsolve(sist,[0,a,alfa])
        a1=(z1+z2)*m/2+(x1+x2)*m
        K=(a1-a_linha)/m
        Gs1B,Gs2A=escorregamentos(a_linha,x1,x2,alfa,alfa_linha,z1,z2,K)
        print("\n\nSituação Z1+Z2<60:\nExige novo valor de entre-eixo(a') e ângulo de pressão(alfa');")
        print("x1:%.5f e x2:%.5f \na':%.3f \nalfa': %.2f \nOs escorregamentos específicos máximos:\nGs1B_max=%.5f e Gs2A_max=%.5f"%(x1,x2,a_linha,alfa_linha,Gs1B,Gs2A))
        
    else:
        Gs1B,Gs2A=escorregamentos(a,xSOL,-xSOL,alfa,alfa,z1,z2,0)
        print("\n\nSituação Z1+Z2>=60:")
        print("\nx1=%0.5f e x2=%0.5f \n\nOs escorregamentos específicos máximos:\nGs1B_max=%.5f e Gs2A_max=%.5f" %(xSOL,-xSOL,Gs1B,Gs2A))
        
else:
    print('***ENTRADA NÃO VÁLIDA***')


    #No caso de nos encontrarmos na gama de valores de z1+z2 pequeno, torna-se interessante perceber
    #que implicação teria não se ter aplicado o procedimento de Henriot, uma vez que ao projetista 
    #a imposição de um novo entre-eixo pode não ser conveniente 
    
    
    #De forma a que se possa tomar uma decisão entre considerar ou não este método,
    #apresentam-se as implicações em termos de aumento do escorregamento específico máximo.
if z1+z2<60: 
    y=input('Pretende conhecer os resultados caso não se aplique o sistema de Henriot?\nSe sim, digite y     ')
    if y=='y' or y=='Y':
        def eq(x):
            a=(z1+z2)*m/2
            s1B,s2A=escorregamentos(a,x,-x,alfa,alfa,z1,z2,0)
            return s1B-s2A
        xsol=fsolve(eq,0.3)
        Gs1B_,Gs2A_=escorregamentos(a,xsol,-xsol,alfa,alfa,z1,z2,0)
        print('\nNeste caso, os escorregamentos específicos máximos:\nGs1B_max=%.5f e Gs2A_max=%.5f'%(Gs1B_,Gs2A_))
        B=Gs1B_/Gs1B*100
        A=Gs2A_/Gs2A*100
        print('Gs1B_max e Gs2A_max serão %d%% dos valores calculados usando o sistema de Henriot' %(B))
        input('\n\n\n***ENTER para sair***')
else: 
    input('\n\n\n***ENTER para sair***')

        
    
