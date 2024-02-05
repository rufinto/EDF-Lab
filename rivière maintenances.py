# Simulateur pour le cas en rivière avec maintenances et sans biocides



###################################################################################################
###########################################   modules   ###########################################
###################################################################################################



import numpy as np
import matplotlib.pyplot as plt



######################################################################################################
###########################################   constantes   ###########################################
######################################################################################################



#géométrie:
L = 5 #longueur totale du tuyau en m
ri = 10e-3 #rayon interne du tube métallique en m
di= 2*ri
S=np.pi*ri**2
re=1.5e-2 #rayon exterieur du tube métallique en m
Se=2*np.pi*re*L #surface d'échange d'un tube
temps_total_simulation = 365*24*3600 #de 60 jours
dz = 3e-2 #le pas d'espace en m
dt = 3600 #pas de temps en s
temps=np.arange(0,temps_total_simulation,dt)
Nt=len(temps)
Z=np.arange(0,L,dz)
Nz=len(Z)

#constantes physiques générales fixes déterminées que je devrais plus avoir à bidouiller:
Rpropre = 5e-4 #résistance thermique du tube métallique m^2*K/W
Rinf = 0.346e-3  #resistance thermique final du biofilm obtenu par les données de Nabo
E = 25e3 #énergie d'activation J/mol
R = 8.314 #constante des gaz parfaits en J/mol/K
Lambda_bio = 0.6 #conductivité thermique du biofilm W/m/K
rho = 1e3 #masse volumique de l'eau kg/m^3
cp = 4184 #capacité calorifique massique de l'eau J/Kg/K
e0 = 1e-6 #valeur initiale de l'épaisseur du biofilm de Nabo en m
einf = Rinf*Lambda_bio #valeur de l'épaisseur finale obtenue par Nabo
k25=0.00967
eta_eau=1e-3
eta_pompe=0.7

#paramètres du cas:
Text = 60 + 273 #température de la fine couche d'eau sur la surface du tuyau en Kelvin (K)
Tin = 25 + 273 #température en entrée de l'eau dans le tuyau en Kelvin (K)
Dm = 0.8 #0.4 à 0.8    #dénit massique de l'écoulement en Kg/s
Dv=Dm/rho #débit volumique
v=Dv/S #vitesse moyenne de l'écoulement
Re=rho*2*ri*v/eta_eau #nombre de Reynolds

#constantes financière:
LAMBDA = 0.3164*(Re)**(-0.25) #coefficient de magie
Delta_P=LAMBDA*(L/(2*di))*rho*v**2 #perte de charge dans un tube
P_pump=Delta_P*Dv/eta_pompe #puissance pompe à un tube
fa=0.094



######################################################################################################
###########################################   simulateur   ###########################################
######################################################################################################



def oiler_maintenance(n):
    #description physique du système:
    def de(T,e):
        result=(k25/Lambda_bio) * np.exp( - (E/R) * ( 1/T - 1/298.5 ) )*(einf-e)*e
        return(result)

    def dT(T,e):
        result=( 1/(Dm*cp) ) * ( Text - T ) * 2*np.pi*ri / ( (e/Lambda_bio) + Rpropre )
        return(result)

    temps=np.arange(0,temps_total_simulation,dt)
    Nt=len(temps)

    Z=np.arange(0,L,dz)
    Nz=len(Z)

    e=np.zeros((Nt,Nz))
    T=np.zeros((Nt,Nz))

    T[::,0]=Tin
    e[0,::]=e0

    Pdis=np.zeros(Nt)

    #description économique du système:

    copex=0
    capex=20234*(n*Se)**0.61



    #résolution oiler avec maintenance:

    maintenace=24

    for z in range(1,Nz):
        T[0,z]=T[0,z-1]+dz*dT(T[0,z-1],e[0,z-1])
    for t in range(1,Nt):
        e[t,0]=e[t-1,0]+dt*de(T[t-1,0],e[t-1,0])


    Pdis[0]=n*Dm*cp*(T[0,-1]-T[0,0])


    for t in range(1,Nt):
        # on génere le profile dans le tube à t
        for z in range(1,Nz):

            T[t,z]=T[t,z-1]+dz*dT(T[t,z-1],e[t,z-1])
            e[t,z]=e[t-1,z]+dt*de(T[t-1,z],e[t-1,z])
        

        Pdis[t]=n*Dm*cp*(T[t,Nz-1]-T[t,0])


        # si à t la puissance dissipée par le radiateur ne suffit pas, une maintenance est déclanchée 
        if Pdis[t]<=20e6 and maintenace>0:
            Pdis[t]=0

            maintenace-=1
            copex+=( 0.35e-9 * 20e6/3600 )+( 75 * (n*Se)**0.6 )/24 #cout de la maintenance à chaques heures

        # quand la maintenace est finie on relance le radiateur
        if maintenace==0:
            e[t,::]=e0
            for t in range(t+1,Nt):
                e[t,0]=e[t-1,0]+dt*de(T[t-1,0],e[t-1,0])

            T[t,::]=np.copy(T[0,::])
            
            Pdis[t]=n*cp*Dm*(T[t,-1]-T[t,0])
            maintenace=24



        #quand le radiateur marche correctement il faut alimenter la pompe ça fait du copex
        if Pdis[t]>0:
            copex+=P_pump*n*130e-6
        

    return(T,e,Pdis,copex,capex)


def test_tubes(a,b):
    Tac=[]
    c=30
    for i in range(a,b+1,c):
        T,e,P,cop,cap=oiler_maintenance(i)
        Tac.append(cap*fa/40 + cop)
        print("computing : test_tubes  "+str((i-a)/(b-a)*100)+"%")  #il print le taux d'avancement du programme parce que c'est vraiment long

    plt.plot([i for i in range(a,b+1,c)],[Tac[i]/1000 for i in range(len(Tac))])
    plt.xlabel('nombre de tubes')
    plt.ylabel('TAC en k€')
    plt.title("évolution du TAC sur 40 ans en fonction du nombre de tubes")
    plt.legend()
    plt.show()



#####################################################################################################
###########################################   Résultats   ###########################################
#####################################################################################################


"""
test_tubes(900,1300)
#ce test montre que la maintenance est intéressante entre 400 et 500 tubes
"""
Ntubes=1100

T,e,P,copex,capex=oiler_maintenance(Ntubes)

P_couteux=np.array([max(0 , 20e6 - P[t])for t in range(Nt)])

plt.plot(temps/(3600*24),e[::,0]/1e-6)
plt.xlabel("temps(j)")
plt.ylabel("e(m e-6)")
plt.title("évolution de l'épaisseur de biofilme à l'entrée de l'échangeur en fonction du temps")
plt.legend()
plt.show()

plt.plot(temps/(3600*24),T[::,-1]-273)
plt.xlabel("temps(j)")
plt.ylabel("température(°C)")
plt.title("évolution de la température de sortie de l'échangeur en fonction du temps")
plt.legend()
plt.show()

plt.plot(temps/(3600*24),P/1000000,label="échangeur classique")
plt.plot(temps/(24*3600),P_couteux/1000000,label="échangeur additionnel")
plt.xlabel("temps(j)")
plt.ylabel("puissance(MW)")
plt.title("évolution de la puissance dissipée par l'échangeur, et de la puissance dissipée dans l'échangeur additionnel, en fonction du temps")
plt.legend()
plt.show()


TAC=capex*fa/40 + copex

print("copex = "+str(copex)+" €")
print("capex = "+str(capex)+" €")
print("TAC = "+str(TAC)+" €")

