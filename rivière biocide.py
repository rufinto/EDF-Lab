# Simulateur pour le cas en rivière avec utilisation de biocides et sans maintenances



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
di=2*ri
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

#constante physique générales:
Rpropre = 5e-4 #résistance thermique du tube métallique m^2*K/W
Rinf = 0.346e-3  #resistance thermique final du biofilm obtenu par les données de Nabo
E = 25e3 #énergie d'activation J/mol
R = 8.314 #constante des gaz parfaits en J/mol/K
Lambda_bio = 0.6 #conductivité thermique du biofilm W/m/K
rho = 1e3 #masse volumique de l'eau kg/m^3
cp = 4184 #capacité calorifique massique de l'eau J/Kg/K
e0 = 1e-6 #valeur initiale de l'épaisseur du biofilm de Nabo en m
einf = Rinf*Lambda_bio #valeur de l'épaisseur finale obtenue par Nabo
k25=0.00967 #coeff cinétique loi de Nebolt
eta_eau=1e-3
eta_pompe=0.7

#paramètres du cas:
Text = 60 + 273 #température de la fine couche d'eau sur la surface du tuyau en Kelvin (K)
Tin = 25 + 273 #température en entrée de l'eau dans le tuyau en Kelvin (K)
Dm = 0.8 #0.4 à 0.8    #dénit massique de l'écoulement en Kg/s
Dv=Dm/rho #débit volumique
v=Dv/S #vitesse moyenne de l'écoulement
Re=rho*2*ri*v/eta_eau #nombre de Reynolds
m_eau=Dm*365*24*3600 #masse d'eau à traiter pendant un an
m_biocide=m_eau*5e-3 #masse de biocide pour traiter l'eau

#avec biocide la valeur finale de e est modifiée
einf=einf*0.4872 

#constantes financière:
LAMBDA = 0.3164*(Re)**(-0.25) #coefficient de magie
Delta_P=LAMBDA*(L/(2*di))*rho*v**2 #perte de charge dans un tube
P_pump=Delta_P*Dv/eta_pompe #puissance pompe à un tube
fa=0.094



######################################################################################################
###########################################   simulateur   ###########################################
######################################################################################################



def oiler():
    
    def de(T,e):
        result=(k25/Lambda_bio) * np.exp( - (E/R) * ( 1/T - 1/298.5 ) )*(einf-e)*e
        return(result)

    def dT(T,e):
        result=( 1/(Dm*cp) ) * ( Text - T ) * 2*np.pi*ri / ( (e/Lambda_bio) + Rpropre )
        return(result)

    e=np.zeros((Nt,Nz))
    T=np.zeros((Nt,Nz))

    T[::,0]=Tin
    e[0,::]=e0


    #résolution oiler:

    for z in range(1,Nz):
        T[0,z]=T[0,z-1]+dz*dT(T[0,z-1],e[0,z-1])
    for t in range(1,Nt):
        e[t,0]=e[t-1,0]+dt*de(T[t-1,0],e[t-1,0])

    for t in range(1,Nt):
        for z in range(1,Nz):
            T[t,z]=T[t,z-1]+dz*dT(T[t,z-1],e[t,z-1])
            e[t,z]=e[t-1,z]+dt*de(T[t-1,z],e[t-1,z])

    return(T,e)



#####################################################################################################
###########################################   Résultats   ###########################################
#####################################################################################################

n=1310

T,e=oiler()

P_dis=n*Dm*cp*(T[::,-1]-T[::,0])

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

plt.plot(temps/(3600*24),P_dis/1000000)
plt.xlabel("temps(j)")
plt.ylabel("puissance(MW)")
plt.title("évolution de la puissance dissipée par l'échangeur en fonction du temps")
plt.legend()
plt.show()




copex=P_pump*n*130e-6 + m_biocide*1
capex=20234*(n*Se)**0.61
TAC=capex*fa/40 + copex


print("copex = "+str(copex)+" €")
print("capex = "+str(capex)+" €")
print("TAC = "+str(TAC)+" €")