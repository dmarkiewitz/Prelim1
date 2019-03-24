using LinearAlgebra

include("Q2_central_difference.jl")

#parameter assignments
αmrna=0.347 #inverse min units
αp=4.81*10^-4 #inverse min units
μcell=0.0173 #inverse min units
τx=2.7 #time constant for transcription
τl=0.8 #time constant for translation
G=0.356 #nmol/gDW
Lx=1000 #nt
Lt=333 #AA
Rx=2.046 #nmol/gDW
Rt=80.06 #nmol/gDW
Lx1=1200 #nt
Lx2=2400 #nt
Lx3=600 #nt
LL1=400 #AA
LL2=800 #AA
LL3=200 #AA
vdotx=3600 #nt/min RNAP elongation rate
vdotl=990 #AA/min Ribosome elongation rate
Kx=0.24 #nmol/gDW mrna saturation constant
Kl=454.64 #nmol/gDW protien saturation constant
Wi1=100
W11=10^-6
W12=10.0
W13=5.0
W22=10^-6
W23=50.0
W33=10^-6
nm1=1.5
km1=0.30 #in nmol/gDW
nm2=1.5
km2=1 #in nmol/gDW
nm31=1.5
nm32=10
km31=1 #in nmol/gDW
km32=10 #in nmol/gDW

#building parameter reference vector
Parameter_reference_vector=["αmrna","αp","μcell","τx","τl","G","Lx","Lt","Rx","Rt","Lx1","Lx2","Lx3","LL1","LL2","LL3","vdotx","vdotl","Kx","Kl","Wi1","W11","W12","W13","W22","W23","W33","nm1","km1","nm2","km2","nm31","nm32","km31","km32"]

#building parameter vector
Parameter_vector=[αmrna,αp,μcell,τx,τl,G,Lx,Lt,Rx,Rt,Lx1,Lx2,Lx3,LL1,LL2,LL3,vdotx,vdotl,Kx,Kl,Wi1,W11,W12,W13,W22,W23,W33,nm1,km1,nm2,km2,nm31,nm32,km31,km32]

#prebuilding arrays and vectors
A=zeros(6,6)#prebuilding dilution matrix
S=zeros(6,6)#prebuilding stochiometric matrix
Im=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;]#identity matrix

#regulation functions
#f are bound fraction functions
#u are regulatory functions
fm1(i)=(i^nm1)/((km1^nm1)+i^nm1)
Um1(i)=(W11+Wi1*fm1(i))/(1+W11+Wi1*fm1(i))

fm2(P1)=(P1^nm2)/((km2^nm2)+P1^nm2)
Um2(P1)=(W22+W12*fm2(P1))/(1+W22+W12*fm2(P1))

fm31(P1)=(P1^nm31)/((km31^nm31)+P1^nm31)
fm32(P2)=(P2^nm32)/((km32^nm32)+P2^nm32)
Um3(P1,P2)=(W33+W13*fm31(P1))/(1+W33+W23*fm32(P2)+W13*fm31(P1))

#misc "functions" i.e constants that decompose into more fundamental parameters
Ke1=vdotx/Lx1
Ke2=vdotx/Lx2
Ke3=vdotx/Lx3
Kex1=vdotl/LL1
Kex2=vdotl/LL2
Kex3=vdotl/LL3
λmrna=αmrna+μcell #inverse min units
λp=αp+μcell #inverse min units

#kinetic limit if transcription
vx1=Ke1*Rx*(G/(τx*Kx+(τx+1)*G))
vx2=Ke2*Rx*(G/(τx*Kx+(τx+1)*G))
vx3=Ke3*Rx*(G/(τx*Kx+(τx+1)*G))

#kinetic rate of transcription
rx1(i)=vx1*Um1(i)
rx2(P1)=vx2*Um2(P1)
rx3(P1,P2)=vx3*Um3(P1,P2)

#Other subfunctions
mRNAStar1(i)=rx1(i)/λmrna
mRNAStar2(P1)=rx2(P1)/λmrna
mRNAStar3(P1,P2)=rx3(P1,P2)/λmrna

#kinetic limit of translation
rl1(i)=Kex1*Rt*(mRNAStar1(i)/(τl*Kl+(τl+1)*mRNAStar1(i)))
rl2(P1)=Kex2*Rt*(mRNAStar2(P1)/(τl*Kl+(τl+1)*mRNAStar2(P1)))
rl3(P1,P2)=Kex3*Rt*(mRNAStar3(P1,P2)/(τl*Kl+(τl+1)*mRNAStar3(P1,P2)))

#setting time step size
τ=1.0 #min

#Building A matrix
A[1,1]=-λmrna
A[2,2]=-λmrna
A[3,3]=-λmrna
A[4,4]=-λp
A[5,5]=-λp
A[6,6]=-λp

#Building S matrix
S=Im

#Buidling Acarrot
Acarrot=exp(A*τ)

#Building Scarrot
Scarrot=inv(A)*(Acarrot-Im)*S

#Reaction vector
r(i,P1,P2)=[ rx1(i); rx2(P1); rx3(P1,P2); rl1(i); rl2(P1); rl3(P1,P2)]

#Setting intial concentration condition
xk=zeros(6,1)

#initializing steps
k=0

#number of steps
ns=1000

#initializing time memory vector
t=collect(0:1:ns)

#concentration of inducer
i=0

#running the ns+1 steps
while k<=length(t)-1
        global k
        if k==0
                global xk
                global x
                x=xk #initializing concentration at each step vecort

        else
                global xk
                global x
                x=hcat(x,xk) #building concentration at each step vecort
        end
        k=k+1
        xk=Acarrot*xk+Scarrot*r(i,xk[4],xk[5])
end

#adding 60 min and entering phase 1
ns=1000+60

#initializing time memory vector
t=collect(0:1:ns)

#running the ns+1 steps
while k<=length(t)-1
        global k
        if k==0
                global xk
                global x
                x=xk #initializing concentration at each step vecort

        else
                global xk
                global x
                x=hcat(x,xk) #building concentration at each step vecort
        end
        k=k+1
        xk=Acarrot*xk+Scarrot*r(i,xk[4],xk[5])
end

#adding inducer in mMol
i=10

#adding 300 min to run
ns=1060+300

#building time vector
t=collect(0:1:ns)

#running the ns+1 steps
while k<=length(t)-1
        global k
        global i
        if k==0
                global xk
                global x
                x=xk #initializing concentration at each step vecort

        else
                global xk
                global x
                x=hcat(x,xk) #building concentration at each step vecort
        end
        k=k+1
        xk=Acarrot*xk+Scarrot*r(i,xk[4],xk[5])
end

x
t

Parameter_reference_vector=["αmrna","αp","μcell","τx","τl","G","Lx","Lt","Rx","Rt","Lx1","Lx2","Lx3","LL1","LL2","LL3","vdotx","vdotl","Kx","Kl","Wi1","W11","W12","W13","W22","W23","W33","nm1","km1","nm2","km2","nm31","nm32","km31","km32"]

Parameter_vector=[αmrna,αp,μcell,τx,τl,G,Lx,Lt,Rx,Rt,Lx1,Lx2,Lx3,LL1,LL2,LL3,vdotx,vdotl,Kx,Kl,Wi1,W11,W12,W13,W22,W23,W33,nm1,km1,nm2,km2,nm31,nm32,km31,km32]

#building x reference vector
x_reference=["mRNA1","mRNA2","mRNA3","P1","P2","P3"]

#prebuilding sensativity matrix (states)
senp1=ones(length(x_reference),length(Parameter_reference_vector),21)
senp2E=ones(length(x_reference),length(Parameter_reference_vector),21)
senp2L=ones(length(x_reference),length(Parameter_reference_vector),21)

#building in scaling factors
for ti in 1020:1040
        for xi in 1:6
                for pi in 1:length(Parameter_reference_vector)
                        senp1[xi,pi,ti-1019]=(Parameter_vector[pi]/x[xi,ti])*senp1[xi,pi,ti-1019]
                end
        end
end

for ti in 1075:1095
        for xi in 1:6
                for pi in 1:length(Parameter_reference_vector)
                        senp2E[xi,pi,ti-1074]=(Parameter_vector[pi]/x[xi,ti])*senp2E[xi,pi,ti-1074]
                end
        end
end

for ti in 1325:1345
        for xi in 1:6
                for pi in 1:length(Parameter_reference_vector)
                        senp2L[xi,pi,ti-1324]=(Parameter_vector[pi]/x[xi,ti])*senp2L[xi,pi,ti-1324]
                end
        end
end

#Building in the sensativity information
for pi in 1:length(Parameter_reference_vector)
        dxdp=Q2_central_difference(Parameter_vector,pi)
        for xi in 1:6
                for ti in 1020:1040
                        senp1[xi,pi,ti-1019]=dxdp[xi,ti]*senp1[xi,pi,ti-1019]
                end
        end
end

for pi in 1:length(Parameter_reference_vector)
        dxdp=Q2_central_difference(Parameter_vector,pi)
        for xi in 1:6
                for ti in 1075:1095
                        senp2E[xi,pi,ti-1074]=dxdp[xi,ti]*senp2E[xi,pi,ti-1074]
                end
        end
end

for pi in 1:length(Parameter_reference_vector)
        dxdp=Q2_central_difference(Parameter_vector,pi)
        for xi in 1:6
                for ti in 1325:1345
                        senp2L[xi,pi,ti-1324]=dxdp[xi,ti]*senp2L[xi,pi,ti-1324]
                end
        end
end

senp1
senp2E
senp2L

#Code for Q.2.c. below

#building the absolute value matrixes
abssenp1=abs.(senp1)
abssenp2E=abs.(senp2E)
abssenp2L=abs.(senp2L)

#parameter assignment
T=20
Δx=1

#building the time averaged sensitivity array
P1Nij=(Δx/(2*T))*(abssenp1[:,:,1]+2*abssenp1[:,:,2]+2*abssenp1[:,:,3]+2*abssenp1[:,:,4]+2*abssenp1[:,:,5]+2*abssenp1[:,:,6]+2*abssenp1[:,:,7]+2*abssenp1[:,:,8]+2*abssenp1[:,:,9]+2*abssenp1[:,:,10]+2*abssenp1[:,:,11]+2*abssenp1[:,:,12]+2*abssenp1[:,:,13]+2*abssenp1[:,:,14]+2*abssenp1[:,:,15]+2*abssenp1[:,:,16]+2*abssenp1[:,:,17]+2*abssenp1[:,:,18]+2*abssenp1[:,:,19]+2*abssenp1[:,:,20]+abssenp1[:,:,21])
#Left empy here so the data in atom can go over this comment
P2ENij=(Δx/(2*T))*(abssenp2E[:,:,1]+2*abssenp2E[:,:,2]+2*abssenp2E[:,:,3]+2*abssenp2E[:,:,4]+2*abssenp2E[:,:,5]+2*abssenp2E[:,:,6]+2*abssenp2E[:,:,7]+2*abssenp2E[:,:,8]+2*abssenp2E[:,:,9]+2*abssenp2E[:,:,10]+2*abssenp2E[:,:,11]+2*abssenp2E[:,:,12]+2*abssenp2E[:,:,13]+2*abssenp2E[:,:,14]+2*abssenp2E[:,:,15]+2*abssenp2E[:,:,16]+2*abssenp2E[:,:,17]+2*abssenp2E[:,:,18]+2*abssenp2E[:,:,19]+2*abssenp2E[:,:,20]+abssenp2E[:,:,21])
#Left empy here so the data in atom can go over this comment
P2LNij=(Δx/(2*T))*(abssenp2L[:,:,1]+2*abssenp2L[:,:,2]+2*abssenp2L[:,:,3]+2*abssenp2L[:,:,4]+2*abssenp2L[:,:,5]+2*abssenp2L[:,:,6]+2*abssenp2L[:,:,7]+2*abssenp2L[:,:,8]+2*abssenp2L[:,:,9]+2*abssenp2L[:,:,10]+2*abssenp2L[:,:,11]+2*abssenp2L[:,:,12]+2*abssenp2L[:,:,13]+2*abssenp2L[:,:,14]+2*abssenp2L[:,:,15]+2*abssenp2L[:,:,16]+2*abssenp2L[:,:,17]+2*abssenp2L[:,:,18]+2*abssenp2L[:,:,19]+2*abssenp2L[:,:,20]+abssenp2L[:,:,21])
#Left empy here so the data in atom can go over this comment


#SVD's
P1=svd(P1Nij)
P2E=svd(P2ENij)
P2L=svd(P2LNij)

#Applying taking absolute value of first column vectors of U for each time
#window and ranking it by hand and doing the analysis and discussion on paper
saup1=abs.(P1.U[:,1])
saup2E=abs.(P2E.U[:,1])
saup2L=abs.(P2L.U[:,1])
