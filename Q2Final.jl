using Plots
using LinearAlgebra

#prebuilding arrays and vectors
A=zeros(6,6)#prebuilding dilution matrix
S=zeros(6,6)#prebuilding stochiometric matrix
Im=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;]#identity matrix

#parameter assignments
αmrna=0.347 #inverse min units
αp=4.81*10^-4 #inverse min units
μcell=0.0173 #inverse min units
λmrna=αmrna+μcell #inverse min units
λp=αp+μcell #inverse min units
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
Ke1=vdotx/Lx1
Ke2=vdotx/Lx2
Ke3=vdotx/Lx3
Kex1=vdotl/LL1
Kex2=vdotl/LL2
Kex3=vdotl/LL3

#regulation parameters
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

#building plot
Q2ai=plot(t,x[4,:],title="Model Ran to Steady State W/O Inducer",color="Blue",xlabel="Time (min)",ylabel="concentration (nmol/gDW)")
plot!(t,x[5,:],title="Model Ran to Steady State W/O Inducer",color="Red",xlabel="Time (min)",ylabel="concentration (nmol/gDW)")
plot!(t,x[6,:],title="Model Ran to Steady State W/O Inducer",color="Green",xlabel="Time (min)",ylabel="concentration (nmol/gDW)")

#saving the plot
savefig("!!!Location!!! \\Q2ai.png")

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

#buidling plot
Q2aii=plot(t[1000:1061],x[4,1000:1061],title="Phase 1",color="Blue",xlabel="Time (min)",ylabel="concentration (nmol/gDW)")
plot!(t[1000:1061],x[5,1000:1061],title="Phase 1",color="Red",xlabel="Time (min)",ylabel="concentration (nmol/gDW)")
plot!(t[1000:1061],x[6,1000:1061],title="Phase 1",color="Green",xlabel="Time (min)",ylabel="concentration (nmol/gDW)")

#saving the plot
savefig("!!!Location!!! \\Q2aii.png")

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

#building the plot
Q2aiii=plot(t[1060:1361],x[4,1060:1361],title="Phase 2",color="Blue",xlabel="Time (min)",ylabel="concentration (nmol/gDW)")
plot!(t[1060:1361],x[5,1060:1361],title="Phase 2",color="Red",xlabel="Time (min)",ylabel="concentration (nmol/gDW)")
plot!(t[1060:1361],x[6,1060:1361],title="Phase 2",color="Green",xlabel="Time (min)",ylabel="concentration (nmol/gDW)")

#saving the plot
savefig("!!!Location!!! \\Q2aiii.png")

#building total plot
Q2atotal=plot(t,x[4,:],title="I1-FFL Q2.a",color="Blue",xlabel="Time (min)",ylabel="concentration (nmol/gDW)",legend=:topleft)
plot!(t[:],x[5,:],title="I1-FFL Q2.a",color="Red",xlabel="Time (min)",ylabel="concentration (nmol/gDW)",legend=:topleft)
plot!(t[:],x[6,:],title="I1-FFL Q2.a",color="Green",xlabel="Time (min)",ylabel="concentration (nmol/gDW)",legend=:topleft)

#saving the plot
savefig("!!!Location!!! \\Q2atotal.png")
