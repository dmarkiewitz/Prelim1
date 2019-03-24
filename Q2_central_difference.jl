using LinearAlgebra

function Q2_central_difference(parameter_vector,sensative_parameter_index)
#preforming a varaiton of 1% of the parameter of intereset
#sensative_parameter_index is the parameter that we are looking at its sensativity's index

        dlength=2*0.01*parameter_vector[sensative_parameter_index]

        for index in 99:2:101
                global xk
                global x
                global i
                global x99
                global x101
                        parameter_vector[sensative_parameter_index]=parameter_vector[sensative_parameter_index]*index/100

                        #retriving the parameter
                        αmrna=parameter_vector[1]
                        αp=parameter_vector[2]
                        μcell=parameter_vector[3]
                        τx=parameter_vector[4]
                        τl=parameter_vector[5]
                        G=parameter_vector[6]
                        Lx=parameter_vector[7]
                        Lt=parameter_vector[8]
                        Rx=parameter_vector[9]
                        Rt=parameter_vector[10]
                        Lx1=parameter_vector[11]
                        Lx2=parameter_vector[12]
                        Lx3=parameter_vector[13]
                        LL1=parameter_vector[14]
                        LL2=parameter_vector[15]
                        LL3=parameter_vector[16]
                        vdotx=parameter_vector[17]
                        vdotl=parameter_vector[18]
                        Kx=parameter_vector[19]
                        Kl=parameter_vector[20]
                        Wi1=parameter_vector[21]
                        W11=parameter_vector[22]
                        W12=parameter_vector[23]
                        W13=parameter_vector[24]
                        W22=parameter_vector[25]
                        W23=parameter_vector[26]
                        W33=parameter_vector[27]
                        nm1=parameter_vector[28]
                        km1=parameter_vector[29]
                        nm2=parameter_vector[30]
                        km2=parameter_vector[31]
                        nm31=parameter_vector[32]
                        nm32=parameter_vector[33]
                        km31=parameter_vector[34]
                        km32=parameter_vector[35]

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
                                if k==0

                                        x=xk #initializing concentration at each step vecort
                                else

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
                                if k==0

                                        x=xk #initializing concentration at each step vecort

                                else

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

                                if k==0

                                        x=xk #initializing concentration at each step vecort

                                else

                                        x=hcat(x,xk) #building concentration at each step vecort
                                end
                                k=k+1
                                xk=Acarrot*xk+Scarrot*r(i,xk[4],xk[5])
                        end

                        if index==99
                                x99=x
                        else
                                x101=x
                        end
        end

        dxdp=(x101-x99)/dlength
        return dxdp
end
