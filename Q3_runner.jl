using JuMP
using Plots
using GLPKMathProgInterface
using DelimitedFiles

#for part a)

#building transcription rate and translation rate equations

#parameter assignments
αmrna=8.35 #inverse hr units
αp=9.9*10^-3 #inverse hr units
μcell=0 #inverse hr units
λmrna=αmrna+μcell #inverse hr units
λp=αp+μcell #inverse hr units
τx=2.7 #time constant for transcription
τl=0.8 #time constant for translation
G=5*10^-3 #μM
Lx=1000 #nt
Lt=330 #AA
Rx=0.15 #μM
Rl=1.6 #μM
Lx1=924 #nt
LL1=308 #AA
vdotx=216000 #nt/hr RNAP elongation rate
vdotl=59400 #AA/hr Ribosome elongation rate
Kx=0.3 #μM mrna saturation constant
Kl=57.0 #μM protien saturation constant
Ke1=vdotx/Lx1
Kex1=vdotl/LL1
nm=1.5
km=0.30 #mM
W1=0.26
W2=300.0

#regulation functions
#f are bound fraction functions
#u are regulatory functions
fm(ii)=(ii^nm)/((km^nm)+ii^nm)
Um(ii)=(W1+W2*fm(ii))/(1+W1+W2*fm(ii))

#kinetic limit if transcription
vx1=Ke1*Rx*(G/(τx*Kx+(τx+1)*G))

#kinetic rate of transcription
rXcarrot(ii)=vx1*Um(ii)

#Other subfunctions
mRNAStar1(ii)=rXcarrot(ii)/λmrna

#kinetic limit of translation
rL(ii)=Kex1*Rl*(mRNAStar1(ii)/(τl*Kl+(τl+1)*mRNAStar1(ii)))

#building V_boundtries by hand
V_boundaries_vector(ii)=[0 Inf;rXcarrot(ii) rXcarrot(ii);0 αmrna; 0 Inf;0 rL(ii); 0 Inf;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5]


#constructing stochiometric matrix
S_matrix=[-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
                0 -924 0 0 0 0 0 1 0 0 0 0 0 0 0;0 1 -1 -1 1 0 0 0 0 0 0 0 0 0 0;0 1848 0 0 616 2 0 0 0 0 0 0 0 0 -1;
                0 0 924 0 0 0 0 0 0 -1 0 0 0 0 0;0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0;0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0;
                0 0 0 0 -308 1 0 0 0 0 0 0 0 0 0;0 0 0 0 -616 0 0 0 0 0 0 0 1 0 0;0 0 0 0 308 -1 0 0 0 0 0 0 0 0 0;
                0 0 0 0 616 0 0 0 0 0 0 0 0 -1 0;0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0;0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0;
                0 0 0 0 0 -1 0 0 0 0 1 0 0 0 0;0 0 0 0 0 1 0 0 0 0 0 -1 0 0 0]

#for part b)

function Max_Translation(ii)

    #constructing stochiometric matrix
    S_matrix=[-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
                    0 -924 0 0 0 0 0 1 0 0 0 0 0 0 0;0 1 -1 -1 1 0 0 0 0 0 0 0 0 0 0;0 1848 0 0 616 2 0 0 0 0 0 0 0 0 -1;
                    0 0 924 0 0 0 0 0 0 -1 0 0 0 0 0;0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0;0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0;
                    0 0 0 0 -308 1 0 0 0 0 0 0 0 0 0;0 0 0 0 -616 0 0 0 0 0 0 0 1 0 0;0 0 0 0 308 -1 0 0 0 0 0 0 0 0 0;
                    0 0 0 0 616 0 0 0 0 0 0 0 0 -1 0;0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0;0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0;
                    0 0 0 0 0 -1 0 0 0 0 1 0 0 0 0;0 0 0 0 0 1 0 0 0 0 0 -1 0 0 0]

    #checking out S_matrix
    S_matrix

    #gets length of v vector
    reactions=15

    m=Int64(length(S_matrix)/reactions)

    #prebuilding C transposed
    cT=zeros(1,reactions)

    #Maximize translation rate v5
    cT[5]=1

    #S_matrix
    S_matrix

    #building transcription rate and translation rate equations

    #parameter assignments
    αmrna=8.35 #inverse hr units
    αp=9.9*10^-3 #inverse hr units
    μcell=0 #inverse hr units
    λmrna=αmrna+μcell #inverse hr units
    λp=αp+μcell #inverse hr units
    τx=2.7 #time constant for transcription
    τl=0.8 #time constant for translation
    G=5*10^-3 #μM
    Lx=1000 #nt
    Lt=330 #AA
    Rx=0.15 #μM
    Rl=1.6 #μM
    Lx1=924 #nt
    LL1=308 #AA
    vdotx=216000 #nt/hr RNAP elongation rate
    vdotl=59400 #AA/hr Ribosome elongation rate
    Kx=0.3 #μM mrna saturation constant
    Kl=57.0 #μM protien saturation constant
    Ke1=vdotx/Lx1
    Kex1=vdotl/LL1
    nm=1.5
    km=0.30 #mM
    W1=0.26
    W2=300.0

    #regulation functions
    #f are bound fraction functions
    #u are regulatory functions
    fm(ii)=(ii^nm)/((km^nm)+ii^nm)
    Um(ii)=(W1+W2*fm(ii))/(1+W1+W2*fm(ii))

    #kinetic limit if transcription
    vx1=Ke1*Rx*(G/(τx*Kx+(τx+1)*G))

    #kinetic rate of transcription
    rXcarrot(ii)=vx1*Um(ii)

    #Other subfunctions
    mRNAStar1(ii)=rXcarrot(ii)/λmrna

    #kinetic limit of translation
    rL(ii)=Kex1*Rl*(mRNAStar1(ii)/(τl*Kl+(τl+1)*mRNAStar1(ii)))

    #building V_boundtries by hand
    V_boundaries_vector=[0 Inf;rXcarrot(ii) rXcarrot(ii);0 αmrna; 0 Inf;0 rL(ii); 0 Inf;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5]

    #check up on information vectors
    cT

    #FBA constraint
    zero_vector_FBA_constraint=zeros(m,1)

    #Model Construction

    P_Model=Model(solver=GLPKSolverLP())

    #Variables

    @variable(P_Model, V_boundaries_vector[i,1] <= V_vector[i in 1:reactions] <= V_boundaries_vector[i,2])

    #Objective

    @objective(P_Model, Max, (cT*V_vector)[1])

    #Constraints

    @constraint(P_Model, S_matrix*V_vector .== zero_vector_FBA_constraint)

    #Output
    Status=solve(P_Model)

    print(P_Model)

    v=zeros(reactions,1)

    #Building solution vector
    for i in 1:reactions
        v[i]=getvalue(V_vector[i])
    end

    #here is the flux's for this system
    v

    #Obtaining optimized value, here it is the maximuim translation rate
    opt_value=(cT*v)[1]

    return opt_value/λp #to get steady state level need to divide by the total degradation rate of the protien
end

#establishing the under and lower bounds on the inducer concentration in μM
I1=0.0001
I2=10

#creating plot
plot(Max_Translation,I1,I2, xscale = :log10, xlabel="Inducer Concentration (mM)", ylabel="Steady-State Protein Conc. (microM)",title="Steady-State Protein Conc. vs Inducer Conc.")

#Saves plot
savefig("!!!Location!!!\\Q3.b.png")

#maximium v5 for any inducer concentration in mM call the following function
function Max_Translation_rate(ii)

    #constructing stochiometric matrix
    S_matrix=[-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
                    0 -924 0 0 0 0 0 1 0 0 0 0 0 0 0;0 1 -1 -1 1 0 0 0 0 0 0 0 0 0 0;0 1848 0 0 616 2 0 0 0 0 0 0 0 0 -1;
                    0 0 924 0 0 0 0 0 0 -1 0 0 0 0 0;0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0;0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0;
                    0 0 0 0 -308 1 0 0 0 0 0 0 0 0 0;0 0 0 0 -616 0 0 0 0 0 0 0 1 0 0;0 0 0 0 308 -1 0 0 0 0 0 0 0 0 0;
                    0 0 0 0 616 0 0 0 0 0 0 0 0 -1 0;0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0;0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0;
                    0 0 0 0 0 -1 0 0 0 0 1 0 0 0 0;0 0 0 0 0 1 0 0 0 0 0 -1 0 0 0]

    #checking out S_matrix
    S_matrix

    #gets length of v vector
    reactions=15

    m=Int64(length(S_matrix)/reactions)

    #prebuilding C transposed
    cT=zeros(1,reactions)

    #Maximize translation rate v5
    cT[5]=1

    #S_matrix
    S_matrix

    #building transcription rate and translation rate equations

    #parameter assignments
    αmrna=8.35 #inverse hr units
    αp=9.9*10^-3 #inverse hr units
    μcell=0 #inverse hr units
    λmrna=αmrna+μcell #inverse hr units
    λp=αp+μcell #inverse hr units
    τx=2.7 #time constant for transcription
    τl=0.8 #time constant for translation
    G=5*10^-3 #μM
    Lx=1000 #nt
    Lt=330 #AA
    Rx=0.15 #μM
    Rl=1.6 #μM
    Lx1=924 #nt
    LL1=308 #AA
    vdotx=216000 #nt/hr RNAP elongation rate
    vdotl=59400 #AA/hr Ribosome elongation rate
    Kx=0.3 #μM mrna saturation constant
    Kl=57.0 #μM protien saturation constant
    Ke1=vdotx/Lx1
    Kex1=vdotl/LL1
    nm=1.5
    km=0.30 #mM
    W1=0.26
    W2=300.0

    #regulation functions
    #f are bound fraction functions
    #u are regulatory functions
    fm(ii)=(ii^nm)/((km^nm)+ii^nm)
    Um(ii)=(W1+W2*fm(ii))/(1+W1+W2*fm(ii))

    #kinetic limit if transcription
    vx1=Ke1*Rx*(G/(τx*Kx+(τx+1)*G))

    #kinetic rate of transcription
    rXcarrot(ii)=vx1*Um(ii)

    #Other subfunctions
    mRNAStar1(ii)=rXcarrot(ii)/λmrna

    #kinetic limit of translation
    rL(ii)=Kex1*Rl*(mRNAStar1(ii)/(τl*Kl+(τl+1)*mRNAStar1(ii)))

    #building V_boundtries by hand
    V_boundaries_vector=[0 Inf;rXcarrot(ii) rXcarrot(ii);0 αmrna; 0 Inf;0 rL(ii); 0 Inf;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5;-10^5 10^5]

    #check up on information vectors
    cT

    #FBA constraint
    zero_vector_FBA_constraint=zeros(m,1)

    #Model Construction

    P_Model=Model(solver=GLPKSolverLP())

    #Variables

    @variable(P_Model, V_boundaries_vector[i,1] <= V_vector[i in 1:reactions] <= V_boundaries_vector[i,2])

    #Objective

    @objective(P_Model, Max, (cT*V_vector)[1])

    #Constraints

    @constraint(P_Model, S_matrix*V_vector .== zero_vector_FBA_constraint)

    #Output
    Status=solve(P_Model)

    print(P_Model)

    v=zeros(reactions,1)

    #Building solution vector
    for i in 1:reactions
        v[i]=getvalue(V_vector[i])
    end

    #here is the flux's for this system
    v

    #Obtaining optimized value, here it is the maximuim translation rate
    opt_value=(cT*v)[1]

    return opt_value
end
