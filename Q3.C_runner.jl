using JuMP
using Plots
using GLPKMathProgInterface
using DelimitedFiles
using LinearAlgebra

#vector contraining rates for the optimized objective function z which is equal to 0.17075390602405563
ov=[0.21091091457007352,0.21091091457007352,0.21091091457007352,0.17075390602405563,0.17075390602405563,
    52.592203055409136,52.59220305540657,194.88168506274815,0.17075390602985863,194.88168506274815,
    52.59220305540657,52.59220305540657,105.18440611081314,105.18440611081314,600.1321823471226]

#assigning the maximimuim optimized value for v5 determined previously in Q3.b / using Q3_runnner
#and using function Max_Translation(ii) with ii equal to 10
max_opt_value=0.17075390602405563

function shadow_price(Ex_Flux_to_Test_Index,Metabolite_associtated_to_flux,max_opt_value,ov,perturbation)
#perturbations should be less than 1
#perturbation should be constant across the trials as to eliminate the need for
#scaling factor

        #for FBA zero constraints
        a=zeros(1,17)
        b=zeros(1,17)

        #for shadow prices / perturbation constraints
        c=ones(1,15)
        d=zeros(1,15)

        a[Metabolite_associtated_to_flux]=-4
        b[Metabolite_associtated_to_flux]=4
        c[Ex_Flux_to_Test_Index]=0
        d[Ex_Flux_to_Test_Index]=ov[Ex_Flux_to_Test_Index]*perturbation

        #assinging indcuer concentration
        ii=10 #mM

        #constructing stochiometric matrix
        S_matrix=[-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;-1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
                        0 -924 0 0 0 0 0 1 0 0 0 0 0 0 0;0 1 -1 -1 1 0 0 0 0 0 0 0 0 0 0;0 1848 0 0 616 2 0 0 0 0 0 0 0 0 -1;
                        0 0 924 0 0 0 0 0 0 -1 0 0 0 0 0;0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0;0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 -308 1 0 0 0 0 0 0 0 0 0;0 0 0 0 -616 0 0 0 0 0 0 0 1 0 0;0 0 0 0 308 -1 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 616 0 0 0 0 0 0 0 0 -1 0;0 0 0 0 1 0 0 0 -1 0 0 0 0 0 0;0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 -1 0 0 0 0 1 0 0 0 0;0 0 0 0 0 1 0 0 0 0 0 -1 0 0 0]

        #preconstructing rows of the S_matrix so that JuMP can process it correctly
        sm=[]
        #building vectors
        sm=push!(sm,transpose(S_matrix[1,1:15]))
        sm=push!(sm,transpose(S_matrix[2,1:15]))
        sm=push!(sm,transpose(S_matrix[3,1:15]))
        sm=push!(sm,transpose(S_matrix[4,1:15]))
        sm=push!(sm,transpose(S_matrix[5,1:15]))
        sm=push!(sm,transpose(S_matrix[6,1:15]))
        sm=push!(sm,transpose(S_matrix[7,1:15]))
        sm=push!(sm,transpose(S_matrix[8,1:15]))
        sm=push!(sm,transpose(S_matrix[9,1:15]))
        sm=push!(sm,transpose(S_matrix[10,1:15]))
        sm=push!(sm,transpose(S_matrix[11,1:15]))
        sm=push!(sm,transpose(S_matrix[12,1:15]))
        sm=push!(sm,transpose(S_matrix[13,1:15]))
        sm=push!(sm,transpose(S_matrix[14,1:15]))
        sm=push!(sm,transpose(S_matrix[15,1:15]))
        sm=push!(sm,transpose(S_matrix[16,1:15]))
        sm=push!(sm,transpose(S_matrix[17,1:15]))

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
        V_boundaries_vector=[0 Inf;rXcarrot(ii) rXcarrot(ii);0 αmrna; 0 Inf;0 rL(ii); 0 Inf;-10^5 d[7]+c[7]*10^5;-10^5 d[8]+c[8]*10^5;-10^5 d[9]+c[9]*10^5;-10^5 d[10]+c[10]*10^5;-10^5 d[11]+c[11]*10^5;-10^5 d[12]+c[12]*10^5;-10^5 d[13]+c[13]*10^5;-10^5 d[14]+c[14]*10^5;-10^5 d[15]+c[15]*10^5]

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

        @constraint(P_Model, a[1] <= (sm[1]*V_vector) <= b[1])
        @constraint(P_Model, a[2] <= (sm[2]*V_vector) <= b[2])
        @constraint(P_Model, a[3] <= (sm[3]*V_vector) <= b[3])
        @constraint(P_Model, a[4] <= (sm[4]*V_vector) <= b[4])
        @constraint(P_Model, a[5] <= (sm[5]*V_vector) <= b[5])
        @constraint(P_Model, a[6] <= (sm[6]*V_vector) <= b[6])
        @constraint(P_Model, a[7] <= (sm[7]*V_vector) <= b[7])
        @constraint(P_Model, a[8] <= (sm[8]*V_vector) <= b[8])
        @constraint(P_Model, a[9] <= (sm[9]*V_vector) <= b[9])
        @constraint(P_Model, a[10] <= (sm[10]*V_vector) <= b[10])
        @constraint(P_Model, a[11] <= (sm[11]*V_vector) <= b[11])
        @constraint(P_Model, a[12] <= (sm[12]*V_vector) <= b[12])
        @constraint(P_Model, a[13] <= (sm[13]*V_vector) <= b[13])
        @constraint(P_Model, a[14] <= (sm[14]*V_vector) <= b[14])
        @constraint(P_Model, a[15] <= (sm[15]*V_vector) <= b[17])
        @constraint(P_Model, a[16] <= (sm[16]*V_vector) <= b[16])
        @constraint(P_Model, a[17] <= (sm[17]*V_vector) <= b[17])

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

    return abs((opt_value-max_opt_value)/((1-perturbation)*ov[Ex_Flux_to_Test_Index])) #to calculate the shadow price
end

#prebuilding shadow storage vector
λ=zeros(1,9)

#obtaining the shadow prices for the exchange fluxes

λ[1]=shadow_price(7,15,max_opt_value,ov,0.99)
λ[2]=shadow_price(8,4,max_opt_value,ov,0.99)
λ[3]=shadow_price(9,14,max_opt_value,ov,0.99)
λ[4]=shadow_price(10,7,max_opt_value,ov,0.99)
λ[5]=shadow_price(11,16,max_opt_value,ov,0.99)
λ[6]=shadow_price(12,17,max_opt_value,ov,0.99)
λ[7]=shadow_price(13,11,max_opt_value,ov,0.99)
λ[8]=shadow_price(14,13,max_opt_value,ov,0.99)
λ[9]=shadow_price(15,6,max_opt_value,ov,0.99)

#displaying the shadow prices
println(λ)

#thus from the print shadow prices one can tell that Exchange flux 15 has the largest shadow price
#and thus the translation rate is most sensative to.
