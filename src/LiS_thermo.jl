### LiS module ###
### Linnette Teo 2022.01.07 ###
### Thermodynamic model to solve for equilibrium solutions for lithium sulfur system ###
### Three defined regions based on differing solids present ###


    
function thermoModel1!(F,x)
    """
    Equations for region 1: only S8 solid
    """

    global dod
    # variables to solve
    CLi = exp(x[1])
    CS8l = exp(x[2])
    CS8m = exp(x[3])
    CS6m = exp(x[4])
    CS4m = exp(x[5])
    CS2m = exp(x[6])
    CSm2 = exp(x[7])
    phi1 = exp(x[8])
    epcS8 = exp(x[9])


    # porosites
#     epcS8 = 0
    epcLi2S = 0
    epsS8 = 0
    epsLi2S = 0
    ec = 1-epcS8-epcLi2S-efc
    es = 1-epsS8-epsLi2S-efs

    # potential equations, setting overpotentials=0
    F[1] = phi1 - U2 + f *(sS8m_2*log(CS8m/1000) + sS8l_2*log(CS8l/1000) )
    F[2] = phi1 - U3 + f *(sS6m_3*log(CS6m/1000) + sS8m_3*log(CS8m/1000) )
    F[3] = phi1 - U4 + f *(sS4m_4*log(CS4m/1000) + sS6m_4*log(CS6m/1000) )
    F[4] = phi1 - U5 + f *(sS2m_5*log(CS2m/1000) + sS4m_5*log(CS4m/1000) )
    F[5] = phi1 - U6 + f *(sSm2_6*log(CSm2/1000) + sS2m_6*log(CS2m/1000) )

    # charge balance
    F[6] = (ec*Lcat + es*Lsep) * 2*(CS8m + CS6m + CS4m + CS2m + CSm2) +
            epcLi2S/VLi2S*Lcat*2 + epsLi2S/VLi2S*Lsep*2 - S_bal*2*dod/100 - inichg

    # electroneutrality
    F[7] = CAref*(1-efc-e_S8init)/ec + 2*(CS8m + CS6m + CS4m + CS2m + CSm2) - CLi

    # S mass balance
    F[8] = (ec*Lcat + es*Lsep) * (CS8l*8 + CS8m*8 + CS6m*6 + CS4m*4 + CS2m*2 + CSm2*1) +
            epcS8/VS8*Lcat*8 + epsS8/VS8*Lsep*8 + epcLi2S/VLi2S*Lcat*1 + epsLi2S/VLi2S*Lsep*1 -S_bal

    # Precipitation
    F[9] = CS8l - KspS8

end


function thermoModel2!(F,x)
    """
    Equations for region 2: no solids
    """

    global dod
    # variables to solve
    CLi = exp(x[1])
    CS8l = exp(x[2])
    CS8m = exp(x[3])
    CS6m = exp(x[4])
    CS4m = exp(x[5])
    CS2m = exp(x[6])
    CSm2 = exp(x[7])
    phi1 = exp(x[8])

    # porosites
    epcS8 = 0
    epcLi2S = 0
    epsS8 = 0
    epsLi2S = 0
    ec = 1-epcS8-epcLi2S-efc
    es = 1-epsS8-epsLi2S-efs

    # potential equations, setting overpotentials=0
    F[1] = phi1 - U2 + f *(sS8m_2*log(CS8m/1000) + sS8l_2*log(CS8l/1000) )
    F[2] = phi1 - U3 + f *(sS6m_3*log(CS6m/1000) + sS8m_3*log(CS8m/1000) )
    F[3] = phi1 - U4 + f *(sS4m_4*log(CS4m/1000) + sS6m_4*log(CS6m/1000) )
    F[4] = phi1 - U5 + f *(sS2m_5*log(CS2m/1000) + sS4m_5*log(CS4m/1000) )
    F[5] = phi1 - U6 + f *(sSm2_6*log(CSm2/1000) + sS2m_6*log(CS2m/1000) )

    # charge balance
    F[6] = (ec*Lcat + es*Lsep) * 2*(CS8m + CS6m + CS4m + CS2m + CSm2) +
            epcLi2S/VLi2S*Lcat*2 + epsLi2S/VLi2S*Lsep*2 - S_bal*2*dod/100 - inichg

    # electroneutrality
    F[7] = CAref*(1-efc-e_S8init)/ec + 2*(CS8m + CS6m + CS4m + CS2m + CSm2) - CLi

    # S mass balance
    F[8] = (ec*Lcat + es*Lsep) * (CS8l*8 + CS8m*8 + CS6m*6 + CS4m*4 + CS2m*2 + CSm2*1) +
            epcS8/VS8*Lcat*8 + epsS8/VS8*Lsep*8 + epcLi2S/VLi2S*Lcat*1 + epsLi2S/VLi2S*Lsep*1 - S_bal

end



function thermoModel3!(F,x)
    """
    Equations for region 3: only Li2S solid
    """

    global dod
    # variables to solve
    CLi = exp(x[1])
    CS8l = exp(x[2])
    CS8m = exp(x[3])
    CS6m = exp(x[4])
    CS4m = exp(x[5])
    CS2m = exp(x[6])
    CSm2 = exp(x[7])
    phi1 = exp(x[8])
    epcLi2S = exp(x[9])

    # porosites
    epcS8 = 0
#     epcLi2S = 0
    epsS8 = 0
    epsLi2S = 0
    ec = 1-epcS8-epcLi2S-efc
    es = 1-epsS8-epsLi2S-efs

    # potential equations, setting overpotentials=0
    F[1] = phi1 - U2 + f *(sS8m_2*log(CS8m/1000) + sS8l_2*log(CS8l/1000) )
    F[2] = phi1 - U3 + f *(sS6m_3*log(CS6m/1000) + sS8m_3*log(CS8m/1000) )
    F[3] = phi1 - U4 + f *(sS4m_4*log(CS4m/1000) + sS6m_4*log(CS6m/1000) )
    F[4] = phi1 - U5 + f *(sS2m_5*log(CS2m/1000) + sS4m_5*log(CS4m/1000) )
    F[5] = phi1 - U6 + f *(sSm2_6*log(CSm2/1000) + sS2m_6*log(CS2m/1000) )

    # charge balance
    F[6] = (ec*Lcat + es*Lsep) * 2*(CS8m + CS6m + CS4m + CS2m + CSm2) +
            epcLi2S/VLi2S*Lcat*2 + epsLi2S/VLi2S*Lsep*2 - S_bal*2*dod/100 - inichg

    # electroneutrality
    F[7] = CAref*(1-efc-e_S8init)/ec + 2*(CS8m + CS6m + CS4m + CS2m + CSm2) - CLi

    # S mass balance
    F[8] = (ec*Lcat + es*Lsep) * (CS8l*8 + CS8m*8 + CS6m*6 + CS4m*4 + CS2m*2 + CSm2*1) +
            epcS8/VS8*Lcat*8 + epsS8/VS8*Lsep*8 + epcLi2S/VLi2S*Lcat*1 + epsLi2S/VLi2S*Lsep*1 - S_bal

    # Precipitation
    F[9] = CLi^2 * CSm2 - KspLi2S

end


function findIdx(dod,step)
    """Finds index of dod in range 1:step:99"""
    if step < 1
        idx = ceil(Int, dod/step) -1
    else 
        idx = ceil(Int, dod/step)
    end
    
    return idx
end
     
    
function runThermoRegion!(dod_arr, x0, thermoModel, sol_mat, tol, step)
    """
    Solves thermoModel for a specific region over a range dod_arr.
    dod_arr has step size of Int step.
    Starts from initial guess x0.
    Store voltage solutions in place in V_arr.
    Store all solutions in place in sol_mat (initialized as 10 columns)
    Solver used is nlsolve, newton method.
    Breaks if solution does not converge OR
    if we hit a tolerance limit (i.e. vol fracs < tol)
    """

    for i in 1:length(dod_arr)
        global dod = dod_arr[i]

        sol = nlsolve(thermoModel, x0, autodiff = :forward, method=:newton,xtol=1e-10,ftol=1e-10)

        # update initial guess to previous DOD solution
        global x0 = sol.zero

        if converged(sol)
            if ((thermoModel == thermoModel1!) && (exp.(x0[9])>tol)) ||
                thermoModel == thermoModel2! ||
                ((thermoModel == thermoModel3!) && (exp.(x0[9])>tol)) ||

                ((thermoModel == thermoModel1_S3!) && (exp.(x0[9])>tol)) ||
                    thermoModel == thermoModel2_S3! ||
                    ((thermoModel == thermoModel3_S3!) && (exp.(x0[9])>tol))

                # store solutions in sol_mat
                dod_ind  = findIdx(dod,step)   # index of dod
                store_sol!(thermoModel, sol_mat, sol, dod_ind)
#                 print(dod, ", ") #'\n')
            else
#                 print("\n hits tol=", exp.(x0[9]), "\n")
                break
            end
        else
#             print("\n not converged for dod=", dod, "\n")
            break
        end
    end
end            

                                         
function store_sol!(thermoModel, sol_mat, sol, dod_ind)
    """
    Stores solutions into sol_mat depending on which model is used.
    Stored at the index Int dod_ind.
    For indexing purposes as different models have variables in different indices.
    Variable order: CLi, CS8l, CS8m, CS6m, CS4m, CS2m, CSm2, phi1, epcS8, epcLi2S, CS3
    """
    ## models without S3 concentration as a variable
    if size(sol_mat)[2]==10
        if thermoModel == thermoModel1!
            sol_mat[dod_ind,:] = append!(exp.(sol.zero),0)
        elseif thermoModel == thermoModel2!
            sol_mat[dod_ind,:] = append!(exp.(sol.zero),[0,0])
        elseif thermoModel == thermoModel3!
            sol_mat[dod_ind,:] = insert!(exp.(sol.zero),9,0)
        else
            print("\n solution array is > Nvar (number of variables)")
        end
    end

    ## with S3 radical concentration in index 11
    if size(sol_mat)[2]==11
        if thermoModel == thermoModel1_S3!
            sol_mat[dod_ind,:] = insert!(exp.(sol.zero),10,0)
        elseif thermoModel == thermoModel2_S3!
            sol_mat[dod_ind,:] = insert!(insert!(exp.(sol.zero),9,0),9,0)
        elseif thermoModel == thermoModel3_S3!
            sol_mat[dod_ind,:] = insert!(exp.(sol.zero),9,0)
        else
            print("\n solution array is > Nvar (number of variables)")
        end
    end


end                                            
                                                                                                
                                                
function solveThermoAll!(sol_mat, x0_99, x0_1; tol=1e-3, step=1)

    """
    Algorithm:
    Solve Region 3 in descending order, starting from 99% DOD till no more solution found.
    Solve Region 1 in ascending order, starting from 1% DOD till no more solution found.                   
    Solve Region 2 in ascending order, starting from where Region 1 ended till where Region 3 begun.

    Stores all solutions in sol_mat, a len(dod_arr) x Nvar=10 matrix
    x0_99: initial guess for 99% dod, 
    x0_1: initial guess for 1% dod
                                                
    also return dod_r2LB and dod_r2UB which are the lower and upper bounds of region 2, inclusive 
                                            
    Default arguments:
    tol: tolerance limit for solver for vol fracs (see runThermoRegion!)                             step: step size for dod_arr              
                                                
    """
    #### RUNNING REGION 3 #####
    #print("\nREGION 3\n")
    global x0 = x0_99
    dod_arr = 99:-step:1
    thermoModel = thermoModel3!
    runThermoRegion!(dod_arr, x0, thermoModel, sol_mat, tol, step)
    dod_r2UB = dod   # dod where region 3 begins and region 2 ends


    #### RUNNING REGION 1 #####
    #print("\nREGION 1\n")
    global x0 = x0_1
    dod_arr = 1:step:99
    thermoModel = thermoModel1!
    runThermoRegion!(dod_arr, x0, thermoModel, sol_mat, tol, step)
    dod_r2LB = dod


    ### RUNNING REGION 2 #####
    #print("\nREGION 2\n")
    dod_arr = dod_r2LB:step:dod_r2UB
    if dod_r2LB==1
        global x0 = x0_1[1:8]
    else         # look at element before  the LB for initial guess
        idx  = findIdx(dod,step)                                   
        global x0 = log.(sol_mat[idx-1,1:8])                                        
    end
    thermoModel = thermoModel2!
    runThermoRegion!(dod_arr, x0, thermoModel, sol_mat, tol, step)


    return dod_r2LB, dod_r2UB
end






######################### S3 RADICAL MODEL ############################

function thermoModel1_S3!(F,x)
    """
    Equations for region 1: only S8 solid
    With S3 radical chemistry
    """

    global dod
    # variables to solve
    CLi = exp(x[1])
    CS8l = exp(x[2])
    CS8m = exp(x[3])
    CS6m = exp(x[4])
    CS4m = exp(x[5])
    CS2m = exp(x[6])
    CSm2 = exp(x[7])
    phi1 = exp(x[8])
    epcS8 = exp(x[9])
    CS3 = exp(x[10])


    # porosites
#     epcS8 = 0
    epcLi2S = 0
    epsS8 = 0
    epsLi2S = 0
    ec = 1-epcS8-epcLi2S-efc
    es = 1-epsS8-epsLi2S-efs

    # potential equations, setting overpotentials=0
    F[1] = phi1 - U2 + f *(sS8m_2*log(CS8m/1000) + sS8l_2*log(CS8l/1000) )
    F[2] = phi1 - U3 + f *(sS6m_3*log(CS6m/1000) + sS8m_3*log(CS8m/1000) )
    F[3] = phi1 - U4 + f *(sS4m_4*log(CS4m/1000) + sS6m_4*log(CS6m/1000) )
    F[4] = phi1 - U5 + f *(sS2m_5*log(CS2m/1000) + sS4m_5*log(CS4m/1000) )
    F[5] = phi1 - U6 + f *(sSm2_6*log(CSm2/1000) + sS2m_6*log(CS2m/1000) )

    # charge balance
    F[6] = (ec*Lcat + es*Lsep) * 2*(CS8m + CS6m + CS4m + CS2m + CSm2 + CS3/2) +
            epcLi2S/VLi2S*Lcat*2 + epsLi2S/VLi2S*Lsep*2 - S_bal*2*dod/100

    # electroneutrality
    F[7] = CAref + 2*(CS8m + CS6m + CS4m + CS2m + CSm2) + CS3 - CLi

    # S mass balance
    F[8] = (ec*Lcat + es*Lsep) * (CS8l*8 + CS8m*8 + CS6m*6 + CS4m*4 + CS2m*2 + CSm2*1 + CS3*3) + epcS8/VS8*Lcat*8 + epsS8/VS8*Lsep*8 + epcLi2S/VLi2S*Lcat*1 + epsLi2S/VLi2S*Lsep*1 -S_bal

    # Precipitation
    F[9] = CS8l - KspS8

    # Radical chemistry
    F[10] = CS3^2/CS6m - KS3
end


function thermoModel2_S3!(F,x)
    """
    Equations for region 2: no solids
    With S3 radical chemistry
    """

    global dod
    # variables to solve
    CLi = exp(x[1])
    CS8l = exp(x[2])
    CS8m = exp(x[3])
    CS6m = exp(x[4])
    CS4m = exp(x[5])
    CS2m = exp(x[6])
    CSm2 = exp(x[7])
    phi1 = exp(x[8])
    CS3 = exp(x[9])

    # porosites
    epcS8 = 0
    epcLi2S = 0
    epsS8 = 0
    epsLi2S = 0
    ec = 1-epcS8-epcLi2S-efc
    es = 1-epsS8-epsLi2S-efs

    # potential equations, setting overpotentials=0
    F[1] = phi1 - U2 + f *(sS8m_2*log(CS8m/1000) + sS8l_2*log(CS8l/1000) )
    F[2] = phi1 - U3 + f *(sS6m_3*log(CS6m/1000) + sS8m_3*log(CS8m/1000) )
    F[3] = phi1 - U4 + f *(sS4m_4*log(CS4m/1000) + sS6m_4*log(CS6m/1000) )
    F[4] = phi1 - U5 + f *(sS2m_5*log(CS2m/1000) + sS4m_5*log(CS4m/1000) )
    F[5] = phi1 - U6 + f *(sSm2_6*log(CSm2/1000) + sS2m_6*log(CS2m/1000) )

    # charge balance
    F[6] = (ec*Lcat + es*Lsep) * 2*(CS8m + CS6m + CS4m + CS2m + CSm2 + CS3/2) +
            epcLi2S/VLi2S*Lcat*2 + epsLi2S/VLi2S*Lsep*2 - S_bal*2*dod/100

    # electroneutrality
    F[7] = CAref + 2*(CS8m + CS6m + CS4m + CS2m + CSm2) + CS3 - CLi

    # S mass balance
    F[8] = (ec*Lcat + es*Lsep) * (CS8l*8 + CS8m*8 + CS6m*6 + CS4m*4 + CS2m*2 + CSm2*1 + CS3*3) +  epcS8/VS8*Lcat*8 + epsS8/VS8*Lsep*8 + epcLi2S/VLi2S*Lcat*1 + epsLi2S/VLi2S*Lsep*1 - S_bal

    # Radical chemistry
    F[9] = CS3^2/CS6m - KS3

end


function thermoModel3_S3!(F,x)
    """
    Equations for region 3: only Li2S solid
    With S3 radical chemistry
    """

    global dod
    # variables to solve
    CLi = exp(x[1])
    CS8l = exp(x[2])
    CS8m = exp(x[3])
    CS6m = exp(x[4])
    CS4m = exp(x[5])
    CS2m = exp(x[6])
    CSm2 = exp(x[7])
    phi1 = exp(x[8])
    epcLi2S = exp(x[9])
    CS3 = exp(x[10])

    # porosites
    epcS8 = 0
#     epcLi2S = 0
    epsS8 = 0
    epsLi2S = 0
    ec = 1-epcS8-epcLi2S-efc
    es = 1-epsS8-epsLi2S-efs

    # potential equations, setting overpotentials=0
    F[1] = phi1 - U2 + f *(sS8m_2*log(CS8m/1000) + sS8l_2*log(CS8l/1000) )
    F[2] = phi1 - U3 + f *(sS6m_3*log(CS6m/1000) + sS8m_3*log(CS8m/1000) )
    F[3] = phi1 - U4 + f *(sS4m_4*log(CS4m/1000) + sS6m_4*log(CS6m/1000) )
    F[4] = phi1 - U5 + f *(sS2m_5*log(CS2m/1000) + sS4m_5*log(CS4m/1000) )
    F[5] = phi1 - U6 + f *(sSm2_6*log(CSm2/1000) + sS2m_6*log(CS2m/1000) )

    # charge balance
    F[6] = (ec*Lcat + es*Lsep) * 2*(CS8m + CS6m + CS4m + CS2m + CSm2 + CS3/2) +
            epcLi2S/VLi2S*Lcat*2 + epsLi2S/VLi2S*Lsep*2 - S_bal*2*dod/100

    # electroneutrality
    F[7] = CAref + 2*(CS8m + CS6m + CS4m + CS2m + CSm2) + CS3 - CLi

    # S mass balance
    F[8] = (ec*Lcat + es*Lsep) * (CS8l*8 + CS8m*8 + CS6m*6 + CS4m*4 + CS2m*2 + CSm2*1 + CS3*3) + epcS8/VS8*Lcat*8 + epsS8/VS8*Lsep*8 + epcLi2S/VLi2S*Lcat*1 + epsLi2S/VLi2S*Lsep*1 - S_bal

    # Precipitation
    F[9] = CLi^2 * CSm2 - KspLi2S

    # Radical chemistry
    F[10] = CS3^2/CS6m - KS3

end



function solveThermoAll_S3!(sol_mat, x0_99, x0_1; tol=1e-3, step=1)

    """
    Algorithm:
    Solve Region 3 in descending order, starting from 99% DOD till no more solution found.
    Solve Region 1 in ascending order, starting from 1% DOD till no more solution found.                   Solve Region 2 in ascending order, starting from where Region 1 ended till where Region 3 begun.

    Stores all solutions in sol_mat, a len(dod_arr) x Nvar=11 matrix
    x0_99: initial guess for 99% dod,
    x0_1: initial guess for 1% dod
                                                
    also return dod_r2LB and dod_r2UB which are the lower and upper bounds of region 2, inclusive 
    """

    #### RUNNING REGION 3 #####
    #print("\nREGION 3\n")
    global x0 = copy(x0_99)
    append!(x0,CS3_0)
    dod_arr = 99:-step:1
    thermoModel = thermoModel3_S3!
    runThermoRegion!(dod_arr, x0, thermoModel,sol_mat, tol,step)
    dod_r2UB = dod   # dod where region 3 begins and region 2 ends


    # #### RUNNING REGION 1 #####
    #print("\nREGION 1\n")
    global x0 = copy(x0_1)
    append!(x0,CS3_0)
    dod_arr = 1:step:99
    thermoModel = thermoModel1_S3!
    runThermoRegion!(dod_arr, x0, thermoModel, sol_mat, tol, step)
    dod_r2LB = dod


    # ### RUNNING REGION 2 #####
    #print("\nREGION 2\n")
    dod_arr = dod_r2LB:step:dod_r2UB
    if dod_r2LB==1
        global x0 = x0_1[1:8]
    else         # look at element before  the LB for initial guess
        idx  = findIdx(dod,step)                                        
        global x0 = log.(sol_mat[idx-1,1:8])                                        
    end
    append!(x0,CS3_0)
    thermoModel = thermoModel2_S3!
    runThermoRegion!(dod_arr, x0, thermoModel, sol_mat, tol, step)


    return dod_r2LB, dod_r2UB
end
