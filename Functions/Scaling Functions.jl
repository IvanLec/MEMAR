function scaling(score_met, score_not_met, rat, NN , L)
    """Function to compute the scaling of a given MEMAR model with respect to different percentage of methylated sites.
    The function takes as input:
        1) score_met: the predicted probability of true methylated sites 
        2) score_not_met: the predicted probability of true non-mehtylated sites 
        3) rat: the different ratios between the methylated and non-methylated class considered 
        4) NN: number of replicates to compute averages and standard deviations 
        5) L: the number of fixed mehtylated sites considered
    """
    #Initialization of the vector containng all AUPRCs for each of the NN replicate considered
    Matrix_A=[]
    #Looping over the amounts of replicates to compute average quantities
    for i in 1:NN
        AUPRC=[] #Intialization of the vector containinf the AUPRC for each given ratio
        #Looping over the different ratio bewteen methylated and non-methylated classes
        for r in rat 
            #Computing the amount of non-methylated sites with respect to the L methylated ones
            L1=ratio(r)*L
            #At each different ratio of methylayed sites, I cosider L fixed, but randomly chosen mehtylated sites, and L1, non-methylated ones, such test_HELA_STM2457
            #L/(L1+L)=r
            push!(AUPRC, plot_PR_score(shuffle!(vcat(score_met...))[1:L], shuffle!(vcat(score_not_met...))[1:L1],  "g"," ";verbose=true)[3])
        end
        #Adding to 
        push!(Matrix_A, AUPRC)
    end
    #Reshaping the vector of vector to have a matrix, as rows each ratio value, and as columns each of the replicates
    M=(reshape(hcat(Matrix_A...), length(rat), NN))
    #Computing the average for each ratio (over rows)
    av=[mean(M[i, :]) for i in 1:length(rat)]
    #Computing the standard deviation for each ratio (over rows)
    st=[std(M[i, :]) for i in 1:length(rat)]
    return av, st 
end

function ratio(k)
    """Function to compute the amount of non-methylated sites with respect to methylated ones
        The function takes as input: 
            1) k: the ratio of non-methylated sitez
    """
    #Example: if the amount of methylated sites is 10% then k=0.1, the function return 9, corresponding 
    #to the amount of times the non-methylated sites must be compared to the methylated sites to have as methylated ratio 0.1
    return floor(Int64, (1-k)/k)
end

function AUC_t(time, Kon, Koff, list, N)
    """ Funciton to compute the Aread ander the precision recall curve, of the two-state model considering different time instances.
    The function takes as input:
        1) time: a vector of time values for the two-state model
        2) Kon:the set of m6A-writers binding constants;
        3) Koff: the set of m6A-ereasers binding constants;
        4) list: the list of motifs in a given dataset to test the model 
        5) N: a vector containing the amount of each motif in a given dataset ordered by list_mots
    """
    #Initialization of the AUPRC and class imbalance vector
    AUC, ratio=[],[]
    #looping over the time 
    for t in time 
        #Computing for a given time-value the probability of the two-state model
        p=ρ(Kon,Koff,t)
        #Computing the area under the precision recall curve for the two-state model at the given time
        push!(AUC, pr(p,list, N)[3])
        #Computing the amount of methylated motifs predicted by the two-state model
        push!(ratio, sum(p.*N)/sum(N))
    end
    return AUC, ratio 
end 

function ρ(kon, koff, t)
    """ Function that computes the two-state model probability
    The fucntion takes as input: 
    1) kon: the set of m6A-writers binding constants;
    2) koff: the ser of m6A-ereasers binding constants;
    3) t: the time 
    The function returns: 
        a vector of probability for each motif at a given time t 
    """
    return (1 .-exp.(-t*(kon.+koff))).*(kon./(kon.+koff))
end

function pr(p, list_mots, N)
    """Function to compute the precition recall curve and the Area under it. 
    The function takes as input: 
        1) p: a vector of probability for each motif at a given time t 
        2) list_mots: the list of motifs in a given dataset
        3) N: a vector containing the amount of each motif in a given dataset ordered by list_mots
    The function returns:
    1) precision: the values of the precision at each given treshold value
    2) recall: the of the recall at each given treshold value
    3) AUPRC: the area under the precision recall curve
"""
    #Computing the True Positive Rate of the two-state model
    tpr=TPR(p, N)
    #Computing the recall of the two-state model
    recall=[sum(tpr[1:i]) for i in 1:length(list_mots)] 
    #Computing the precision of the two-state model
    precision=[sum(p[1:i].*N[1:i])/sum(N[1:i]) for i in 1:length(list_mots)]
    #Adding as default the first point of the precition recall curve, corresponding to (0, 1)
    precision, recall=vcat([1, precision]...), vcat([0, recall]...)
    #Computing the area under the precision recall curve
    AUC=trapz(recall,precision) 
    return precision, recall, AUC
end

function TPR(rho, N)
    """The function compute the True positive rate for the given set of N methylated motifs.
    The function takes as input: 
        1) rho: a vector of probability for each motif at a given time t 
        2) N: a vector containing the amount of each motif in a given dataset 
    The function returns:
            1) the true positive rate of the given two-state model predictions
    """
    return N.*rho/sum(rho.*N)
end

function probability_8DRACH(MEMARC,Dict_DRACH,list_mots_D )
    """Function  to compute all MEMAR-C scores, for 8-mers
    The function takes as input: 
        1) MEMARC: a datatframe containing the liste of motifs and MEMAR-C predicted probabilities in decreasing order fro all 8-mers 
        2) Dict_DRACH: a dictionary having as keys the list of DRACH motifs.
        3) list_mots_D: the DRACH motifs ordered by MEMAR-C marginalized probabilities
    The function returns:
        1) p0: the probability of all 8-mers DRACH motifs following their score 
        2) list_mots: the list of all 8-mers DRACH motifs following their score  
        3) pos: the corresponding position of list_mots over MEMARC DataFrame 
        4) ind: the indexes of the 8-mers DRACH motifs, following list_mots_D order
    """
    list_mots,pos, ind=[],[],[]
    for i in 1:size(MEMARC)[1]
        if haskey(Dict_DRACH, MEMARC[i, 1][3:7])==true
            push!(list_mots, MEMARC[i, 1])
            push!(pos, i)
        end
    end
    p0=MEMARC[:, 2][pos];
    for i in 1:18
        temp=[]
        [(list_mots[j][3:7]== list_mots_D[i] ? push!(temp, j) : nothing ) for j in 1:length(list_mots)] 
        push!(ind, temp)
    end
    return p0, list_mots, pos, ind 
end


function diff_8mers(rho_t ,ind, N)
    """Function to compute the proprortion of methylated motifs of the two-state model probabilities 
    The function takes as input: 
        1) rho_t: the two state model probability of methylation for each 8-mers DRACH motif at a given time t 
        2) ind: the index of all 8-mers containg a common DRCAH motif 
        3) N: the amount of each DRACH motifs (ordered by MEMAR-C predicted probabilities)
        The function returns:
            The proportion of the DRACH methylated motifs, following MEMAR-C predicted prbabilities
        """
    return [sum((rho_t[ind[i]].*N[ind[i]])) for i in 1:length(ind)]/sum([sum((rho_t[ind[i]].*N[ind[i]])) for i in 1:length(ind)])
end

function eq_prob_MarG(test, train, list_mots_D)
    """
    Function to train MEMAR-C on altered control/methylation datasets.
    The function takes as input:
        1) test: all sequences in control/altered methylation  conditions 
        2) train: the set of sequences used for training for MEMAR-C in standard conditions
        3) list_mots_D: the list of DRACH motifs ordered by MEMAR-C probability in standard conditions
    The function return: 
        A vector of three elements, containing predicted probabilities of all DRACH motifs of MEMAR-C models in a subpartitions of three independent datasets
    """
    #The train is done amongst only the transcripts which were used for training
    gene_train=vcat(intersect(keys(test), keys(train))...)
    #Length of the training sequences in the given control/altered methylation dataset
    N = length(gene_train)
    #Sub-partition into three independent dataset
    gene_partition=[gene_train[1:round(Int, N/3)],gene_train[round(Int, N/3)+1:round(Int, 2*N/3)], gene_train[round(Int, 2*N/3)+1:end]]
    #Returning trained models for each of the three sub-partition
    return [eq_prob_seed_Marg(test, gene_partition[i] , list_mots_D) for i in 1:3 ]
end


function eq_prob_seed_Marg(train, gene_partition, list_mots_D)  
    """ Function to train MEMAR-C model in a given dataset, returning the probability of 8-mers DRACH motifs.
    The function takes as input:
        1) train: a dataset
        2) gene_partition: the set of training genes of the given dataset
        3) list_mots_D: the list of DRACH motifs ordered by MEMAR-C probability in standard conditions
    The function returns:
        1) score_DRACH: the predicted MEMAR-C probability of DRACH motifs
    """     
    #Distinguishing all motifs into methylated and non-methylated, saving the corresponding stoichiometry
    train_met, ratio_train , gene_train_met, train_not, gene_train_not=crete_set_all(train, gene_partition)
    #Training the methylated dataset
    M5=strucutre_model(train_met , ratio_train ; Lseq=8, tolerance=0.00001,learning_rate=2, NP=1000, λ=0.001)
    #Training the non-methylated dataset selecting a fixed amount of non-methylated sites
    M6=strucutre_model(shuffle!(train_not[1:1200000]), ones(1200000); Lseq=8, tolerance=0.00001,learning_rate=2, NP=1000, λ=0.001)  
    Ava=[] #Initializing the average percentage of methylated sites
    [push!(Ava, train[a].perc_met) for a in gene_partition] #Saving all methylated percentage for each training gene
    Ava=mean(Ava) #Computing the consatnt prior as the average amongst all training geens
    #Scoring all the 8-mers
    #For eahc DRACH motifs, the marginalized probability to be methylated, following MEMAR-C, is computed
    score_DRACH=[marginalized_prob_MET(M5, M6, Ava/100, b[1:2]*b[4:5]) for b in list_mots_D]
    return score_DRACH
end

function marginalized_prob_MET(M_met, M_not_met, prior, seq)
    """Function to marginalize the probability from 8-mers DRACH motifs, to DRACH motifs.
    The funciton takes as input:
        1) M_met: MEMAR trained model amongst methylated sites
        2) M_not_met: MEAMR trained over non-methylated sites
        3) prior: the constant prior corresponding to the average percentage of methylated sites in the training data 
        4) seq: the 8-mer sequence contaning in thecentral position a DRACH motif
    The function returns:
        1) M: the probability of methylated DRACH motifs after marginalization
    """
    #Marginalizing the probability of a DRACH motifs being in the methylated dataset
    A1=sum([exp(energy_tri(a[1:2]*seq*a[3:4], M_met.h, M_met.J, M_met.K))/exp(M_met.logZ) for a in mots_list4])
    #Marginalizing the probability of a DRACH motifs being in the non-methylated dataset
    A2=sum([exp(energy_tri(a[1:2]*seq*a[3:4], M_not_met.h, M_not_met.J, M_not_met.K))/exp(M_not_met.logZ) for a in mots_list4])
    #score methylated
    m=A1*prior
    #score non methylated
    u=A2*(1-prior)
    #Normalzation constant
    s=m+u
    #Methylated probability
    M=m/s
    return M
end

function perc_met_dataset(Met, Not_met)
    """Functio to compute the percetage of methylated sites.
    The function takes as input:
        1) Met: a vector of vectors containing the probability to be methylated of true methylated sites
        2) Not_met: a vector containing the probability to be methlated of truly not-methylated adenosines
        The function returns:
            the percentage of methylated sites
    """
    #Collecting the elements of the vector of vectors into one single one
    met=vcat(Met...)
    not_met=vcat(Not_met...)
    #Computing the percentage of methylated sites amongst all adenosine sites
    return length(met)/(length(met)+length(not_met))
end

function Infer_uniform(p_infernece, p_target, c; F=[ones(3) for _ in 1:18], lower=vcat([0.01]), upper=vcat([10.0]), k_factor=1 )
    """Function to compute the uniform-rscaling parameter of the two-state model to minimize the L1-error with a target probability.
    The function takes as input:
        1) p_inference: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in standard conditions
        2) p_target: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in a given target dataset
        3) c: the time of the two-state model 
        An addtional parameter are passed:
        1) F: represents a rescaling of the writer-binding constant for each DRACH motifs, if the initial conditions are not represented by MEMAR-C in standard conditions
        The function returns: 
        1) result: the uniform-rescaled parameter minimizing the L1-error
        2) lower: a vector containing the one-element lower bound of the uniform rescaling, initialized to [0.1]
        3) upper: a vector containing the one-element upper bound of the uniform rescaling, initialized to [10.0]
        4) k_factor: a value which modify the ereasers binding consatnt, it is initialized to 1

        The function returns:
        1) result: the uniform-rescaling paramter minimizing the L1-error in the domain between lower and upper
    """
    #Initialization fo the vector of rescaling paramters of the two-state model
    result=[]
    #Looping over the three-subpartitions of the target-dataset
    for j in 1:3 
        #Defying the loss-function of the target probability and the two-state model, as the sum for each of the 18 DRACH motifs
        loss(α)=sum([abs(ρ(p_infernece[i]*F[i][j]*α[1], (1-p_infernece[i])*k_factor, c)-p_target[j][i]) for i in 1:18])
        #Computing the minimum of the loss-function 
        push!(result ,  Optim.minimizer(optimize(loss, lower, upper,[1.2], Fminbox(LBFGS())))[1])
    end
    return result
end

function parameter_minimization(p_infernece, p_target, c, i, j; F=[ones(3) for _ in 1:18] , xmin=0, pas=0.000001, xmax=2 , k_factor=1)
    """Function to compute the non-uniform rescaling parameter of the two state model to mathc a tager probability of methylated sites from a given dataset.
    The function takes as input:
        1) p_inference: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in standard conditions
        2) p_target: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in a given target dataset
        3) c: the time of the two-state model 
        4) i: i-th DRACH motifs amongst the 18 possible
        5) j: index representing the independent subset amongst three of the target dataset
    Some additional arguments are passed:
        1) F: represents a rescaling of writer-binding constant for a given DRACH motifs, if the initial conditions are not represented by MEMAR-C in standard conditions
        2) xmin: the initial possible value of the uniform-rescaling constant, it is initialized to '0'
        3) pas: the step size from 'xmin' to change the uniform-rescaling, it is initialized to '0.000001'
        4) xmax: the final possible value of the uniform-rescaling constant, it is initialized to '2'
        5) k_factor: a value which modify the ereasers binding consatnt, it is initialized to 1
    The function returns:
        the non-uniform rescaling paramater, from xmin to xmax, that minimize the L1-error between the two-state model uniformly-rescaled and the target MEMAR-C probability
    """
    #If the two-state model in initial conditions is different from standard conditions, then the writer binding constant are rescaled according to 'F'
    a=p_infernece[i]*F[i][j]
    #If the two-state model in initial conditions is different from standard conditions, then the ereaser binding constant are rescaled according to 'F'
    b=(1-p_infernece[i])*k_factor
    #Saving the given target probability of the j-th subpartition, and of the i-th DRACH motif
    d=p_target[j][i]
    #Obtaining the value which minimized the L1 error between the target probability and the two-state model one of the i-th DRACH motifs in the given j-th subset
    return vcat(xmin:pas:xmax...)[argmin([LL(x, a, b, c, d ) for x in vcat(xmin:pas:xmax...)])]
end


function LL(x, a, b, c, d)
    """Function to compute the L1-error between two probabilities:
    The function takes as input: 
        1) x: the writer-binding constant 
        2) a: the non-uniform rescaling parameter for the given motif probability
        3) b the ereaser-binding constant 
        4) c: the time of the two-state model 
        5) d: the target probability
    The function returns: 
        The L1-error  between the two-state model probability and the target one 
        """
    return abs(a*x/(a*x+b)*(1-exp(-c*(a*x+b)))-d)
end

function parameter_p(p_infernece, c, i, j, x_min; F=[ones(3) for _ in 1:18]  ,verb=false, k_factor=1)
    """Function to compute the probability of the-two state model after uniform/non-uniform rescaling of the givent i-th DRACH motifs, for the givent j-th subpartition of the dataset.
    The function takes as input: 
        1) p_inference: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in standard conditions
        2) c: the time of the two-state model 
        3) i: i-th DRACH motifs amongst the 18 possible
        4) j: index representing the independent subset amongst three of the target dataset
        5) x_min: the rescaling parameter minimizing the L1-erro with respec to a target probability
        Some additional optional arguments are passed: 
        1) F: a rescaling of the writer-binding constant for a given DRACH motifs, if the initial conditions are not represented by MEMAR-C in standard conditions
        2) verb: a boolean variable to specify if x_min represents a uniform-rescaling or non-uniform one
        3) k_factor: a value which modify the ereasers binding consatnt, it is initialized to 1
        The function returns:
            The two state-probability of the i-th DRACH motif form the rescaled two-state model 
    """
    #Modyfing the writer-binding constants
    a=p_infernece[i]*F[i][j]
    #Computin the ereaser binding constant 
    b=(1-p_infernece[i])*k_factor
    #Returning the two-state probability either with the rescaled writer binding paramters either unifromly or non-unfiromly
    (verb==false ? (return  fun(x_min[i][j],a, b, c)) : (return fun(x_min[j],a, b, c)))
end



function fun(x, a, b, c)
    """ Function that computes the two-state model probability rescaled unifromly or non-uniformly.
     The function takes as input:
        1) x: the rescaling of the writer binding constants 
        2) a: the set of m6A-writers binding constants;
        3) b: the ser of m6A-ereasers binding constants;
        4) c: the time 
    The function returns: 
        a vector of probabilities for each motif at a given time t (having rescaled writer-binding constants)
        """
    return a*x/(a*x+b)*(1-exp(-c*(a*x+b)))
end



function explanation_degree_param(N, C, tresh_c=[0.005,0.015,0.03])
    """Function to compute the compataibility between MEMAR-C over a dataset and the two-state model predictions 
    The function takes as input: 
        1) N: two-state model predicted probability
        2) C: MEMAR-C probability on a given dataset
        An optional argument is passed: 
        1) tresh_c: the values to distinghuis models into storngly compatible, compatible,
                    incompatible, and strongly incompatible, intialized to : [0.005,0.015,0.03]
        The function returns: 
        1) alp: a vector containing the intensity of the color which will be used in the plot to visualize the compatibility between the two models
    """
    #computing lower bound of two-state model probabilities
    a=mean.(N).-std.(N)
    #computing upper bound of two-state model probabilities
    b=mean.(N).+std.(N)
    #computing lower bound of MEMAR-C probabilities
    c=mean.(C).-std.(C)
    #computing upper bound of MEMAR-C probabilities
    d=mean.(C).+std.(C)
    #Initializing vector of colors, all to 0, corresponding to strong compatibility
    alp=zeros(18)
    #if the distance between the intervals is higher or equal the the first treshold I set all sites to compatible
    alp[findall(x->x==true, [interval_distance((a[i], b[i]), (c[i], d[i])) >=tresh_c[1] for i in 1:18])].=0.1
    #if the distance between the intervals is higher or equal the the second treshold I set all sites to incompatible
    alp[findall(x->x==true, [interval_distance((a[i], b[i]), (c[i], d[i])) >=tresh_c[2] for i in 1:18])].=0.5
    #if the distance between the intervals is higher or equal the the third treshold I set all sites to incompatible
    alp[findall(x->x==true, [interval_distance((a[i], b[i]), (c[i], d[i])) >=tresh_c[3] for i in 1:18])].=1
    return alp
end

function interval_distance((a, b), (c, d)) 
    """Function to compute the distance between predicted probabilities with their standard deviations
    The function takes as input: 
        1) (a, b):  the lower and upper bound of the first model 
        2) (c, d): the lower und upper bound of the second model 
    The function returns 
        The distance between the two intervals.
    """
    #Returning the distance between the two intervals if they do not intersect, otherwise returning 0 as distance
    return  max(0.0, a - d, c - b)
end


function Infer_uniform_FTO(p_infernece, p_target, c; F=[ones(3) for _ in 1:18] ,  lower=vcat([1.0]) , upper=vcat([100.0]), k_factor=1 )
    """Function to compute the uniform-rscaling parameter of the two-state model to minimize the L1-error with a target probability.
    The function takes as input:
        1) p_inference: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in standard conditions
        2) p_target: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in a given target dataset
        3) c: the time of the two-state model 
        An addtional parameter are passed:
        1) F: represents a rescaling of the writer-binding constant for each DRACH motifs, if the initial conditions are not represented by MEMAR-C in standard conditions
        The function returns: 
        1) result: the uniform-rescaled for ereaser-binding constant parameter minimizing the L1-error
        2) lower: a vector containing the one-element lower bound of the uniform rescaling for the ereaser, initialized to [1.0]
        3) upper: a vector containing the one-element upper bound of the uniform rescaling for the ereaser, initialized to [100.0]
        4) k_factor: a value which modify the ereasers binding consatnt, it is initialized to 1
        The function returns:
        1) result: the uniform-rescaling paramter minimizing the L1-error in the domain between lower and upper
    """                        
    #Initialization fo the vector of rescaling paramters of the two-state model
    result=[]
    #Looping over the three-subpartitions of the target-dataset
    for j in 1:3
        #Definition of the loss-function of the target probability and the two-state model, as the sum for each of the 18 DRACH motifs
        loss(α)=sum([abs(ρ(p_infernece[i]*F[i][j], (1-p_infernece[i])*k_factor*α[1], c)-p_target[j][i]) for i in 1:18])
        #Computing the minimum of the loss-function 
        push!(result ,  Optim.minimizer(optimize(loss, lower, upper,[50.0], Fminbox(LBFGS())))[1])
    end
    return result
end


function parameter_minimization_FTO(p_infernece, p_target, c, i, j; F=[ones(3) for _ in 1:18] , xmin=50, pas=0.001, xmax=100, k_factor=1)
    """Function to compute the non-uniform rescaling parameter of the two state model to mathc a tager probability of methylated sites from a given dataset.
    The function takes as input:
        1) p_inference: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in standard conditions
        2) p_target: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in a given target dataset
        3) c: the time of the two-state model 
        4) i: i-th DRACH motifs amongst the 18 possible
        5) j: index representing the independent subset amongst three of the target dataset
    Some additional arguments are passed:
        1) F: represents a rescaling of writer-binding constant for a given DRACH motifs, if the initial conditions are not represented by MEMAR-C in standard conditions
        2) xmin: the initial possible value of the uniform-rescaling constant, it is initialized to '0'
        3) pas: the step size from 'xmin' to change the uniform-rescaling, it is initialized to '0.000001'
        4) xmax: the final possible value of the uniform-rescaling constant, it is initialized to '2'
        5) k_factor: a value which modify the ereasers binding consatnt, it is initialized to 1
    The function returns:
        the non-uniform rescaling paramater, from xmin to xmax, that minimize the L1-error between the two-state model uniformly-rescaled and the target MEMAR-C probability
    """
    #If the two-state model in initial conditions is different from standard conditions, then the writer binding constant are rescaled according to 'F'
    a=p_infernece[i]*F[i][j]
    #If the two-state model in initial conditions is different from standard conditions, then the ereaser binding constant are rescaled by k_factor
    b=(1-p_infernece[i])*k_factor
    #Saving the given target probability of the j-th subpartition, and of the i-th DRACH motif
    d=p_target[j][i]
    #Obtaining the value which minimized the L1 error between the target probability and the two-state model one of the i-th DRACH motifs in the given j-th subset
    return vcat(xmin:pas:xmax...)[argmin([L_FTO(x, a, b, c, d ) for x in vcat(xmin:pas:xmax...)])]
end

function L_FTO(x, a, b, c, d)
    """Function to compute the L1-error between two probabilities:
    The function takes as input: 
        1) x: the ereaser-binding constant 
        2) a: the non-uniform rescaling parameter for the given motif probability
        3) b the ereaser-binding constant 
        4) c: the time of the two-state model 
        5) d: the target probability
    The function returns: 
        The L1-error  between the two-state model probability and the target one 
        """
        #Returning the L1-error between two-state model probability with the rescaled eraser binding ocnstant and the target probability
    return abs(a/(a+b*x)*(1-exp(-c*(a+x*b)))-d)

end

function parameter_p_FTO(p_infernece, c, i, j, x_min; F=[ones(3) for _ in 1:18]  ,verb=false, k_factor=1)
    """Function to compute the probability of the-two state model after uniform/non-uniform rescaling of the givent i-th DRACH motifs, for the givent j-th subpartition of the dataset.
    The function takes as input: 
        1) p_inference: MEMAR-C predicted probabilities amongst marginalized DRACH motifs in standard conditions
        2) c: the time of the two-state model 
        3) i: i-th DRACH motifs amongst the 18 possible
        4) j: index representing the independent subset amongst three of the target dataset
        5) x_min: the rescaling parameter minimizing the L1-erro with respect to a target probability
        Some additional optional arguments are passed: 
        1) F: a rescaling of the writer-binding constant for a given DRACH motifs, if the initial conditions are not represented by MEMAR-C in standard conditions
        2) verb: a boolean variable to specify if x_min represents a uniform-rescaling or non-uniform one
        3) k_factor: a value which modify the ereasers binding consatnt, it is initialized to 1
        The function returns:
            The two state-probability of the i-th DRACH motif form the rescaled two-state model 
    """
    #Modyfing the writer-binding constants
    a=p_infernece[i]*F[i][j]
    #Computin the ereaser binding constant 
    b=(1-p_infernece[i])*k_factor
    #Returning the two-state probability either with the rescaled ereasers binding paramters either unifromly or non-unfiromly
    (verb==false ? (return  fun_FTO(x_min[i][j],a, b, c)) : (return fun_FTO(x_min[j],a, b, c)))
end


function fun_FTO(x, a, b, c)
    """ Function that computes the two-state model probability rescaled unifromly or non-uniformly.
    The function takes as input:
       1) x: the rescaling of the ereaser binding constants 
       2) a: the set of m6A-writers binding constants;
       3) b: the ser of m6A-ereasers binding constants;
       4) c: the time 
   The function returns: 
       a vector of probabilities for each motif at a given time t (having rescaled writer-binding constants)
       """
       #Two-state model probability at th egiven time, having rescaled the ereaser bidning paramter
    return a/(a+b*x)*(1-exp(-c*(a+x*b)))
end



function s_a_f(N, p_, pp_)
    """ Function to compute the difference in the proportion of DRACH methylated motifs, from 2 two-state model predicted probabilities.
    The function takes as input: 
        1) N: The amount of DRACH methylated motif in the given dataset 
        2) p_: the probability of the second two-state model 
        3) pp_: the prpbability of the first two-state model 
    The function returns:
        The mean and the standard deviations of the difference between the proportion of the 2nd two-state model with respect to the 1st one
    """
    #Intialization
    A=[]
    #Considering a matrix form of the porbabilites of the first model (as rows the different subpartions)
    p_vec=hcat(p_...)
     #Considering a matrix form of the porbabilites of the second model (given the different subpartions)
    pp_vec=hcat(pp_...)
    #Saving the difference in the proportion of the methylated motifs
    [push!(A, diff_f(N, p_vec[i, :])- diff_f(N, pp_vec[i, :]) ) for i in 1:size(p_vec)[1]]
    #Matrix of the difference between proprotion of DRACH methlated motifs of the two motifs, having as columns the three different subpartions
    MA=hcat(A...)
    #Retunring the average difference of proportion, and the standard deviation, for each DRACH motif, considering the three different subpartitions.
    return mean.(eachrow(MA)), std.(eachrow(MA))
end

function diff_f(N, p_)
    """Function to compute the predicted proportion of methylated motifs given a two-state model probability. 
    The function takes as input: 
        1) N: the amount of DRACH motifs,
        2) p_: the probability of a two-state model 
        The function returns: 
        The porportion of methylated DRACH motifs following the two-state model predictions
    """
    return p_.*N./sum(p_.*N)
end