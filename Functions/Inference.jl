function crete_set_all(set, gene_label; verbose=true, window=4)
    """Function to create all methylated and non-methylated sites from experimental sites, taking 8-mers around an adenosine.
    The function takes as entry:
        1) set: A dictionary containing all the information of the genes
        2) gene_label: the set of genes considered to create the sites
    An optional argument is passes:
        3) verbose: a boolean variable. It takes teh default value 'true' to create the methylated site sets, whilest 'false' for non methylated ones.
    The function returns:
        1) Met_set: Methylated dataset
        2) ratio: stoichiometry for each methylated sequence in the dataset (if verbose is 'true')
        3) gene_met: list of the origin gene for each methylated sequence in the dataset
        4) Not_met_set: Non-ethylated dataset
        5) gene_not: list of the origin gene for each non-methylated sequence in the dataset
    """
    Met_set,ratio,Not_met_set,gene_met,gene_not=[],[],[],[],[]
    #Looking at the condition to save the stoichiometry 
    if verbose == true
        #Methylation dataset
        for a in gene_label #loop for all the genes
            #gene length
            L=length(set[a].sequence)
            #loop for all methylated sites
            for k in 1:length(set[a].pos_met)
                #Finding sites only if it is possible to define a central motif of length 2*window+1 at a given sites, extremes are left out
                if window<set[a].pos_met[k]<=L-window
                    #saving motifs with central adenosine (removing the A)
                    push!(Met_set,set[a].sequence[set[a].pos_met[k]-window:set[a].pos_met[k]-1]*set[a].sequence[set[a].pos_met[k]+1:set[a].pos_met[k]+window])
                    push!(ratio, set[a].ratio[k])
                    push!(gene_met,a)
                end
            end
        end
    else #the stoichiometry is not available, this value is not saved into the variable 'ratio'
        for a in gene_label #loop for all the genes
            #gene length
            L=length(set[a].sequence)
            #loop for all methylated sites
            for k in 1:length(set[a].pos_met)
                #Finding sites only if it is possible to define a central motif of length 2*window+1 at a given sites, extremes are left out
                if window<set[a].pos_met[k]<=L-window
                    #saving motifs with central adenosine (removing the A)
                    push!(Met_set,set[a].sequence[set[a].pos_met[k]-window:set[a].pos_met[k]-1]*set[a].sequence[set[a].pos_met[k]+1:set[a].pos_met[k]+window])
                    push!(gene_met,a)
                end
            end
        end
    end
    #non-Methylation dataset 
    for a in gene_label #loop for all the genes
        #gene length
        L=length(set[a].sequence)
        #loop for all non-methylated sites
        for k in 1:length(set[a].pos_not)
            #Finding sites only if it is possible to define a central motif of length 2*window+1 at a given sites, extremes are left out
            if window<set[a].pos_not[k]<=L-window
                #saving motifs with central adenosine (removing the A)
                push!(Not_met_set,set[a].sequence[set[a].pos_not[k]-window:set[a].pos_not[k]-1]*set[a].sequence[set[a].pos_not[k]+1:set[a].pos_not[k]+window])
                push!(gene_not,a)
            end
        end
    end
    if verbose == true  #the stoichiometry variable 'ratio' is returned
        return Met_set,ratio,gene_met, Not_met_set, gene_not
    else
        return Met_set, gene_met, Not_met_set, gene_not 
    end
end



function strucutre_model(train_met, weights; Lseq=8, tolerance=0.0001, max_iter=10000, learning_rate=2, gamma=0.8, NP=1000, λ=0)
    """Funcion to inference datasets made a set of sequences around central adenosine. 
    It takes as input: 
        1) train_met: a vector containing all sequences around central adenosines
        2) weights: a vector which account for the stoichiometry of methylated motifs (when available)
        Some additional quantities are defined:
        3) Lseq: the length of sequences, the deafault value is 8 
        4) tolerance: the condition to stop the gradient computation 
        5) max_iter: the maximum number of computation to extremise the function trhough the gradient computation
        6) learning_rate: a parameter to speed the convergence of the algorithm 
        7) gamma: a parameter to speed up algorithm computation through a 'memory' term (momentum)
        8) NP: number of iterations after which an information of the convergence point is given
        9) λ: L2-regularization strength parameter
    The function returns:
        1) A Structured of the type MaxEnt3, correspnding to the set of forces 'h', 'J', and 'K' (1,2,3-mers) and the log of the partition function 'lZ'.
    """
    #Computation of the forces form the given training dataset, and the logarithm of the partition function
    h,J,K, lZ=Inference(train_met, weights, Lseq, tolerance, max_iter, learning_rate, gamma, NP, λ);
    return MaxEnt3(h,J,K,lZ) 
end

function Inference(Sequences, weights, Lseq, tolerance, max_iter, learning_rate, gamma, NP, λ)
    """Function to perform the inference of the parameters of the model from a given training dataset. It takes as input:
        1) Sequences: the set of training sequences
        2) weights: the stoichiometry (if available) for each sequence 
        3) Lseq: the length of the sequences 
        4) tolerance: the condition to stop the gradient computation 
        5) max_iter: the maximum number of computation to extremise the function trhough the gradient computation
        6) learning_rate: a parameter to speed the convergence of the algorith 
        7) gamma: a parameter to speed up algorithm computation through a 'memory' term (momentum)
        8) NP: number of iterations after which an information of the convergence point is given
        9) λ: L2-regularization strength parameter
    The function returns:
        1)param_infer: the set of forces 'h', 'J', and 'K' (1,2,3-mers) and the log of the partition function 'lZ'.
    """
    #Computing the observables (1-mers, 2-mers, 3-mers)
    n_obs=count_Mots(Sequences, weights, Lseq)
    #Masking all nucleotides base on gauge 
    mask, list=masking(Lseq)
    n_obs=n_obs[mask.==true]
    #Initial starting point of the optimization
    x=rand(length(n_obs))
    #Extreming function
    param_infer=Optimization(x, list, max_iter, tolerance, n_obs, Lseq,learning_rate, gamma, NP, λ )
    return param_infer
end

function masking(L)
    """A function to mask all the dependent variables (so all the parameters that, due to the model constaints, are not needed to be directly inferrred.)
    It takes as argument:
        1) L: the size of the sequence to infer 
    The function returns: 
        1) mask: a vector containing the masked variables linked to 1-mers, 2-mers and 3-mers
        2) a vector containing the list of dependent parameters for each set
    """

    #Masking 1-mers
    mask_1mer=vcat(ones(1, L*length(dna_alphabet))...) #there are L*||alphabet|| parameters for nucleotides
    #all 'T' nucleotides are choose to be dependent variables
    mask_1mer[findall(x->x=="T", repeat(dna_alphabet, L))].=0 
    L1=dna_alphabet[findall(x->x!="T",dna_alphabet)] 

    #Masking 2-mers
    mask_2mer=vcat(ones(1, (L-1)*length(mots_list2))...) #there are (L-1)*||alphabet^2|| parameters for di-nucleotides
    #all motifs containing at least one 'T' nucleotides and the motifs 'AG' and 'GC' are choosen to be dependent variables
    mask_2mer[findall(x->(count('T', x)>0 || x=="AC" || x=="GC") ==true, repeat(mots_list2, L-1))].=0
    L2=mots_list2[findall(x->(count('T', x)==0 && x!="AC" && x!="GC"),mots_list2)]

    #Masking 3-mers
    mask_3mer=vcat(ones(1, (L-2)*length(mots_list3))...) #there are (L-2)*||alphabet^3|| parameters for tri-nucleotides
    #all motifs containing at least one 'T' nucleotides, 'CA', and the motifs "CCC" and "AAA" are choosen to be dependent variables
    mask_3mer[findall(x->(count('T', x)>0 || count("CA", x)>0 || x=="CCC" || x=="AAA") ==true, repeat(mots_list3, L-2))].=0 
    L3=mots_list3[findall(x->(count('T', x)==0 && count("CA", x)==0 && x!="CCC" && x!="AAA") ,mots_list3)]
    #Putting inside one vector all the masked motifs
    mask=vcat(mask_1mer,mask_2mer,mask_3mer...) 
    return mask, [L1, L2, L3]
end

function count_Mots(Sequences, weigths, Lseq) 
    """Function to compute the average weigthed frequency at each positon for a given dataset. It takes as input:
        1)Sequences: the set of training sequences
        2) weights: the stoichiometry (if available) for each sequence
        3) Lseq: the length of the sequences 
    The function returns:
        1) A vector containint the average weigthed frequency
    """
    M = permutedims(hcat([split(s, "")[1:end] for s in Sequences]...)) #Converting into a matrix of strings
    N=sum(weigths) #Weigthed number of sequences
    counting=[]
    #1-mer
    #looking for the nucleotide in the matrixes at all the rows, at the indicated position, summing the weights for all the found nucleotides
    [[push!(counting, (isempty(findall(x->x==a, M[:, pos])) ? 0 : sum(weigths[(findall(x->x==a, M[:, pos]))]))) for a in dna_alphabet] for pos in 1:Lseq] 
    #2-mers
    for pos in 1:Lseq-1
        #Cosnidering all the di-nucleotides (unifying two consecutive columns)
        M2= [row[pos] * row[pos+1] for row in eachrow(M)]
        #looking for the di-nucleotide in the matrixes at the indicated position, summing the weights for all the found di-nucleotides
        [push!(counting, (isempty(findall(x->x==ab, M2)) ? 0 : sum(weigths[(findall(x->x==ab, M2))]))) for ab in mots_list2]
    end
    #3-mers
    for pos in 1:Lseq-2
        #Cosnidering all the tri-nucleotides (unifying three consecutive columns)
        M3= [row[pos] * row[pos+1]* row[pos+2] for row in eachrow(M)]
        #looking for the tri-nucleotide in the matrixes at the indicated position, summing the weights for all the found tri-nucleotides
        [push!(counting, (isempty(findall(x->x==abc, M3))) ? 0 : sum(weigths[(findall(x->x==abc, M3))])) for abc in mots_list3]
    end
    return counting./N #Computing frequency with respect to weights
end

function Optimization(x, list, max_iter, tolerance, n_obs, Lseq, learning_rate, gamma, NP,λ)
    """Function to infer the model parameters, optimizing the computation through gradient computation with momenutm and automati differentation.
    The function takes as arguments:
        1) x: a vector containing all the independent parameters of the model (which need to be directly inferred)
        2) list: the list of all the indepent parameters
        3) max_iter: the maximum number of computation to extremise the function trhough the gradient computation
        4) tolerance: the condition to stop the gradient computation 
        5) n_obs: the average weitghed frequency of independent motifs
        6) Lseq: the length of the sequences 
        7) learning_rate: a parameter to speed the convergence of the algorith 
        8) gamma: a parameter to speed up algorithm computation through a 'memory' term (momentum)
        9) NP: number of iterations after which an information of the convergence point is given
        10) λ: L2-regularization strength parameter
    The function returns:
        1) the set of forces 'h', 'J', and 'K' (1,2,3-mers) and the log of the partition function 'lZ'.
    """
    #Definition of the function to minimize
    f(x)= compute_logZ(x, Lseq, list)-sum(n_obs.*x)+λ*sum(x.^2) #Definition of the function to extremize
    t=vcat(zeros(1, length(x))...) #initialization of empty vector
    for i in 1:max_iter 
        grad = ReverseDiff.gradient(f, x) #gradient computation
        t=(t * gamma)+ learning_rate * (grad) #gradient with momentuum
        x_new = x - t # Update all elements of x simultaneously
        if i%NP==0 #if 'NP' iterations have passed:
            println("At iteration $i the total error is $(sum(abs.((-grad) ))) with maximum error of $(maximum(abs.((grad) )))")
            flush(stdout)
        end 
        #Checking if the maximum gradient component is lower than the tolerence one
        if maximum(abs.((grad) )) < tolerance 
            println("Converged at iteration $i, where the maximum error is $(maximum(abs.((-grad) )))")
            x = x_new
            break
        end
        # Update x for next iteration
        x = x_new
        if i == max_iter #If the maximum iter has been reached withouth convergence:
            println("Reached maximum iterations $i without convergence.")
        end
    end
    #Saving the values of the forces into dictionaries for all independent variables
    h, J, K=vector_to_dict(x, Lseq)
    #Addinf for all dictionaries the list of dependent variables (due to the model constaints)
    h_,J_,K_=add_masked_keys(h, dna_alphabet, Lseq), add_masked_keys(J, mots_list2, Lseq-1), add_masked_keys(K, mots_list3, Lseq-2)
    return h_,J_,K_, compute_logZ(x, Lseq, list)
end

function compute_logZ(x, Lseq, list)
    """Function to compute the logarithm of the partition function. It takes as input:
    1) x: a vector containing all the independent parameters of the model (which need to be directly inferred)
    2) Lseq: the length of the sequence
    3) list: the list of all the indepent parameters
    """
    #Tranform vector of parameters into vocabularies
    return log(compute_Z(x, Lseq, list))
end

function delta(x,list, mot)
    """Delta Dirac function. It takes as input:
    1) x: a vector containing all the independent parameters of the model (which need to be directly inferred)
    2) list: the list of all the indepent parameters
    3) mot: the specific motif: sum_{i ∈ list} δ(x[i]-mot)
    The function returns: 
    1) The value of the parameter for the given motif sum_{i ∈ list} δ(x[i]-mot)=x[mot] iff mot ∈ list
    """
    #Delta Dirac for motives
    t=findfirst(x->x==mot, list)
    #If not found it return 0
    isnothing(t)==true ? (return 0) : (return x[t])
end

function compute_Z(x, Lseq, list) 
    """Compute the partition function of the model. It takes as input:
        1) x: a vector containing all the independent parameters of the model (which need to be directly inferred)
        2) Lseq: the length of the sequence
        3) list: the list of all the indepent parameters.
    The function returns 
        1) Z: the partition function of the model
    """
    #Initialization of the transfer matrix as identity one
    T=I
    #Throught the trasnfer matrix, a product of matrices is performed, from the last index (Lse1-1) to the first (1).
    indexes=reverse(vcat(1:(Lseq-1)...))
    for n in indexes #looping over the indexes
        T=transfer_matrix(x, n, Lseq, list)*T #product of transfer matrixes
    end 
    Z=sum(diag(T)) #The paritition function is the trace of the resulting trasnfer matrix
    return Z
end

function transfer_matrix(x, n, L, list) 
    """Function to compute the transfer matrix. It takes as input:
        1) x: a vector containing all the independent parameters of the model (which need to be directly inferred)
        2) n: the index of the transfer matrix 
        3) L: the length of the sequence
        4) list: the list of all the indepent parameters
    The function returns: 
        T: the transfer matrix at the given index
    """
    L1=list[1] #List of mono-nucleotide
    L2=list[2] #List of di-nucleotide
    L3=list[3] #List of tri-nucleotide
    L_alph=16
    #Initialization of the transfer matrix to 0.0 elements
    T=[]
    #Initialization to matrix of zeros 
    for i in 1:L_alph
        t=[]
        [push!(t, 0.0) for j in 1:L_alph] #Initialization empty vector
        push!(T, t)
    end
    T=reshape(vcat(T'...), L_alph, L_alph) #reshaping into an empty matrix
    for i in 1:L_alph
        for j in 1:L_alph
                #Differentiating between the last Transfer matrix and the rest
                if n!=L-1
                    if mots_list2[i][2]==mots_list2[j][1]
                        T[i, j]=exp(
                            (delta(x[(n-1)*length(L1)+1:length(L1)*n], L1, "$(mots_list2[i][1])")+delta(x[n*length(L1)+1:length(L1)*(n+1)], L1, "$(mots_list2[i][2])")+delta(x[(n+1)*length(L1)+1:length(L1)*(n+2)], L1, "$(mots_list2[j][2])"))/3+
                            (delta(x[(n-1)*length(L2)+1+length(L1)*L:length(L2)*n+length(L1)*L], L2, "$(mots_list2[i])")+delta(x[n*length(L2)+1+length(L1)*L:length(L2)*(n+1)+length(L1)*L], L2, "$(mots_list2[j])"))/2+
                            delta(x[(n-1)*length(L3)+1+length(L1)*L+length(L2)*(L-1):length(L3)*n+length(L1)*L+length(L2)*(L-1)],L3, "$(mots_list2[i]*mots_list2[j][2])"))
                    end
                else
                    T[i, j]=exp(
                        ((delta(x[(n-1)*length(L1)+1:length(L1)*n], L1,"$(mots_list2[i][1])"))+2*(delta(x[(n)*length(L1)+1:length(L1)*(n+1)], L1,"$(mots_list2[i][2])"))+2*(delta(x[1:length(L1)], L1,"$(mots_list2[j][1])"))+delta(x[length(L1)+1:length(L1)*2], L1,"$(mots_list2[j][2])"))/3+
                        (delta(x[(n-1)*length(L2)+1+length(L1)*L:length(L2)*n+length(L1)*L],L2, "$(mots_list2[i])")+delta(x[1+length(L1)*L:length(L2)+length(L1)*L],L2,"$(mots_list2[j])"))/2)
                end
        end
    end
    return T
end



function vector_to_dict(x, Lseq)
    """Function to add the vector of the forces into thres dictionaries for 1,2,3-mers. It takes as arguments:
        1) x: a vector containing all the independent parameters of the model (which need to be directly inferred)
        2) Lseq: the length of the sequence
    The function returns:
        1) h: the dictionary containing all independent 1-mers forces
        2) J: the dictionary containing all independent 2-mers forces
        3) K: the dictionary containing all independent 3-mers forces
    """
    h=[] #initialization to empty vector
    c=1 #counting indexes for independent variables
    #1-mers
    for i in 1:Lseq #looping over all possible positions of nucleotides
        h_=Dict() #empty dictionary
        for a in dna_alphabet #looping over all nucleotides
            if a!="T" #looking for all motifs different from dependent nucleotide parameters
                h_[a]=x[c]
                c+=1
            end
        end
        push!(h, h_)
    end
    #2-mers
    J=[] #initialization to empty vector
    for i in 1:Lseq-1 #looping over all possible positions of di-nucleotides
        J_=Dict() #empty dictionary
        for a in mots_list2 #looping over all di-nucleotides
            if count('T', a)==0 && a!="AC" && a!="GC" #looking for all motifs different from dependent di-nucleotide parameters
                J_[a]=x[c]
                c+=1
            end
        end
        push!(J, J_)
    end
    #3-mers
    K=[] #initialization to empty vector
    for i in 1:(Lseq-2)  #looping over all possible positions of tri-nucleotides
        K_=Dict() #empty dictionary
        for a in mots_list3 #looping over all di-nucleotides
            if count('T', a)==0 && count("CA", a)==0 && a!="CCC" && a!="AAA"  #looking for all motifs different from dependent tri-nucleotide parameters
                K_[a]=x[c]
                c+=1
            end
        end
        push!(K, K_)
    end
    return h,J,K
end

function add_masked_keys(param, list, L)
    """Function which add to the dictionaries containing all independent forces all the dependent ones.
    The function takes as arguments:
    1) param: the vocabulary of independent forces
    2) list: the of indepent parameters
    3) L: the position of the motifs 
    """
    #looking if the motifs are not present, and in this case adding a zero value (to the zero-Gauge)
    [[(haskey(param[pos], mot) ? nothing : param[pos][mot]=0) for mot in list] for pos in 1:L]
    return param
end