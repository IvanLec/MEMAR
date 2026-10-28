function mapping_position(set, gene_label)
    """The function maps the position of transcripts into an ideal sequence made of 33.3̄% 5'-UTR, 33.3̄% CDS and 33.3̄% 3'-UTR
    The function takes as entries:
        1)set: a dictionary containint test set
        2)genes: the set of genes for a subset0 of the set dictionary for which score want to be computed
    The function returns:
        -MAP: a dictionary containing mapped position into this ideal sequence for all coding mRNA
    """
    MAP=Dict()
    for a in gene_label #loop over all genes 
        if typeof(set[a].cds)==CDS #checking only coding mRNA
            #=When a gene is coding, I map the position of the coding regions into a sequence 
            made of 1/3 5'-UTR, 1/3-CDS, 1/3 3'-UTR =#
            utr5=set[a].cds.start-1 #length of 5'-UTR
            cds=set[a].cds.stop-set[a].cds.start+1 #length of CDS
            utr3=length(set[a].sequence)-set[a].cds.stop+1 #length of 3'-UTR
            #Creating the mapping with the same number of points between the new extrems and the previous ones
            map=vcat(vcat(range((1/3)/utr5, 1/3, utr5)), vcat(range(1/3+(1/3)/cds, 2/3, cds)),  vcat(range(2/3+(1/3)/utr3, 1, utr3)))
            #Adding to the dictionary the mapping vector
            MAP[a]=map
        end
    end
    return MAP
end 

function score_ps_ds_pm(M_met, M_not_met, set, genes, prior; verb=false, map=[], d_met=[], d_not=[], binning=0.001,N_window=8)
    """The function computes the score on a set of genes. It takes as input:
        1)M_met: The model learnt for methylated sequences
        2)M_not_met: The model learnt for non-methylated sequences
        3)set: a dictionary containint test set
        4)genes: the set of genes for a subset0 of the set dictionary for which score want to be computed
        5)prior: the prior for sequences
        Some additional infromation can be used
        6)verb: it is a boolean variable to include or not the distance m6A distribution, it is initialized as false
        7)d_met: the ditribution of m6A sites, it is a vector conaint the set of point [x,y], and it is initialized to zero
        8)d_not: the ditribution of non methylated A sites, it is a vector conaint the set of point [x,y], and it is initialized to zero
        9)binning: the size of the discretization for the distribution of m6A and A sites
        10)N_window: the size of short motives to score around adenosines. This is set by default equal to 8
    The function returns:
        1)score_met: the score assigned by the model to all m6A sites in the reference experiment
        2)score_not: the score assigned by the model to all non modified A sites in experiments
    """ 
    score_met,score_not_met=[],[] #Initialization of empty scores
    for a in genes #Loop for all the genes to score
        #We compute all the score for each site of the sequence
        S=score_genome_complete(M_met, M_not_met, prior[a].perc_met/100, set[a].sequence, N_window; verbose=verb, Mapping=map, dist_met=d_met, dist_not=d_not, bin=binning, gene=a)
        push!(score_met, S[set[a].pos_met]) #Selecting score of m6A sites
        push!(score_not_met, S[set[a].pos_not]) #Selecting score of A sites
    end
    return score_met,score_not_met
end


function score_genome_complete(M_met, M_not_met,prior, seq, N; verbose=false, Mapping=[], dist_met=[], dist_not=[], bin=0.001 , gene=[])
    """
    The function computes a score for all sites of a given sequence. It takes as imput:
        1)M_met: The model learnt for methylated sequences
        2)M_not_met: The model learnt for non-methylated sequences
        3)prior: the prior for the given sequence
        4)seq: the sequence to score
        5)N: the size window of short motives around adenosines
        6)vebose: a boolean variable either to consider distance distribution of m6A sites or A. It is initilaized to false
        7)Mapping: dictionary containing the mapping of the given gene
        8)dist_met: the ditribution of m6A sites, it is a vector conaint the set of point [x,y], and it is initialized to zero
        9)dist_not: the ditribution of non methylated A sites, it is a vector conaint the set of point [x,y], and it is initialized to zero
        10)bin: the size of the discretization for the distribution of m6A and A sites
        11)gene: the dictionary name related to the considered gene to test
    The function returns:
        1) x a vector of the same length of the sequence, with model prediction score of methylated sites
    """
    w=convert(Int64, (N)/2) #taking half window size
    x=[] #inizialization of score
    for i in 1:length(seq)
        if seq[i]=='A' #looking for all i sites
            #if the position i is not at extremities (so regions where a central motive of length N cannot be considered) the score is 0, otherwise the model score is computed
            (i<w+1 || i > length(seq)-w) ? push!(x, 0.0) : push!(x, normalized_prob_MET(M_met, M_not_met, prior, seq[i-w:i-1]*seq[i+1:i+w], verbose, Mapping, dist_met, dist_not, bin,i, gene))
        else
            #if it is not an A sites the score assigned is zero
            push!(x, 0.0)
        end
    end
    return x
end

function normalized_prob_MET(M_met, M_not_met, prior, seq, verbose, mapping, dist_met, dist_not, bin, ii, gene)
    """
    Function to score the sequence following MERAM models. it takes as input:
        1)M_met: The model learnt for methylated sequences
        2)M_not_met: The model learnt for non-methylated sequences
        3)prior: the prior for the given sequence
        4)seq: a small motive to assign a score
        5)verbose: a booelean to compute the score modulated by distance distribution
        6)Mapping: dictionary containing the mapping of the given gene at the position of the motive seq
        7)dist_met: the ditribution of m6A sites, it is a vector conaint the set of point [x,y], and it is initialized to zero
        8)dist_not: the ditribution of non methylated A sites, it is a vector conaint the set of point [x,y], and it is initialized to zero
        9)bin: the size of the discretization for the distribution of m6A and A sites
        10)ii: position of the mapping
        11)gene: the dictionary name related to the considered gene to test
    """
    A1=exp(energy_tri(seq, M_met.h, M_met.J, M_met.K))/exp(M_met.logZ)
    A2=exp(energy_tri(seq, M_not_met.h, M_not_met.J, M_not_met.K))/exp(M_not_met.logZ)
    if verbose==false 
        #Non modulating by distance distribution
        #score methylated
        m=A1*prior
        #score non methylated
        u=A2*(1-prior)
        #Normalzation constant
        s=m+u
        #Methylated probability
        M=m/s
    else
        #Modulating by distance distribution
        #score methylated
        m=A1*prior*(dist_met[2][argmin(abs.(mapping[gene][ii].-dist_met[1]))])*bin
        #score non methylated
        u=A2*(1-prior)*(dist_not[2][argmin(abs.(mapping[gene][ii].-dist_not[1]))])*bin
        #Normalzation constant
        s=m+u
        #Methylated probability
        M=m/s
    end
    return M
end

function avarege_distribution(set, genes, F, bin; label=1)
    """Function to compute the average localization profiles of experimental methylated sites or non-methylated adenosines. it takes as entry:
        1) set: a dictionary containint test set
        2) genes: keys of the dictionary of the considered testing transcripts 
        3) F: the window size to smooth the resulting distribution 
        4) bin: the binning size to compute the histogram, therefore the density profile, of sites
        5) label: a binary value which indicates with 1 (default values) to select methylated sites, or with 0, non-methylated sites 
    The function return 
        1) xmean: the spaced bins over which the density distrubution has been compute_F1_threshold
        2) ymean: the avareged density distribution over xmean
    """
    #Mappingthe position of sequences into an ideal sequence made of 33.3̄% 5'-UTR, 33.3̄% CDS and 33.3̄% 3'-UTR
    MAP=mapping_position(set, vcat(genes...))
    dist_met=[]
    #Saving wither methylated sites (label=1) or non-methylated sites (label=0)
    (label==1 ?  [push!(dist_met, MAP[a][set[a].pos_met]) for a in keys(MAP)] :  [push!(dist_met, MAP[a][set[a].pos_not]) for a in keys(MAP)])
    #Computing the density distribution through an histogram
    yy, xx=hist(vcat(dist_met...), histtype="step" ,bins=bin, density=true)
    #Taking the bins and the avrege density over a window of length 'F'
    xmean, ymean=xx[1:end-1], movmean(yy, F);
    return xmean, ymean
end

function avarege_distribution_DRACH(set, genes, F, bin; label=1)
    """Function to show the localization distribution of methylated DRACH motifs. It takes as entry:
        1) set: a dictionary containint test set
        2) genes: keys of the dictionary of the considered testing transcripts 
        3) F: the window size to smooth the resulting distribution 
        4) bin: the binning size to compute the histogram, therefore the density profile, of sites
        5) label: a binary value which indicates with 1 (default values) to select methylated sites, or with 0, non-methylated sites 
    The function return 
        1) xmean: the spaced bins over which the density distrubution has been compute_F1_threshold
        2) ymean: the avareged density distribution over xmean of DRACH motifs
    """
    #Mappingthe position of sequences into an ideal sequence made of 33.3̄% 5'-UTR, 33.3̄% CDS and 33.3̄% 3'-UTR
    MAP=mapping_position(set, vcat(genes...))
    dist_met=[]
    #Computing the density distribution through an histogram looking for only DRACH methylated motifs
    (label==1 ?  [push!(dist_met, MAP[a][intersect(mask_DRACH(set[a].sequence, Dict_DRACH), set[a].pos_met)]) for a in keys(MAP)] :  [push!(dist_met, MAP[a][intersect(mask_DRACH(set[a].sequence, Dict_DRACH), set[a].pos_not)]) for a in keys(MAP)])
    yy, xx=hist(vcat(dist_met...), histtype="step" ,bins=bin, density=true)
    #Taking the bins and the avrege density over a window of length 'F'
    xmean, ymean=xx[1:end-1], movmean(yy, F);
    return xmean, ymean
end


function AUPRC_hyperparams(train_met, ratio_train, train_not,ind_not,constant_prior, validation; pas=0.005)
    """Function to compute difference AUPRC values, through MEMAR, to show the variability in performance based on different tolerance value (to stop gradient computation) or L2-regularization strengh (to avoid overfitting)
    The function takes as input:
            1) train_met: dataset of methylated sites 
            2) ratio_train: stoichiometry of methylated sites 
            3) train_not: dataset of non-methylated sites 
            4) ind_not: indexes of 2 000 000 randomly chosen non-methylated sites used to train
            5) constant_prior: ratio between methylated and non-methylated sites
            6) validation: set of validating genes to compute performance
            7) MAP: Dictionary containint the mappin towords an ideal sequence made of 33.3̄% 5'-UTR, 33.3̄% CDS and 33.3̄% 3'-UTR
            8) pas: parameter to change treshold value to compute precision recall curve. Its default value is 0.005
        The fucntion returns:
            1) AU: A mtrix containing AUPRCs computed at diferent tolerence values (rows) and L2-regularization paramter values (columns)

    """
    AU=[]
    for tol in [0.001, 0.0001, 0.00001, 0.000001] #Set of considered for the gradient tolerance
        AA=[]
        for Λ in [0, 0.001, 0.01, 0.1] #Set of considered for the L2-regularization strength
            #Training at each couple of hyperparameter a model for methylated and non-methylated dataset
            M5=strucutre_model(train_met , ratio_train ; Lseq=8, tolerance=tol,learning_rate=2, NP=1000, λ=Λ)
            M6=strucutre_model(train_not[ind_not], vcat(ones(1,length(ind_not ))...); Lseq=8, tolerance=tol,learning_rate=2, NP=1000, λ=Λ)
            #Scoring validation sites with the simplest model, with constant prior 
            score_met_ps,score_not_met_ps=score_ps_ds_pm(M5, M6, validation, keys(validation), constant_prior)
            #Computing the Precision-Recall Curve
            x, y=PR(vcat(score_met_ps...), vcat(score_not_met_ps...); pass=pas)
            #Filtering if NaN values are retrieved over extreme treshold values
            ind=length(x)
            for i in 1:length(y)
                if isnan(y[i])==true
                    ind=i-1
                    break
                end
            end
            #Saving before extreme values (if present)
            x=x[1:ind]
            y=y[1:ind]
            #x=x[1:findlast(x->x>0, x)]
            #y=y[1:findlast(x->x>0, x)]
            #Computing AUPRC
            area = abs(sum((y[i] + y[i+1]) * (x[i+1] - x[i]) / 2 for i in 1:length(x)-1))
            push!(AA, area)
        end
        push!(AU, AA)
    end
    return AU
end


function mask_DRACH(seq, Dict_DRACH) 
    """Function which localize the position of only DRACH sites in a sequence. It takes as input:
        1) seq: sequence to mask 
        2) Dict_DRACH: Dictionary of DRACH motifs
        The function return: 
        1)The position of all DRACH motifs found in increasing order
    """
    positions_masked=[] 
    #Position of all found DRACH (collection of ranges from initial to end position of the DRACH)
    indexes=vcat([findall(a, seq, overlap=true) for  a in keys(Dict_DRACH)]...)
    #Filtering from the collection of DRACH ranges, to single positions (considering only the third components, which correspond to the A sites)
    ind=[indexes[i][3] for i in 1:length(indexes)]
    #Sorting in increasing order
    return sort(ind)
end

function score_ps_ds_pm_DRACH(M5, M6,set,genes, prior; verb=false, map=[], d_met=[], d_not=[], binning=0.001,N_window=8, VERB=false)
    """
    The function computes the score on a set of genes. It takes as input:
        1)M_met: The model learnt for methylated sequences
        2)M_not_met: The model learnt for non-methylated sequences
        3)set: a dictionary containint test set
        4)genes: the set of genes for a subset0 of the set dictionary for which score want to be computed
        5)prior: the prior for sequences
        Some additional infromation can be used
        6)verb: it is a boolean variable to include or not the distance m6A distribution, it is initialized as false
        7)d_met: the ditribution of m6A sites, it is a vector conaint the set of point [x,y], and it is initialized to zero
        8)d_not: the ditribution of non methylated A sites, it is a vector conaint the set of point [x,y], and it is initialized to zero
        9)binning: the size of the discretization for the distribution of m6A and A sites
        10)N_window: the size of short motives to score around adenosines. This is set by default equal to 8
        11)VERB: to save the stoichiometry (if available) of only the methylated DRACH motifs
    The function returns:
        1)score_met: the score assigned by the model to all DRACH m6A sites in the reference experiment
        2)score_not: the score assigned by the model to all DRACH non modified A sites in experiments
        3)r: the stoichiometry (if available) of only the methylated DRACH motifs
    """ 
    score_met, score_not_met, r=[],[],[]
    for label in genes #loop for all considered genes
        #Scoring with the decided model strategy test transcripts
        score=score_genome_complete(M5, M6, prior[label].perc_met/100, set[label].sequence, N_window; verbose=verb, Mapping=map, dist_met=d_met, dist_not=d_not, bin=binning, gene=label)
        #Finding the positions of all methylated sites inside DRACH motifs
        DRACH_met= intersect(mask_DRACH(set[label].sequence, Dict_DRACH), set[label].pos_met)
        #Saving for all the DRACH motifs the corresponding stoichiometry (iff the booelan variable "VERB" is "true)
        if VERB==true
            idx = findall(x -> x in DRACH_met, set[label].pos_met)
            push!(r, set[label].ratio[idx])
        end
        #Find all the non-methylated motifs inside DRACH motifs
        DRACH_not= intersect(mask_DRACH(set[label].sequence, Dict_DRACH), set[label].pos_not)
        #Add model scores of the position of the selected sites
        push!(score_met, score[DRACH_met])
        push!(score_not_met, score[DRACH_not])
        end
        #Return the methylated and non-methylated score and the stichiometry
        return score_met,score_not_met, r
end

function energy_tri(seq, h, J, K)
    """Function to compute the full energy of the model based on 1,2,3-mers statistics at each position of a given sequence of length N. It takes as input:
        1) A sequence of length n
        2) h: a vector for each nucleotide position of the sequence, containing a dictionary for all 1-mers 
        3) J: a vector for each di-nucleotide position of the sequence, containing a dictionary for all 2-mers 
        4) K: a vector for each tri-nucleotide position of the sequence, containing a dictionary for all 3-mers 
    The function returns:
        1) The energy of the model for the given sequence.
    """
    #Initialization
    e=0  
    #Summing all 1-mers contributions
    [e+=h[i]["$(seq[i])"] for i in 1:length(seq)]
    #Summing all 2-mers contributions
    [e+=J[i][seq[i:i+1]] for i in 1:length(seq)-1]
    #Summing all 3-mers contributions
    [e+=K[i][seq[i:i+2]] for i in 1:length(seq)-2]
    return e 
end

function score_DRACH(M1, M2, Dict_DRACH, prior)
    """
    Function which maginalize model probability over 4-mers to compute the absolute score given to DRACH motifs. It takes as inputs:
    1)M1: methylated model
    2)M1: non methylated model
    3)Dict_DRACH: dictionary containing as keys all the DRACH motifs
    4)prior: constant prior value to score DRACH motifs
    The function returns: 
    1)Dict_score_D : A dictionary containg as keys all the DRACH motifs, and as values the score assigned by the model
    """
    Dict_score_D=Dict()
    for a in keys(Dict_DRACH) #Considering all DRACH motifs
        #Initialization
        A1, A2=0,0
        for b in mots_list4 #Marginalization computation over extremes of 8-mers
            A1+=exp(energy_tri(b[1:2]*a[1:2]*a[4:5]*b[3:4], M1.h, M1.J, M1.K))/exp(M1.logZ)
            A2+=exp(energy_tri(b[1:2]*a[1:2]*a[4:5]*b[3:4], M2.h, M2.J, M2.K))/exp(M2.logZ)
        end
        #Methylated score
        m=A1*prior
        #non methylated score
        u=A2*(1-prior)
        #Normalization constant
        s=m+u
        #Final methylation probability
        M=m/s
        Dict_score_D[a[1:2]*a[4:5]]=M 
    end
    return Dict_score_D 
end

function score_8_mers(M5, M6, prior)
    """Function to score all 8-mers. It takes as entry: 
        1)M5: methylated model
        2)M6: non methylated model
        3)prior: constant prior value to score DRACH motifs
        The function returns: 
        1)score : the score assigned by the model for each motif
    """
    #Initialization of an empty vector
    score=[]
    for seq in mots_list8 
        A1=exp(energy_tri(seq, M5.h, M5.J, M5.K))/exp(M5.logZ)
        A2=exp(energy_tri(seq, M6.h, M6.J, M6.K))/exp(M6.logZ)
        #Methylated score
        m=A1*prior
        #non methylated score
        u=A2*(1-prior)
        #Normalization constant
        s=m+u
        #Final methylation probability
        M=m/s
        push!(score, M)
    end
    return score
end

function filtering_non_nucelotides(test)
    """
    Function to find amongst the transcripts if there is one containing non-nucleotide charcters.
    A message is shown in this case with the details, and the gene is saved to be later excluded from the testing dataset.
    It takes as entry: 
    1) test: A dictionary containing the sequences of the genes
    The function returns:
    1)genes_not: a list of all the genes with non-nucleotide letters
    """
    #Initialization to empty vector
    genes_not=[] 
    #loop over all the genes (which are the keys of the vocabulary 'test')
    for a in keys(test)
        #looking for all the positions which do not encode for a nucleotide letter
        pos=findall(x->x!='A'&& x!='C' && x!= 'G' && x!='T', test[a].sequence)
        #looking if at least one position is found (so if the vector 'pos' is not empty)
        if isempty(pos)==false  
            #In this case a message with the details are shown
            println("The gene $a contain a (some) non-nucleotide(s) character(s) at position(s) $pos corresponding to the chacracter(s): $(test[a].sequence[pos]) ")
            #This genes is saved in the list of genes to be talen out of the test 
            push!(genes_not, a)
        end
    end
    return genes_not 
end

function filter_DRACH(set; verbose=false, ratio=[])
    """
    Function to find all DRACH motifs. It takes as entry:
    1) set: containing the list of 8-long motifs to filter (withouth the central adenosine) 
    A set of optional parameter:
    2) vebrose: a boolena variable (the default value is 'false')
    3)ratio: an empty vector to save the toichiometry of DRACH motifs (the default value is an empty vector)
    The function returns:
    1)set_D: if verbose is false, consisitng on the set of DRACH motifs
    2)R: a vector of stoichiometry for the corresponding DRACH motifs (iff verbose is 'true')
    """
    #initialization of the empty vector
    set_D=[]
    #initialization of the a dictionary of DRACH motifs without central adenosine 
    Dict_DRCH=Dict()
    [Dict_DRCH[a[1:2]*a[4:5]]=[]   for a in DRACH]
    #Verying if the central motif is a DRACH, and savining the motif if this is the case
    [(haskey(Dict_DRCH, seqs[3:6]) ? push!(set_D, seqs) : nothing) for seqs in set]
    #Looking at the variable to save the stoichiometry
    if verbose==true
    #Verying if the central motif is a DRACH, and savining the stoichiometry
        R=[]
        [(haskey(Dict_DRCH, set[i][3:6]) ? push!(R, ratio[i]) : nothing) for i in 1:length(set)]
        return set_D, R 
    else
    return set_D
    end
end

function null_model(seqs)
    """Function to define the Baseline DRACH model. It assigns a probability of 1 to DRACH motifs, 0 otherwise. 
    It takes as entry:
    1) A sequence of 8-mers without central adenosine
    It returns:
    2) The probability of the baseline model 
    """
    result=[]
    #Constructing a Dictionary of DRACh motifs without the central adenosine
    Dict_DRCH=Dict()
    [Dict_DRCH[a[1:2]*a[4:5]]=[]   for a in DRACH]
    #Veryfing the condition and assigning the porbability
    (haskey(Dict_DRCH, seqs[3:6]) ? result=1 : result=0)
    return result
end

function copy_common_only(common_genes, T, T_ref)
    """Function to add gene which were not found in altered datasets, which were present in unaltered (control) condition over the same cell line,
     with at least one detected methylated motif. The function takes as input:
        1) common_genes: a set of genes name which are different form train, and validation set, common beween prior dataset and test set.
        2) T: test set in altered condition (a dictionary containing testing genes as keys)
        3) T_ref: test set in unaltered condition (a dictionary containing testing genes as keys)
     The function returns: 
        1) V: the resulting test set in unaltered conditions, with added genes, present in unaltered conditions, with all sites indicated as non-mehtylated.
     """
     #initialization of the resulting dictionary
    V = Dict()
    for a in common_genes #looping over all the set of genes from unaltered conditions
        if haskey(T,a)==true #if the gene is present in the altered conditions too, than I simply copy it into "V'
            V[a]=deepcopy(T[a])
        else #otherwise, I copy the information from the unaltered conditions
            tt=deepcopy(T_ref[a])
            #redefinition of the positions, with none methylated (since it was not present in the altered condition)
            tt.pos_met=[] 
            #redefinition of all sites as non-methylated 
            tt.pos_not=findall('A', tt.sequence)
            #updating the percentage of methylation to zeros
            tt.perc_met=0 
            #copying the gene with the adapted information for the altered dataset
            V[a]=tt 
        end
    end
    return V
end

function divide_score(S, keys, set)
    """Function to divide scoring in methylated and non-methylated sites. It takes as arguments:
        1) S: the score of genes 
        2) keys: the names of the genes 
        3) set: the information of the detected m6A sites for all the genes
    The function return: 
        1) the score of methylated and non-methylated motifs for each gene 
    """
    score_met,score_not_met=[],[]   #Initialization of empty scores
    for i in 1:length(keys)
        push!(score_met, S[i][set[keys[i]].pos_met]) #Selecting score of m6A sites
        push!(score_not_met, S[i][set[keys[i]].pos_not]) #Selecting score of A sites
    end
    return score_met, score_not_met 
end

function divide_score_DRACH(S, keys, set)
    """Function to divide scoring in methylated and non-methylated DRACH sites. It takes as arguments:
        1) S: the score of genes 
        2) keys: the names of the genes 
        3) set: the information of the detected m6A sites for all the genes
    The function return: 
        1) the score of methylated and non-methylated motifs for each gene 
    """
    score_met,score_not_met=[],[]   #Initialization of empty scores
    for i in 1:length(keys)
        #looking at methylated DRACH sites
        DRACH_met= intersect(mask_DRACH(set[keys[i]].sequence, Dict_DRACH), set[keys[i]].pos_met)
        #Find all the non-methylated motifs inside DRACH motifs
        DRACH_not= intersect(mask_DRACH(set[keys[i]].sequence, Dict_DRACH), set[keys[i]].pos_not)
        #Add model scores of the position of the selected sites
        push!(score_met, S[i][DRACH_met]) #Selecting score of m6A sites
        push!(score_not_met, S[i][DRACH_not]) #Selecting score of A sites
    end
    return score_met, score_not_met 
end


function score_ava(score_met,score_not_met, ind; verb=false, pass=0.001, rd=[])
    """Function to compute the average quantities over a a subset of the overall genes. 
    The function takes as input: 
        1) score_met: the prediction of methylation for all methylated sites 
        2) score_not_met: the prediction of methylation for all non-methylated sites 
        3) ind: a set of inxdexes to divide data in different indpendent subsets
    Some additional parameters can be pased to the function: 
        4) verb: a boolean variable to control the number of wanted saved information  (default value is false)
        5) pass: the value of the step for the precision recall curve 
        6) rd: a set of indexes of chosen of non-methylated sites to allow fair comaprison at fixed 90% non-methylated adenosine. The size corresonds to the number of sub-partitions corresponding to: length(ind))
    The function returns:
        1) imb: the overall ratio between methylated and non-methylated sites for a given subset
        2) t: the treshold value at the maximum F1 value of the PRC for each dataset
        3) a: the AUPRCs 
        4) f: the maximum F1 scores

    """
    #initialization to empty vectors tresholds, AUPRC, F1 maximum score and dataset imbalance
    t,a,f, imb=[],[],[], []
    c=1 #counting index
    for i in ind #looping over the number of independent datasets
        if verb==false #if verbose is fals, the imbalance is taken as it is, otherwise, only 10% o m6A sites
        _,tre, A, F=plot_PR_score(vcat(score_met[i]...), vcat(score_not_met[i]...),  "black", " "; verbose=true, pas=pass)
        else
            if isempty(rd)==true #looking if it has been defined a chosen set on non-methylated adenosine
            _,tre, A, F=plot_PR_score(vcat(score_met[i]...), shuffle!(vcat(score_not_met[i]...))[1:9*length(vcat(score_met[i]...))],  "black", " "; verbose=true, pas=pass)
            else 
                _,tre, A, F=plot_PR_score(vcat(score_met[i]...), vcat(score_not_met[i]...)[rd[c]],  "black", " "; verbose=true, pas=pass)
                c+=1 #If the the set of non-methylated adenosine the coundint index is updated of one.
            end
        end
        #adding the computed treshold, AUPRC, max F1 score and imbalance
        push!(t, tre)
        push!(a, A)
        push!(f, F)
        push!(imb, length(vcat(score_met[i]...))/(length(vcat(score_met[i]...))+length(vcat(score_not_met[i]...))))
    end
    #returning the wanted quantities
    return imb, t, a, f 
end

function ava_dist_DRACH(test, labels, S_D, ind, sous)
    """ The function comupute the frequency of DRACH methylated motis in a given test set, ordering in an decresing ordered based on model score.
    It takes as input: 
        1) test: a given test dataset with detected m6A information overa a group of genes
        2) labels: the group of testes genes
        3) S_D: dictionary contianing as keys the DRACH motifs without central adenosine, and as values, the corresponding model prediction of methylation probability. 
        4) ind: the index of DRACH motifs orderd in decreasing order by model prediction
        5) sous: the indexes of the independent sous-set of testing geens to compute average information
    The function returns: 
        1) St: the freqeuncy ordered by models prediction (in decreasing order) for each independent testing-subset of genes
    """
    ST=[] #initialixation of the empty vector
    for i in 1:length(sous) #looping over the independent sub-set of testing genes
        #creating for testing genes all methylated and non -methylated 8-mers
        test_motives, _, _,test_motives_not, _=crete_set_all(test, vcat(keys(labels)...)[sous[i]]); 
        #Filtering 8-mers over the ones which include DRACH motifs
        test_motives_D,_=filter_DRACH(test_motives; verbose=true, ratio=vcat(ones(1, length(test_motives))...))
        #Computing frequency over the set of filtered 8-mers in the given subset and adding it to the resulting 'ST' vector
        push!(ST,freq_met_motifs(S_D,test_motives_D, ind))
    end
    return ST 
end
function freq_met_motifs(S_D,test_D, ind)
    """ The function comupute the frequency of DRACH methylated motis in a given test set, ordering in an decresing ordered based on model score.
    It takes as input: 
        1) S_D: dictionary contianing as keys the DRACH motifs without central adenosine, and as values, the corresponding model prediction of methylation probability. 
        2) test_D: all 8-mers containing DRACH motifs in a given test set 
        3) ind: the index of DRACH motifs orderd in decreasing order by model prediction
    The function returns: 
        1) freq: the freqeuncy ordered by models prediction (in decreasing order)
    """
    #initialization of an empty dictionary for DRACH motifs
    Dict_mers=Dict()
    [Dict_mers[a]=0 for a in keys(S_D)] #initialization to zero
    L=length(test_D) #length of all motifs
    [Dict_mers[test_D[i][3:6]]+=1/L for i in 1:L] #counting the total frequency of a motif in the methylated datase (over DRACH motifs)
    #[Dict_mers[test_D[i][3:6]]+=1 for i in 1:L]
    freq=vcat(values(Dict_mers)...)[ind] #ordering frequency based on predicted probability by the model
    return freq
end

function ava_avarege_distribution(test,labels,av, bb, sous)
    """Function to compute the localization distribution of methylated DRACH motifs average over an independent subpartition of genes. It takes as entry:
    1) test: a dictionary containint test set with methylation information
    2) labels: keys of the dictionary of the considered testing transcripts 
    3) av: the window size to smooth the resulting distribution 
    4) bb: the binning size to compute the histogram, therefore the density profile, of sites
    5) sous: the indexes of the independent sous-set of testing geens to compute average information
The function returns: 
    1) xmean: the spaced bins over which the density distrubution has been compute_F1_threshold for each subsett
    2) ymean: the avareged density distribution over xmean of DRACH motifs for each subset
"""
    #initialization of empty vector
    X,Y=[],[]
    for i in sous #looping over the sub-partitions of genes
        x, y=avarege_distribution_DRACH(test, vcat(keys(labels)...)[i], av, bb; label=1) #computing the distribution of methylated DRACH
        push!(X, x)
        push!(Y, y)
    end
    return X, Y 
end

function stat_m6A(MAP,Test, tre, S, sous)
    """Function to compute some statistics relatively to a given test set, averaged over sub-partitions. The function takes as entry:
        1) MAP: Dictionary containint the mappin towords an ideal sequence made of 33.3̄% 5'-UTR, 33.3̄% CDS and 33.3̄% 3'-UTR
        2) Test: a dictionary containint test set with methylation information
        3) tre: average treshold at maximum F1 score dictionary 
        4) S: a dictionary with genes as keys, and predicted probabilities at each site as values  
        5) sous: the indexes of the independent sous-set of testing geens to compute average information
    The function returns;
        1) R: the percentage of DRACH methylated sites for each  sub-set
        2) P: percentage of predicted DRACH methylated sites for each  sub-set
        3) D: the distribution of predicted DRACH methylated sites for each  sub-set
    """
    #Initialization of empty vectors
    R, P, D=[],[],[]
    for i in 1:length(sous) #looping over the sub-partitions of genes
        pp,p_ref,dist_met_pred_ps_ds_pm=[],[],[]
        for label in vcat(keys(MAP)...)[sous[i]] #looping over all the genes for the given subpartition
            mm=mask_DRACH(Test[label].sequence, Dict_DRACH) #finding all DRACH site
            v=intersect(findall(x->x>=tre[i], S[label]), mm) #finding all predicted DRACH methylated sites
            push!(pp,length(v)/length(mm)*100) #percentage of predicted DRACH methylated sites over all DRACH motifs
            push!(dist_met_pred_ps_ds_pm, MAP[label][v]) #mapping the position of predicted DRACH methylated sites into a metage prfile of 33.3̄% 5'-UTR, 33.3̄% CDS and 33.3̄% 3'-UTR
            int=length(intersect(Test[label].pos_met,mm )) #counting all the number of experimentally detected methylated sites
            (isempty(int)==false ? push!(p_ref, int/length(mm)*100) : push!(p_ref, 0)) #computing the percentage of DRACH methylated sites
        end
        push!(R,p_ref )
        push!(P,pp )
        push!(D, dist_met_pred_ps_ds_pm)
    end
    return R, P,D
end

function average_density_difference(prediction, experiment, N, binning)
    """Function to compute the  density difference (in independent sub-partitions of genes) between the amount of methylated sites between model predictions and experimental maps. 
    The function takes as imput:
        1) prediction: the amount of predicted methylated sites over a given treshold 
        2) experient: the amount of detected methylated sites from expeirments 
        3) N: the number of sub-partitions of genes to compute the average difference 
        4) binning: the bins size of the distribution to compute differences
    The function returns: 
        1) denisty_differe: the derired density difference quantity 
    """
    density_difference=[] #initialiation of an empty vector 
        for i in 1:N #looping over the sub-partitions of genes
            d, _=hist(abs.(prediction[i].-experiment[i]), bins=binning,density=true) #computing the density differnece
            push!(density_difference, d)
        end
    return density_difference
end


function average_dist_local(distance, N, binning)
    """Function to compute the distribution of methylated sites in the trascript (in independent sub-partitions of genes).
    The function takes as imput:
        1) distance: the position of methylated sites in a metage profile made of 33.3̄% 5'-UTR, 33.3̄% CDS and 33.3̄% 3'-UTR
        2) N: the number of sub-partitions of genes to compute the average difference 
        3) binning: the bins size of the distribution
    The function returns: 
        1) denisty: the derired density quantity 
    """
    density=[] #initialiation of an empty vector 
        for i in 1:N #looping over the sub-partitions of genes
            d, _=hist(vcat(distance[i]...), bins=binning,density=true)
            push!(density, d)
        end
    return density 
    end

    function plot_local(dd, y_ref, x_ref, col; y_max=7.1)
        """Function to plot the the distribution of methylated sites, mapped over a metage profile made of 33.3̄% 5'-UTR, 33.3̄% CDS and 33.3̄% 3'-UTR, between prediction and experimental results. 
        It takes as input: 
            1) dd: the predicted methylated sites
            2) y_ref: the experimental methylated sites profile
            3) x_ref: the binning (equivalent between experiments and predictions)
            4) col: the colour desired to plot the rersulting curve
            An optional argyent: 
            5) y_max: parameter to regulate the size of the y axis
        the functin retuns: 
            1) the desired plot with difference frequency distribution.
        """
        #Computing the mean and the standard deviation over subpartitions, and averaging over a fix window of 20 results 
        M=movmean(mean.(eachrow(( hcat(dd...)))),20)
        σ=movmean(std.(eachrow(( hcat(dd...)))),20)
        #Plotting mean curve and the standard deviation filling with the chosen colour
        plot(x_ref[1], M, color=col,label="Predicted", linewidth=1, ls="--") 
        fill_between(x_ref[1], M.-σ, M .+ σ, color=col, alpha=0.3)
        #Repeating the same as before for experimental results
        M=movmean(mean.(eachrow(( hcat(y_ref...)))),20)
        σ=movmean(std.(eachrow(( hcat(y_ref...)))),20)
        plot(x_ref[1], M, color=col,label="Experimental", linewidth=1) 
        fill_between(x_ref[1], M.-σ, M .+ σ, color=col, alpha=0.5)
        #settings
        xticks(fontsize=12)
        yticks(fontsize=14)
        legend()
        axvline(x=1/3, color="grey", ls="--", alpha=0.3 , zorder=0)
        axvline(x=2/3, color="grey", ls="--" ,alpha=0.3, zorder=0)
        ylabel(L"density", size=15)
        xlabel(L"position", size=15)
        ylim(0,y_max)
        return gcf()
    end

    function stat_plot(MOD, REF, lab, sizes, col)
        """Function to plot the methylated DRACH motifs requency different between prediction and experimental sites, ordered by model score (in decreasing order). The function takes as input: 
            1) MOD: the frequency of precticted methylated sites
            2) REF: the frequency of experimental methylated sites
            3) lab: label name for the altered condition
            4) sizes: window size 
            5) col: colout o the plotted elements 
        The function returns: 
            1) the plot of the difference methylated frequency
        """
        close("all")
        figure(figsize=sizes)
        #computing and plotting the difference, by averaging over independent sub-partitions of testing genes
        bar(vcat(1:N...), mean.(eachrow(( hcat(MOD...).-hcat(REF...)))),yerr=std.(eachrow(( hcat(MOD...).-hcat(REF...)))), edgecolor="black", linewidth=0.5,capsize=2,color=col, ecolor=col,alpha=0.5, label=lab)
        ylim(-0.06, 0.17)
        legend(fontsize=12)
        title(L"DRACH \ frequency", size=20)
        xticks(vcat(vcat(1:N...)),vcat(MOT5[1:N]), rotation=90)
        ylabel(L"\Delta f_{methylated}", size=15)
        xlabel(L"motifs", size=15)
        xticks(fontsize=12)
        yticks(fontsize=14)

        return gcf()
    end

    function plot_dens_diff(binning, col, lab, dd; y_max=0.8)
        """ Function to plot the difference amount of methylated sites (between predictions and experiments). 
        It takes as input: 
        1) binning: the bins size of the distribution
        2) col: colout o the plotted elements 
        3) lab: label name for the altered condition
        4) dd: the density difference of the amount predicted methkated sites and the experimental one
    and an optional argument:
         5) y_max: parameter to regulate the size of the y axis'
    The functinn returns: 
        1) the corresponding plot with the density difference
        """
        bar(binning[1:end-1].+0.5, mean.(eachrow(( hcat(dd...)))), yerr= std.(eachrow(( hcat(dd...)))), width=1,color=col, ecolor=col, alpha=0.5, label=lab)
        xlabel("\$| \\Delta_\\% DRACH\$ | \n (prediction-experiment)", size=12)
        ylabel(L"density", size=15)
        #ylim(-0.1,1)
        xticks(vcat(0:10:30...), [ "0", "10", "20", "30"])
        xlim(-0.5,30)
        ylim(-0.05,y_max )
        xticks(fontsize=12)
        yticks(fontsize=14)
        legend(fontsize=12)
        title(L"DRACH \ abundance",size=20)
        return gcf()
    end

    function create_seq_dict(seqs)
        """Function to create a dictionary containing the name of the genes as keys, and the sequences as values.
        The function takes as input: 
            1) seqs: the set of sequence 
            The function returns:
            1) Dicit: the resulting dictionary
        """
        #Initialization of empty dictionary
        Dicti=Dict()
        [Dicti[seqs[i][1]]=seqs[i][2] for i in 1:length(seqs)]
        return Dicti 
    end

    function format_Deep(D, genes, Dict_seq)
        """Function which provides a dictionary, containing as keys the selected genes, and as value a tuple containing as first entry
        a vector with the position of the assigned methylation probability, and as second entry, the vector of the scores. 
        It takes as input: 
            1) D: the output of DeepSRAMP predictions, where probability are assigned to a sequence aligned to a reference
            2) genes: a list of selcted genes
            3) Dict_seq: a dictionary containing the  reference seq to which score were assigned
        """
        A=Dict()
        for a in genes #looping over the list of genes
            poss=findall(x->x==a, D[:, 1]) #finding the position of sites to whcih score were assigned for the given gene
            map_pos=[count('-', Dict_seq[a][1:D[p,2].+1]) for p in poss] #mapping position from the aligned sequence to the reference one without any gaps
            A[a]=(D[poss,2].+1 .- map_pos, D[poss, 3]) #Saving the new position and the score assigned
        end
        return A
    end

    function scoring_Deep_met(common_genes,Dict_Deep, Test_HEK)
        """Function to score methylated sites from DeepSRAMP results.
        The function takes as input:
            1) common_genes: the set of genes to format scores 
            2) Dict_Deep: a dictionary containing the  reference seq to which score were assigned by DeepSRAMP
            3) Test_HEK: the corresponding test set with detected methylated sites. 
        The function returns: 
            1)score_met_Deep_D: the score of experimental methylated DRACH motifs assigned by DeepSRAMP
        """
        score_met_Deep_D=[] #initialization of the empty vector for methylated siets
        for a in common_genes #looping over the set of genes
            pos_mask=mask_DRACH(Test_HEK[a].sequence, Dict_DRACH) #finding all DRACH motifs
            S=[]
            DRACH_met=intersect(pos_mask, Test_HEK[a].pos_met) #Finding all experimentally detected methylated DRACH motifs
            N=length(DRACH_met) #saving the amount of DRACH methylated motifs
            t=intersect(Dict_Deep[a][1], DRACH_met) #Finding the position of the DRACH methylated sites which were assinged a score by DeepSRAMP
            #For each DRACH methylated sites, for which it was assigned a score to DeepSRAMP, i find the position of the site in the list
            #and the corresponding probability is saved 
            push!(S, [Dict_Deep[a][2][findfirst(x->x==i, Dict_Deep[a][1])] for i in t]) 
            for i in 1:N-length(vcat(S...)) #If any DRACH methylated sites were not assigned a score by DeepSRAMP, than a score of 0 is assigned
                push!(S, 0)
            end
            push!(score_met_Deep_D, vcat(S...))

        end
        return score_met_Deep_D
    end
    
    function scoring_Deep_not(common_genes,Dict_Deep, Test_HEK)
        """Function to score non-methylated sites from DeepSRAMP results.
        The function takes as input:
            1) common_genes: the set of genes to format scores 
            2) Dict_Deep: a dictionary containing the  reference seq to which score were assigned by DeepSRAMP
            3) Test_HEK: the corresponding test set with detected methylated sites and non-methylated ones. 
        The function returns: 
            1)score_met_Deep_D: the score of experimental non-methylated DRACH motifs assigned by DeepSRAMP
        """
        score_not_met_Deep_D=[] #initialization of the empty vector for non-methylated siets
        for a in common_genes #looping over the set of genes
            pos_mask=mask_DRACH(Test_HEK[a].sequence, Dict_DRACH)  #finding all DRACH motifs
            S=[]
            DRACH_not_met=intersect( pos_mask, Test_HEK[a].pos_not) #Finding all experimentally detected non-methylated DRACH motifs
            N=length(DRACH_not_met) #saving the amount of DRACH non-methylated motifs
            t=intersect(Dict_Deep[a][1], DRACH_not_met)
            #For each DRACH methylated sites, for which it was assigned a score to DeepSRAMP, I find the position of the site in the list
            #and the corresponding probability is saved 
            push!(S, [Dict_Deep[a][2][findfirst(x->x==i, Dict_Deep[a][1])] for i in t])
            #If any DRACH non-methylated sites were not assigned a score by DeepSRAMP, than a score of 0 is assigned
            push!(S, vcat(zeros(1, N-length(vcat(S...)))...))
            push!(score_not_met_Deep_D, vcat(S...))
        end
        return score_not_met_Deep_D
    end

    function methylated_mots(Dict_DRACH, MAP_mRNA_GLORI, test)
        """Function to give the list of all DRACH methylated motifs for each testing sequence and the corresponding stoichiometry information (concerning only GLORI-seq data).
        The function takes as input: 
            1) Dict_DRACH: a dictionary having as keys all DRACH motifs 
            2) MAP_mRNA_GLORI: a dictionary containing all testing coding mRNAs
            3) test: a dictionary containint test set with methylation information
        The function returns:
            1) F_mots: list of all DRACH methylated motifs for each testing sequence
            2) r: stoichiometry information for each DRACH methylated motif in each testing sequence.
        """
        #initialization empty dictionary and list
        r=Dict()
        F_mots=[]
        [r[a]=[] for a in keys(Dict_DRACH)] #Adding all keys corresponding to DRACH motifs
        for label in keys(MAP_mRNA_GLORI) #looping over all testing coding mRNA
            #looking at the positon of all methylated DRACH site in the given sequence
            DRACH_met= intersect(mask_DRACH(test[label].sequence, Dict_DRACH), test[label].pos_met) 
            #finding the index of the position of the DRACH methylated motif inside the list of all methylated motifs
            idx = findall(x -> x in DRACH_met, test[label].pos_met)
            temp=[]
            #Saving the DRACH methylated motif for each found site
            [push!(temp, test[label].sequence[test[label].pos_met[i]-2:test[label].pos_met[i]+2]) for i in idx]
            #Pushing the corresponding stoichiometry
            [push!(r[test[label].sequence[test[label].pos_met[i]-2:test[label].pos_met[i]+2]], test[label].ratio[i]) for i in idx]
            push!(F_mots,temp)
        end
        return F_mots, r 
    end

    function stochiometry(F_mots, set)
        """Function to compute the predicted soichiometry for a given model.
        The function takes as input: 
            1) F_mots: list of all DRACH methylated motifs for each testing sequence
            2) set: the score of all experimentally methylated sites from a given model
        The function returns:
            1) stoc: the predicted soichiometry
        """
        #Iniialization of an empty dictionary
        stoc=Dict()
        [stoc[a]=[] for a in keys(Dict_DRACH)] #Adding all DRACH motifs as keys
        for i in 1:length(F_mots) #looping over each testing sequence
            if isempty(F_mots[i])==false #looking at sequence which have at least one methylated DRACH
                for j in 1:length(F_mots[i]) #looping over all DRACH methylated sites of the sequence
                    push!(stoc[F_mots[i][j]],set[i][j])
                end
            end
        end
        return stoc
    end