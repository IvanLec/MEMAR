#Structure to save detected methylation sites over IAV segments
mutable struct m6A
    segment::Int64 #Segment name 
    strand::String #Positive or Negative
    position::Int64 #Position in the given segment of the detected m6A
    sequence::String #Nucleotide sequence of the segment
 end 
#Structure to collect additional information relative to  detected methylation sites over IAV segments
struct info_IAV
    sequence #Nucelotide sequence
    cds #Position of codon START and STOP
    pos_met #postion of the detected m6A sites
    pos_not #position of all non-detected m6A sites
    perc_met #percentage of methylation (between m6A and A sites)
end

 function create_m6A(indexes::Vector, rows, i::Int64) 
    """Function to create a dataset of m6A sites from the file containing experimental results for IAV. The function takes as input:
        1) indexes: is a vector containing the rows with the information ordered following the structure "m6A";
        2) rows: containing the file rows with experimental m6A information 
        3) i: the selected row to create the structure
    The function returns:
        1) m6A: the structure relativly to all the infroamtion contained in row 'i'
    """
    #Since we have define the Alphabet with DNA letters, withouth any loss of generalization, nucleotide 'U' will be transformed into 'T' 
    return m6A(rows[i,indexes[1]], rows[i,indexes[2]], rows[i,indexes[3]], replace(rows[i,indexes[4]], "U"=>"T")) 
end

function find_string(String, list; th=1)
    """Function to map the sub-string from the original experimental seqeunces to one with highest resemblence. Position of methylated motifs are therefore mapped into the new position.
    The function takes as arguments: 
        1) String: the sequence to perfom the map;
        2) list: the m6A structures with all the information relativly to the methylated motifs with central adenosine and their position in the original experimental sequence;
    the function can take an optional parameter:
        3) th: which corresponds to the number of accepted mutations between the substring in the original sequence and the mapped one, its default value is 1.
    The function returns:
        1) D: a dictionary containing all the infromation of the mapped methylated motifs relativly to the mapped sequence.
    """
    D=Dict() #Looking for all 'A' sites 
    pos_a=findall('A', String)
    #Since submotifs are of length 13, we filter adenosines retaining only those which can be central to these motifs length
    pos_A=filter(x->x-6>=1 && x+6<=length(String), pos_a)
    for i in 1:length(list) #looping over the detected methylated motifs in 'list'
            for j in pos_A #looping over 'A' positions
                C=count(x -> x[1] != x[2], zip(list[i].sequence, String[j-6:j+6])) #counting difference of submotifs of length 13 between original and reference sequence
                if C<=th #only iff the nucleotide mutation are lower than the tresholds value we consider to have found the submotifs (otherwise we continue with the next position)
                   if C==0 #if there is a full (100%) identity
                        D[list[i]]=(j, String[j-6:j+6],  C) #then in the dictionary this is added , giving the new position, the substring and the 0 nucleotide difference
                   else
                    #otherwise we look for the nucleotide mutation, saving its position
                    diff_position =[k for k in 1:13 if String[j-6:j+6][k] != list[i].sequence[k]][1]
                    D[list[i]]=(j,  String[j-6:j+6], C,  "$(String[j-6:j+6][diff_position]) -> $(list[i].sequence[diff_position])", diff_position)
                    #=we therefore save in the dictionary this is added this methylated motif, providing the new positon, the new substring, the nucleotide difference
                    and the mutationt from the reference to the original sequence, and its relative position in the submotif=#
                   end
                end
            end
        end
    return D 
end

function cmp_inv_IAV(Duck)
    """Function to consider the vRNA (-sstr) of IAV, which is the complementary inverse of the mRNA (+sstr). The function takes as aargument:
        1)Duck: the set of sequneces (with their labels) of the reference sequence. 
        The function returns: 
        1) A:  the set of inverted and ocmplemented sequneces (with their labels) of the reference sequence. 
    """
    A=[] #empty initialization
    [push!(A, (Duck[i][1], inverse_complement(Duck[i][2]))) for i in 1:length(Duck)]
    return A 
end

function inverse_complement(seq::String)
    """Function to inverse complement a sequence. It takes as input:
        1) seq: the nucleotide sequence (in DNA alphabet and uppercase)
    The function returns:
        1) inv_c: the inverse complement of the input sequence.
    """ 
    #Defyning the dictionary of complemented nucleotides:
    comp = Dict('A'=>'T', 'T'=>'A', 'C'=>'G', 'G'=>'C')
    #reversing the sequence
    rev = reverse(seq)
    #complementing the reverse sequence for each position
    inv_c = join([comp[i] for i in rev])
    return inv_c
end

function AUC_IAV(score, pos,seq, TRE; Dict=Dict(),l=0, pass=0.0001, c, legen)
    """Function to plot the PRC over infleuza A and the F1 at a given treshold value.
    The function takes as input: 
        1)score: a vector contaiing for each position the porbability to methylate given by a model
        2)pos: the position of all the detected methylated sites
        3)seq: the concerned sequence
        4) TRE: the treshold value to compute the F1 score
    Some additional parameters can be provided:
        5) Dict: a dictionary contianing all the submotifs to filter out results (default is a empty dictionary)
        6) l: the window size from a central adenosine to filter the motifs (for DRACH is 2, for RAC is 1). The default value is 0
        7) pass: the step to change the treshold for the plot of the precision recall curve 
        8) c: a string containing the coulor
        9) legen: the legend to be shown for the plot
    The function returns:
        1) the plot of the Precision Recall with maximum F1 score and the F1 score at the given treshold value 'TRE'
    """
    #Initialization of a vector to contain a flag=1 for methylated sites
    mask=zeros(length(seq))
    [mask[i]=1.0 for i in pos]
    #Looking for all possible A sites
    pos_A=findall(x->x=='A', seq) 
    if l==0 #if we do not need to filter out motifs
        plot_PR_score(score[intersect(findall(x->x==1, mask), pos_A)], score[intersect(findall(x->x==0, mask), pos_A)], c, legen, pas=pass)
        #Computing the F1 value at the given treshold
        #recall
        X=length(intersect(findall(x->x>=TRE,score), pos_A,findall(x->x==1, mask) ))/length(intersect(findall(x->x==1, mask), pos_A))
        #precision
        Y=length(intersect(findall(x->x>=TRE,score), pos_A,findall(x->x==1, mask) ))/length(intersect(findall(x->x>=TRE,score), pos_A))
        println("F1 at treshold: ($(round(X, digits=3)), $(round(Y, digits=3)))")
        println("Imbalance: $(round((length(intersect(findall(x->x==1, mask), pos_A))/length(pos_A)*100), digits=3))%")
    else #if we need to filter out over some motifs
        #looking for all the position of the concerned motifs in the dictionary
        range_Dict=vcat([findall(a, seq) for a in keys(Dict)]...)
        pos_Dict=[range_Dict[i][l+1] for i in 1:length(range_Dict)]
        plot_PR_score(score[intersect(findall(x->x==1, mask), pos_Dict)], score[intersect(findall(x->x==0, mask), pos_Dict)],c, legen, pas=pass)
        #recall
        X=length(intersect(findall(x->x>=TRE,score), pos_Dict,findall(x->x==1, mask) ))/length(intersect(findall(x->x==1, mask), pos_Dict))
        #precision
        Y=length(intersect(findall(x->x>=TRE,score), pos_Dict,findall(x->x==1, mask) ))/length(intersect(findall(x->x>=TRE,score), pos_Dict))
        println("F1 at treshold: ($(round(X, digits=3)), $(round(Y, digits=3)))")
        println("Imbalance: $(round((length(intersect(findall(x->x==1, mask), pos_Dict))/length(pos_Dict)*100), digits=3))%")
    end
    #Plotting the maximum F1 score
    scatter(X, Y, marker="D", edgecolor="black", linewidth=0.5, color="firebrick", zorder=2)
    text(X+0.03, Y-0.03, "F1 at chosen \n treshold:$(round(digits=3, (2 .* (Y.*X)./(Y.+X))))", fontsize=12, color="black", zorder=2)
    xticks(fontsize=14)
    yticks(fontsize=14)
    legend(fontsize=12)
    ylabel(L"precision", size=15)
    xlabel(L"recall", size=15)
    ylim(-0.05,1.05)
    xlim(-0.05,1.05)
    return gcf()
end 


function filter_score(score, seq, ww, Dict_filtr)
    """Function to filter the score only over selected motifs. It takes as input: 
        1)score: a vector contaiing for each position the porbability to methylate given by a model
        2)seq: the concerned sequence
        3) ww: the window size from a central adenosine to filter the motifs (for DRACH is 2, for RAC is 1)
        4) Dict_filtr: a dictionary contianing all the submotifs to filter out results 
    The function returns: 
        1) the score with 0 to the positions of all motifs which are not contained in Dict_filtr
        2) filtr: the filtering vector
    """
    #Finding all A positions
    poss=findall('A', seq)
    #initialization the filter to ones
    filtr=ones(length(seq)...)
    for i in 1:length(poss)#looping oer all the position
        if poss[i]>ww && poss[i]<length(seq)-(ww+1) #Checking a windowd of 2ww+1 can be centered at this position
            (haskey(Dict_filtr, seq[poss[i]-ww:poss[i]+ww])==false ?  filtr[poss[i]]=0 : nothing) #Checking if the motif does not belong to the dictionary and putting a 0
        end 
    end
    return score.*filtr, filtr
end

function plot_score(score, seq, pos_met; Dict_filtr=[], lw=1, ww=0, sizes=(5,5), tresh=0)
    """Function to plot the methylation probability of positive strands. It takes as input:
        1)score: a vector contaiing for each position the porbability to methylate given by a model
        2)seq: the concerned sequence
        3)pos_met: the position of all the detected methylated sites
    It takes as optional parameters:
        4) Dict_filtr: a dictionary contianing all the submotifs to filter out results (default to empty vector)
        5) lw: parameter to se the line width (deafult it 1)
        6) ww: the window size from a central adenosine to filter the motifs (for DRACH is 2, for RAC is 1). Deafult is zero
        7) sizes: a tuple containing the size of the window figure
        8) tresh: the value of the treshold above which sites are predicted methylated (0 by default)
    The function returns: 
        1) The upper window plot containing the landscape of probability.
    """
    #Finding all A positions
    poss=findall('A', seq)
    (ww!=0 ? ((score, filtr)= filter_score(score, seq, ww, Dict_filtr)) : nothing)
    figure(figsize=(sizes))
    if ww!=0 #Filtering  experimetnal methylated sites over the chosen set
        range_Dict=vcat([findall(a, seq) for a in keys(Dict_filtr)]...)
        pos_Dict=[range_Dict[i][ww+1] for i in 1:length(range_Dict)]
        pos_met=intersect(pos_met, pos_Dict)
    end
    #plotting experimental methylated sites 
    [plot([pos_met[i],pos_met[i]],[0,2], linestyle=":", color="lightsteelblue", linewidth=lw, alpha=1) for i in 1:length(pos_met), alpha=0.7]
    #Plotting probability profiles over sequences
    for i in 1:length(score) #looping over all sites
        if score[i]>=tresh #Colouring all sites with a probability higher than the treshold
            plot([i, i], [0.0, score[i]], color="goldenrod",  alpha=0.75, linewidth=lw ) 
        else #if  isempty(intersect(i, findall(x->x==1, filtr)))==false  #&& score[i]>0.01
            plot([i, i], [0.0, score[i]], color="lightgrey",  alpha=0.9, linewidth=lw )
        end
    end
    #Horizontal line of reference
    axhline(y=0,color="lightgrey",linewidth=lw)
    return gcf()
end


function plot_score_neg(score, seq, pos_met; Dict_filtr=[], ww=0, lw=1,sizes=(5,5), tresh=0)
    """Function to plot the methylation probability of negative strands. It takes as input:
    1)score: a vector contaiing for each position the porbability to methylate given by a model
    2)seq: the concerned sequence
    3)pos_met: the position of all the detected methylated sites
It takes as optional parameters:
    4) Dict_filtr: a dictionary contianing all the submotifs to filter out results (default to empty vector)
    5) lw: parameter to se the line width (deafult it 1)
    6) ww: the window size from a central adenosine to filter the motifs (for DRACH is 2, for RAC is 1). Deafult is zero
    7) sizes: a tuple containing the size of the window figure
    8) tresh: the value of the treshold above which sites are predicted methylated (0 by default)
The function returns: 
    1) The lower window plot containing the landscape of probability.
"""
    #Finding all A positions
    poss=findall('A', seq)
    (ww!=0 ? ((score, filtr)= filter_score(score, seq, ww, Dict_filtr)) : nothing)
    if ww!=0 #Filtering  experimetnal methylated sites over the chosen set
        range_Dict=vcat([findall(a, seq) for a in keys(Dict_filtr)]...)
        pos_Dict=[range_Dict[i][ww+1] for i in 1:length(range_Dict)]
        pos_met=intersect(pos_met, pos_Dict)
    end
    #plotting experimental methylated sites 
    [plot([pos_met[i],pos_met[i]],[0,-2], linestyle=":", color="lightsteelblue", linewidth=lw, alpha=1) for i in 1:length(pos_met), alpha=0.7]
    #Plotting probability profiles over sequences
    for i in 1:length(score) #looping over all sites
        if score[i]>=tresh #Colouring all sites with a probability higher than the treshold
            plot([i, i], -[0.0, score[i]], color="lightseagreen",  alpha=0.75, linewidth=lw ) 
        else #if isempty(intersect(i, findall(x->x==1, filtr)))==false #&& score[i]>0.01
            plot([i, i], -[0.0, score[i]], color="lightgrey",  alpha=0.9, linewidth=lw )
        end
    end
    xlim(0,length(seq)+1)
    yticks([-1, -0.5, 0, 0.5, 1], [1, 0.5, 0, 0.5, 1])
    ylim(-1.05,1.05)
    xticks(fontsize=14)
    yticks(fontsize=14)
    xlabel(L"sites", size=15)
    ylabel(L"probability", size=15)
    return gcf()
end