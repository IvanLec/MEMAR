#List of DRACH motifs 
const DRACH=["TAACA", "AAACA", "GAACA", "TGACA", "AGACA", "GGACA", "TAACC", "AAACC", "GAACC", "TGACC", "AGACC", "GGACC", "TAACT", "AAACT", "GAACT", "TGACT", "AGACT", "GGACT"];

#Size of the window before and after the A 
const window=4 

#Alphabet, choosen as DNA letters for simplicity (the choice does noy change the mdelization, it remains equivalent between RNA/DNA sequences)
const dna_alphabet=["A", "C", "G", "T"]; 

#Construction of all di-nucleotides
const mots_list2 = vcat([join.([collect(Iterators.product([dna_alphabet for _ in 1:k]...))...]) for k in 2:2]...);
#Construction of all tri-nucleotides
const mots_list3 = vcat([join.([collect(Iterators.product([dna_alphabet for _ in 1:k]...))...]) for k in 3:3]...);
#Construction of all 4-nucleotides
const mots_list4 = vcat([join.([collect(Iterators.product([dna_alphabet for _ in 1:k]...))...]) for k in 4:4]...);
#Construction of all 8-nucleotides
const mots_list8 = vcat([join.([collect(Iterators.product([dna_alphabet for _ in 1:k]...))...]) for k in 8:8]...);

#Strcuture to save start and stop codons for coding transcripts
mutable struct CDS 
    start
    stop
end

#Structure to store the forces and the log partition function of inferred models
mutable struct MaxEnt3  
    h #forces for 1-mers
    J #forces for 2-mers
    K #force for 3-mers
    logZ
end 

#Strcuture needed to store all the information for genes from GLORI-seq experiments
mutable struct Gene_info_t 
    sequence #Mature transcript for a given gene
    pos_met #Positions of methylated As in mature transcripts
    ratio #frequency of methylation (stoichiometry)
    expression #gene expression level in experiments
    cds #start-stop for coding RNA (or 0 element)
    perc_met #percentage between m6A and A 
    splicing #junction position for mature RNA
    strand #strand of the gene over the mapped genome
    pos_not #Positions of all A non methylated
    selected_not #randomly selected sites amongst non methyleted As
end

#Strcuture needed to store all the information for genes from miCLIP experiments
mutable struct Gene_info_CLIP #Structure 
    sequence #Mature transcript for a given gene
    pos_met #Positions of methylated As in mature transcripts
    cds #start-stop for coding RNA (or 0 element)
    perc_met  #percentage between m6A and A 
    splicing #junction position for mature RNA
    strand #strand of the gene over the mapped genome
    pos_not #Positions of all A non methylated
end
