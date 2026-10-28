mutable struct Gene_info_t #Structure to save mapped m6As in transcripts
    sequence
    pos_met #mapped methylation position
    ratio #stoichiometry
    expression #gene expression
    cds #START and STOP codon indexes
    perc_met #percentage of methylated sites in the given transcript
    splicing #Splicing sites
    strand #strand polarity
    pos_not #position of all non-methylated adenosines
    selected_not #random selection of non-mehtylated adenosines
end
const start_codon="ATG"
const stop_codon=["TAG", "TGA", "TAA"];

mutable struct CDS #Structure to save START and STOP codon
    start
    stop
end

#Chromosome length of Reference Genome
const CHR_TO_LENGTH=Dict("chr1"=>248956422,"chr2"=>242193529,"chr3"=>198295559,"chr4"=>190214555,"chr5"=>181538259,"chr6"=>170805979,"chr7"=>159345973,"chrX"=>156040895,"chr8"=>145138636,"chr9"=>138394717,"chr10"=>133797422,"chr11"=>135086622,"chr12"=>133275309,"chr13"=>114364328,"chr14"=>107043718,"chr15"=>101991189,"chr16"=>90338345,"chr17"=>83257441,"chr18"=>80373285,"chr20"=>64444167,"chrY"=>57227415,"chr19"=>58617616,"chr22"=>50818468,"chr21"=>46709983);

const Length_winodw=150 #arbitrary window to select motifs

const genes_position=4; #The dataframes in position 4 from experimental data contain the labels of the transcripts where methylations have been found 


function Annotations_gene(file)
    """Function to read NCBI annotations file of human transcripts
    The function takes as input:
        -file: path of the file to read
    The function returns:
        -CDS_intervals: vector of the same number of human transcipts with either a CDS structure containg start and stop codon potion or empty vector, indicating non-coding RNA
    """
    #Creating empty data
    Data=split(file, "\n")
    #Finding all sequences presented, which are separated by "//" in the text file
    indexes=vcat(0, findall(x->x=="//", Data)...);
    Transcipts=[Data[indexes[i]+1: indexes[i+1]] for i in 1:length(indexes)-1]
    gene_vec=finding_gene(Transcipts)
    return gene_vec
end

function finding_gene(Transcripts)
    """Function to find transcripts
    The function takes as entry:
        -Transcripts: data of human transcripts
    The function returns:
        -gene_vec: all transcript labels in the annotation file of human data
    """
    #Creating empty vector
    gene_vec=[]  
    for i in 1:length(Transcripts)
        for j in 1:length(Transcripts[i]) 
            #finding the label gene (Spaces to distinguish correctly "gene" row specification)
            if count("   gene     ", Transcripts[i][j])==1 
                push!(gene_vec, split(Transcripts[i][j+1], '"')[2])
            end
        end
    end
    #output
    return gene_vec
end

function longest_transcripts(Transcripts, gene_vec)
    """Function to select amongst all transcripts the longest ones (for each gene, indeed, multiple versions are presents).
    The function takes as entry:
        -Transcripts: data of human transcripts
        -gene_vec: all transcript labels in the annotation file of human data
    The function returns:
        -D: Dictionary having as keys the label of the genes, and as values the longest transcripts
        -ind: the index corresponding to the position of the found longest transcripts in the human transcripts list
    """
    D=Dict()
    ind=Dict()
    for i in 1:length(Transcripts)
         #Checking if this gene is already present in the dictionary
        if haskey(D, gene_vec[i])==true
             #if it is, and the old one is shorter than the new we substitute with the new version
            if length(D[gene_vec[i]])<length(Transcripts[i][2])
                D[gene_vec[i]]=Transcripts[i][2]
                ind[gene_vec[i]]=i #substituting the new index
            end
        else
            #if the gene was not present, we save it for the first time
            D[gene_vec[i]]=Transcripts[i][2] 
            ind[gene_vec[i]]=i #saving for the first time the index
        end
    end
    return D, ind
end

function Annotations_CDS(file)
    """Function to read NCBI annotations file of human transcripts
    The function takes as input:
        -file: path of the file to annotation to read
    The function returns:
        -CDS_intervals: a list containing either a vector of length two containg START and STOP codon position, either a 0-element vector
         which indicates that no CDS annotation were found (therefore: non-coding RNA) 
    """
    #Creating empty data
    Data=[]
    Data=split(file, "\n")
    #Finding all sequences presented, which are separated by "//" in the text file
    indexes=vcat(0, findall(x->x=="//", Data)...);
    Transcipts=[Data[indexes[i]+1: indexes[i+1]] for i in 1:length(indexes)-1]
    CDS_vec=finding_CDS(Transcipts)
    CDS_intervals=parse_CDS(CDS_vec)
    return CDS_intervals 
end

function Annotations_genome(file)
    """Function to read the annotatio file of the human geneomes.
    The function takes as input: 
        -file: path of the file to annotation genome to read
    The function returns:
        -T: all annoted genomes
    """
    Data=split(file, "\n")
    T=[]
    for i in 1:length(Data)-1
        if Data[i][1]!='#'
            push!(T, Data[i])
        end
    end
    return T 
end

function finding_CDS(Transcripts)
    """Function to find CDS annotation.
    The function takes as entry:
        -Transcripts: data of human transcripts
    The function returns:
        -CDS_vec: a list containing either a vector of length two containg START and STOP codon position, either a 0-element vector
         which indicates that no CDS annotation were found (therefore: non-coding RNA) 
    """
    #Creating empty vector
    CDS_vec=[]  
    for i in 1:length(Transcripts)
        flag=0 #initialization of a flag to point if a transcrit is a non-coding RNA
        for j in 1:length(Transcripts[i]) 
            if count("CDS   ", Transcripts[i][j])!=0 #finding the label CDS (Spaces to distinguish "CDS" label from possible amminoacid motifs)
                push!(CDS_vec, Transcripts[i][j])
                flag+=1
            end
        end
        if flag==0
            push!(CDS_vec, flag) #return the integer number 0 to signal non-coding RNAs
        end
    end
    #output
    return CDS_vec 
end


function parse_CDS(CDS_vec)
    """Function to read a file ".txt" line and to convert strings representative of START and STOP codon position into a vector of numbers
    The function takes as entry:
        -CDS_vec: a vector containing either a vector of two containg start ands stop codon position, either a 0-elemnent vector
         which indicates no CDS annotation can be found, (non-coding RNA)
        -CDS_intervals: vector of the same number of human transcripts with either a CDS structure containg start and stop codon potion or empty vector, indicating non-coding RNAs
    """
    CDS_intervals=[] #Empty vector
    for i in 1:length(CDS_vec)
        if CDS_vec[i]==0 #when equals to the flag
            push!(CDS_intervals, CDS_vec[i])
        else
            matches = collect(eachmatch(r"\d+", CDS_vec[i]))
            numbers = [parse(Int, m.match) for m in matches]
            push!(CDS_intervals, CDS(numbers[1], numbers[2]))
        end
    end
    return CDS_intervals 
end

function exons(human_transcripts,indexes_long,Annotations_dna)
    """Function to find the position of the exons:
    The function takes as entries:
        -human transcripts: data of human transcripts
        -indexes_long: the position of the longest transcripts in the data
        -Annotations_dna: the annotation infromation of the genome transcribing for the given trascript.
        The function returns:
            -Code_genes: the list of exons corresponding to the given genome's rasncript
            -GENE_to_REF: the code to identify the genome from the given transcript name
        """
    Code_genes=Dict() #initialization of an empty vocabulary
    #finding the code of each transcript to identify the corresponding genome in the DNA annotation file
    [(Code_genes[("ID=rna-"*split( human_transcripts[a][1], " ")[1])]=[]) for a in indexes_long];
    for i in 1:length(Annotations_dna) #looping over all annotated genomes
        K=split(split(Annotations_dna[i], "\t")[9],";")[1] #finding the identification code of the given genome
        if haskey(Code_genes, K)==true 
            j=i+1 
            while split(Annotations_dna[j], "\t")[3]=="exon"
                #Saving the list of exons
                push!(Code_genes[K],[parse(Int64,split(Annotations_dna[j], "\t")[4]),parse(Int64,split(Annotations_dna[j], "\t")[5])])
                j+=1
            end
        end
    end
    #Saving the genome code from the transcript name
    GENE_to_REF=Dict()
    [(GENE_to_REF[gene_vec[a]]=("ID=rna-"*split( human_transcripts[a][1], " ")[1])) for a in indexes_long];
    #outputs
    return Code_genes, GENE_to_REF
end


function mapping(df, Code_genes,GENE_to_REF, p_cov)
    """Function to map experimental datasets into longest trasncripts.
    The function takes as input:
        -df: the dataframe containing the exeperimental results
        -Code_genes: the list of exons corresponding to the given genome's rasncript
        -GENE_to_REF: the code to identify the genome from the given transcript name
        -p_cov: the position stoichiometry of methylated adenosines
        The functions returns:
        -positions_A: the mapped position of methylated adenosine in the transcript
        -positions_cov: the stoichiometry of each methylated adenosine
        -GENE_to_STRAND: a dictionary containing as keys the transcript names and a values the corresponding DNA strand transribing for it
        -GENE_to_CHR: a dictonary containg as keys transcript names and a values the chromosone containing the genome transcribing for it
    """
    #initialization
    positions_A=Dict()
    positions_cov=Dict()
    GENE_to_STRAND=Dict()
    GENE_to_CHR=Dict()
    for i in 1:size(df)[1] #For all methylatd adenosines
        if haskey(GENE_to_REF,df[i,4])==true #cheking if the gene name is present in the gene list from annotation file?
            if df[i,3]=="+" #checking if the genome coding for the transcript has positive polarity
                LL=0 #initialization of the length to check for exon bounds 
                for j in 1:length(Code_genes[GENE_to_REF[df[i,4]]]) #for each exon coding for the longest trasncript
                    a=Code_genes[GENE_to_REF[df[i,4]]][j][1] #saving, for the given exon, the lower bound
                    if a<= df[i,2] <=Code_genes[GENE_to_REF[df[i,4]]][j][2] #cheking the mapped position is inside the exon bounds
                        #Cheking if the given transcript has already had a methylated mapped adenosine
                        if haskey(positions_A,df[i,4])==true 
                            #In that case, I simply add the position of the methylated adenosine mapped over the longest isofrom postion in the list of methylated adenosine for that transcript
                            push!(positions_A[df[i,4]], df[i,2]-a+1+LL)
                            #Saving the stoichiometr (eventually converting the string number into a float64)
                            push!(positions_cov[df[i,4]], (typeof(df[i,p_cov])==Float64 ? df[i,p_cov] : parse(Float64,df[i,p_cov]))  )
                            
                        else
                            #If it is the first time for a mapped position in a given transcript I save the strand of DNA transcribing for it, and the chromosome
                            GENE_to_CHR[df[i,4]]=df[i,1]
                            GENE_to_STRAND[df[i,4]]=df[i,3]
                            positions_cov[df[i,4]]=[(typeof(df[i,p_cov])==Float64 ? df[i,p_cov] : parse(Float64,df[i,p_cov]))  ]
                            positions_A[df[i,4]]=[df[i,2]-a+1+LL]
                            
                        end
                    end
                    #adding the length till the exon considered 
                    LL+=Code_genes[GENE_to_REF[df[i,4]]][j][2]-a+1
                end
            else 
                #in case of the negative DNA polarity, we repeat the same procedure, although, inverting the sequence to map positions
                LL=0
                R=length(Code_genes[GENE_to_REF[df[i,4]]])
                for j in 1:R #for each exon 
                    #checking exon bounds
                    a=Code_genes[GENE_to_REF[df[i,4]]][j][1]
                    if a<= df[i,2] <=Code_genes[GENE_to_REF[df[i,4]]][j][2]
                        if haskey(positions_A,df[i,4])==true
                            #saving mapped position and stoichiometry
                            push!(positions_A[df[i,4]], CHR_TO_LENGTH[df[i,1]].-df[i,2] .- (CHR_TO_LENGTH[df[i,1]]-Code_genes[GENE_to_REF[df[i,4]]][j][2]).+1+LL )
                            push!(positions_cov[df[i,4]], (typeof(df[i,p_cov])==Float64 ? df[i,p_cov] : parse(Float64,df[i,p_cov]))  )
                        else
                            #saving strand polarity and chromosome
                            GENE_to_CHR[df[i,4]]=df[i,1]
                            GENE_to_STRAND[df[i,4]]=df[i,3]
                            positions_cov[df[i,4]]=[(typeof(df[i,p_cov])==Float64 ? df[i,p_cov] : parse(Float64,df[i,p_cov]))  ]
                            positions_A[df[i,4]]=[CHR_TO_LENGTH[df[i,1]].-df[i,2] .- (CHR_TO_LENGTH[df[i,1]]-Code_genes[GENE_to_REF[df[i,4]]][j][2]).+1+LL ]
                            
                        end
                    end
                    #adding the length till the exon considered 
                    LL+=Code_genes[GENE_to_REF[df[i,4]]][j][2]-a+1
                end
            end
        end
    end
    #outputs
    return positions_A,positions_cov,GENE_to_STRAND,GENE_to_CHR
end

function mapping_details(df,positions_A,Code_genes,GENE_to_REF, Annotations_rna,indexes_long,gene_vec, p_exp)
    """Function that give additionally detailes over the m6A mapped transcripts.
        The functions takes as entry: 
        -df: the dataframe containing the exeperimental results
        -positions_A: the mapped position of methylated adenosine in the transcript
        -Code_genes: the list of exons corresponding to the given genome's rasncript
        -GENE_to_REF: the code to identify the genome from the given transcript name
        -Annotations_rna: annotation file details of human transcripts 
        -indexes_long: the position of the longest transcripts in the data
        -gene_vec: all transcript labels in the annotation file of human data
        -p_exp: position of the gene expression in the experimental data
    The function returns: 
        -GENE_to_SPLICING: position of the splicing sites in the longest trasncript (for each gene)
        -Dict_gene_expression: experiemntal expression of each gene
        -Dict_CDS: codon START and STOP position for each transcript

    """
    #initialization
    GENE_to_SPLICING=Dict()
    #looping over all mapped transcripts
    for K in keys(positions_A)
        LL=[] #saving all exon lengths for the mature transcript
        for j in 1:length(Code_genes[GENE_to_REF[K]]) 
            a=Code_genes[GENE_to_REF[K]][j][1]
            if j!=length(Code_genes[GENE_to_REF[K]]) 
                push!(LL, (isempty(LL)==true ? 0 : LL[end]) +Code_genes[GENE_to_REF[K]][j][2]-a+1) 
            end
        end
        GENE_to_SPLICING[K]=LL 
    end    
    genes=unique(df[:,genes_position]); #The function unique allows to give the list of all genes present from experimental results withouth any redundence. 
    indexes_unique_genes = [findfirst(==(x), df[:, 4]) for x in genes]; 
    #For ach unique gene, the gene expression is saved
    Dict_gene_expression=Dict()
    [Dict_gene_expression[df[indexes_unique_genes[i], 4]]=df[indexes_unique_genes[i], p_exp] for i in 1:length(indexes_unique_genes)];
    CDS_interval=Annotations_CDS(Annotations_rna);
    Dict_CDS=Dict()
    #For ach unique gene, the bounds of the coding region is saved
    [Dict_CDS[gene_vec[indexes_long[i]]]=CDS_interval[indexes_long[i]] for i in 1:length(indexes_long)]
    return GENE_to_SPLICING,Dict_gene_expression,Dict_CDS
end

function mapping_to_structure(positions_A,positions_cov,Dict_gene_expression,Dict_CDS,GENE_to_SPLICING,GENE_to_STRAND,long_transcripts, verbose)
    """Function that saves all experimental m6A sites into the mapped ones over the longest transcripts into a structure.
    The function takes as entry: 
        -positions_A: the mapped position of methylated adenosine in the transcript
        -position_cov: the position stoichiometry of methylated adenosines
        -Dict_gene_expression: experiemntal expression of each gene
        -Dict_CDS: codon START and STOP position for each transcript
        -GENE_to_SPLICING: position of the splicing sites in the longest trasncript (for each gene)
        -GENE_to_STRAND: a dictionary containing as keys the transcript names and a values the corresponding DNA strand transribing for it
        -long_transcripts: gene longest isoforms
        -veborse: boolean variable

    """
    GENE_to_INFO=Dict()
    for a in vcat(keys(positions_A)...)
        ind=sortperm(positions_A[a])
        pos_not=setdiff!(findall('A',long_transcripts[a]), positions_A[a])
        T=shuffle!(pos_not[findall(x-> Length_winodw<x<=length(long_transcripts[a])-Length_winodw, pos_not)])
        #For each transcript saving a strcuture containg all traucture information
        GENE_to_INFO[a]=Gene_info_t(long_transcripts[a], 
        sort(positions_A[a]),
        positions_cov[a][ind],
        Dict_gene_expression[a],
        Dict_CDS[a],
        length(positions_A[a])/count("A",long_transcripts[a])*100, 
        GENE_to_SPLICING[a], 
        GENE_to_STRAND[a],
        pos_not,
        sort(T[1:minimum([9*length(positions_A[a]),length(T)])]))
    end
    #checking that all methylated adenosines are mapped into adenosines positions
    if verbose==true
        for KK in keys(GENE_to_INFO)
            if count("A", GENE_to_INFO[KK].sequence[GENE_to_INFO[KK].pos_met])!= length(GENE_to_INFO[KK].pos_met) || length(GENE_to_INFO[KK].selected_not)<9*length(GENE_to_INFO[KK].pos_met)
                delete!(GENE_to_INFO, KK)
            end
        end
    else
        for KK in keys(GENE_to_INFO)
            if count("A", GENE_to_INFO[KK].sequence[GENE_to_INFO[KK].pos_met])!= length(GENE_to_INFO[KK].pos_met) 
                delete!(GENE_to_INFO, KK)
            end
        end
    end
    return GENE_to_INFO
end

function mapping_wrapped(human_transcripts,indexes_long,Annotations_dna,Annotations_rna,gene_vec,Dataset_m6A,long_transcripts, p_cov, p_exp; verb=false)
    Code_genes, GENE_to_REF=exons(human_transcripts,indexes_long,Annotations_dna);
    positions_A,positions_cov,GENE_to_STRAND,GENE_to_CHR=mapping(Dataset_m6A, Code_genes,GENE_to_REF, p_cov)
    GENE_to_SPLICING,Dict_gene_expression,Dict_CDS=mapping_details(Dataset_m6A,positions_A,Code_genes,GENE_to_REF, Annotations_rna,indexes_long,gene_vec, p_exp)
    GENE_to_INFO=mapping_to_structure(positions_A,positions_cov,Dict_gene_expression,Dict_CDS,GENE_to_SPLICING,GENE_to_STRAND,long_transcripts, verb)
    return GENE_to_INFO
end

function pandas_to_julia(df_py)
    df_jl = DataFrame()

    for col in df_py.columns
        df_jl[!, Symbol(col)] = df_py[col].tolist()
    end

    return df_jl
end

function split_train_val_test(Dataset_1, Dataset_2; train=[], validate=[])
    FILE1=convert_to_panda(Dataset_1)
    FILE2=convert_to_panda(Dataset_2)
    Dict_sites=Dict()
    Dict_sites, header=myscript.readfile(FILE1, Dict_sites)
    Dict_sites, header=myscript.readfile(FILE2, Dict_sites)
    if isempty(train)==true 
        train, validate, test=myscript.split_TVT(Dict_sites)
    end 
    Dataset_train=myscript.write_set(train, Dict_sites, header);
    Dataset_validate=myscript.write_set(validate, Dict_sites, header);
    Validation=pandas_to_julia(Dataset_validate)
    Train=pandas_to_julia(Dataset_train);
    Test1=deepcopy(Dataset_1)
    Test2=deepcopy(Dataset_2)
    genes=vcat(train..., validate...)
    delete!(Test1, findall(x->isempty(findall(y->y==x, genes))==false, Dataset_1[:, 4]))
    delete!(Test2, findall(x->isempty(findall(y->y==x, genes))==false, Dataset_2[:, 4]))
    return Train, Validation, Test1, Test2
end

function convert_to_panda(Dataset_m6A)
    FILE1=Dataset_m6A;
    cols = names(FILE1)                     
    data = [FILE1[!, c] for c in cols]      
    FILE1_py = pd.DataFrame(pycall(pybuiltin("zip"), PyObject, cols, data))
    FILE1_py = pd.DataFrame(Dict(c => FILE1[!, c] for c in cols), columns=cols);
    return FILE1_py
end