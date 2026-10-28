const CHR_TO_LENGTH=Dict("chr1"=>   249250621, 	"chr2"=> 	243199373 	,"chr3"=> 	198022430 ,"chr4"=> 	191154276 ,"chr5"=> 	180915260 ,"chr6"=> 	171115067 ,"chr7"=> 	159138663 ,"chrX"=> 	155270560 ,"chr8"=> 	146364022 ,"chr9"=>	141213431 ,	"chr10"=> 	135534747 ,"chr11"=> 	135006516 ,	"chr12"=> 	133851895 ,"chr13"=> 	115169878,"chr14"=> 	107349540 ,"chr15"=> 	102531392 ,"chr16"=> 	90354753 ,"chr17"=> 	81195210,"chr18"=> 	78077248 ,"chr20"=> 	63025520 ,	"chrY"=> 	59373566,"chr19"=> 	59128983 ,	"chr22"=> 	51304566 ,	"chr21"=>	48129895);
const start_codon="ATG"
const stop_codon=["TAG", "TGA", "TAA"];

mutable struct CDS #Structure to save start and stop codon
    start
    stop
end
#Defining the position fo some labels
const column_strand=5
const column_pos=3
const column_gene=16
const column_chr=2
const column_ref=11;
const genes_position=

mutable struct Gene_info_CLIP #Structure 
    sequence
    pos_met
    cds
    perc_met
    splicing
    strand
    pos_not
end

const Length_winodw=4

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
    """Function to find CDS annotation
    The function takes as entry:
        -Transcripts: data of human Transcripts
    The function returns:
        -CDS_vec: a vector containing either a vector of two containg start ands stop codon position, either a 0-elemnent vector
         which indicates no CDS annotation can be found, so non-coding RNA flag
    """
    #Creating empty vector
    gene_vec=[]  
    for i in 1:length(Transcripts)
        for j in 1:length(Transcripts[i]) 
            if count("   gene     ", Transcripts[i][j])==1 #finding the label gene (Spaceces to distinguish possible amminoacid combinations)
                push!(gene_vec, split(Transcripts[i][j+1], '"')[2])
            end
        end
    end
    return gene_vec
end

function longest_transcripts(Transcripts, gene_vec)
    """Function to select among all transcripts (which for each gene multiple versions are presents) the longest ones.
    The function takes as entry:
        -Transcripts: data of human Transcripts
    The function returns:
        -D: Dictionary containing the label of the genes, and the longest versions
    """
    D=Dict()
    ind=Dict()
    for i in 1:length(Transcripts)
        if haskey(D, gene_vec[i])==true #Checking if this gene is already present in the dictionary
            if length(D[gene_vec[i]])<length(Transcripts[i][2]) #if it is, and the old one is shorter than the new we substitute with the new version
                D[gene_vec[i]]=Transcripts[i][2]
                ind[gene_vec[i]]=i #substituting the new index
            end
        else
            D[gene_vec[i]]=Transcripts[i][2] #if the gene was not present, we save it for the first time
            ind[gene_vec[i]]=i #saving for the first time the index
        end
    end
    return D, ind
end

function Annotations_genome(file)
    Data=split(file, "\n")
    T=[]
    for i in 1:length(Data)-1
        if Data[i][1]!='#'
            push!(T, Data[i])
        end
    end
    return T 
end

function exons(human_transcripts,indexes_long,Annotations_dna)
    Code_genes=Dict()
    [(Code_genes[("ID=rna-"*split( human_transcripts[a][1], " ")[1])]=[]) for a in indexes_long];
    for i in 1:length(Annotations_dna)
        K=split(split(Annotations_dna[i], "\t")[9],";")[1]
        if haskey(Code_genes, K)==true 
            j=i+1 
            while split(Annotations_dna[j], "\t")[3]=="exon"
                push!(Code_genes[K],[parse(Int64,split(Annotations_dna[j], "\t")[4]),parse(Int64,split(Annotations_dna[j], "\t")[5])])
                j+=1
            end
        end
    end
    GENE_to_REF=Dict()
    [(GENE_to_REF[gene_vec[a]]=("ID=rna-"*split( human_transcripts[a][1], " ")[1])) for a in indexes_long];
    return Code_genes, GENE_to_REF
end


function mapping(df, Code_genes,GENE_to_REF)
    #Reading experimental results as DataFrame
    #Creating a dictiony with the chromosome length (relative to GRCh38)
    positions_A=Dict()
    positions_cov=Dict()
    GENE_to_STRAND=Dict()
    GENE_to_CHR=Dict()
    for i in 1:size(df)[1]
        if haskey(GENE_to_REF,df[i, column_gene])==true #is there the gene name from data in the annotation file?
            #I check i the gene maps into the positive or negative strand of the human genome
            if df[i, column_strand]=="+" 
                LL=0
                for j in 1:length(Code_genes[GENE_to_REF[df[i, column_gene]]]) #for each exon which is present  
                    a=Code_genes[GENE_to_REF[df[i, column_gene]]][j][1] #I save the starting position of the exon
                    if a<= df[i, column_pos] <=Code_genes[GENE_to_REF[df[i, column_gene]]][j][2] #I check if the mapping position is inside the given range
                        if haskey(positions_A,df[i, column_gene])==true #If I have already information about this gene I save only the position, mapping it into the transcipt
                            #The real position in the transcript will be the mapping on the geno minus the staring position of the exon +1, plus all the already visited exons
                            push!(positions_A[df[i, column_gene]], df[i, column_pos]-a+1+LL)
                        else
                            #If it the first time I save strand infromation and the chromosome
                            GENE_to_STRAND[df[i, column_gene]]=df[i, column_strand]
                            positions_A[df[i, column_gene]]=[df[i, column_pos]-a+1+LL]
                            GENE_to_CHR[df[i, column_gene]]=df[i, column_chr]
                        end
                        #Inside the found genes, I count the actual number of mapping sites retrived inside the correct ranges
                    end
                    #if the position is not find I simply in the exon I add the exon length
                    LL+=Code_genes[GENE_to_REF[df[i, column_gene]]][j][2]-a+1
                end
            else
                #In the negative case, due to inverse characteristics between transcript and genomes the mapping is sligthly differret
                LL=0
                R=length(Code_genes[GENE_to_REF[df[i, column_gene]]])
                for j in 1:R #for each exon which 
                    a=Code_genes[GENE_to_REF[df[i, column_gene]]][j][1]
                    if a<= df[i, column_pos] <=Code_genes[GENE_to_REF[df[i, column_gene]]][j][2]
                        if haskey(positions_A,df[i, column_gene])==true
                            push!(positions_A[df[i, column_gene]], CHR_TO_LENGTH[df[i, column_chr]].-df[i, column_pos] .- (CHR_TO_LENGTH[df[i, column_chr]]-Code_genes[GENE_to_REF[df[i, column_gene]]][j][2]).+1+LL )
                        else
                            GENE_to_STRAND[df[i, column_gene]]=df[i, column_strand]
                            positions_A[df[i, column_gene]]=[CHR_TO_LENGTH[df[i, column_chr]].-df[i, column_pos] .- (CHR_TO_LENGTH[df[i, column_chr]]-Code_genes[GENE_to_REF[df[i, column_gene]]][j][2]).+1+LL ]
                            GENE_to_CHR[df[i, column_gene]]=df[i, column_chr]
                        end
                    end
                    LL+=Code_genes[GENE_to_REF[df[i, column_gene]]][j][2]-a+1
                end
            end
        end
    end
    return positions_A,positions_cov,GENE_to_STRAND,GENE_to_CHR
end

function Annotations_CDS(file)
    """Function to read NCBI annotations file of human transcripts
    The function takes as input:
        -file: path of the file to read
    The function returns:
        -CDS_intervals: vector of the same number of human transcipts with either a CDS structure containg start and stop codon potion or empty vector, indicating non-coding RNA
    """
    #Creating empty data
    Data=[]
    #Opening file in reading modes
    #=open(file, "r") do file
        for line in eachline(file)
            push!(Data,line) #reading all lines
        end
    end=#
    Data=split(file, "\n")
    #Finding all sequences presented, which are separated by "//" in the text file
    indexes=vcat(0, findall(x->x=="//", Data)...);
    Transcipts=[Data[indexes[i]+1: indexes[i+1]] for i in 1:length(indexes)-1]
    CDS_vec=finding_CDS(Transcipts)
    CDS_intervals=parse_CDS(CDS_vec)
    return CDS_intervals 
end

function finding_CDS(Transcripts)
    """Function to find CDS annotation
    The function takes as entry:
        -Transcripts: data of human Transcripts
    The function returns:
        -CDS_vec: a vector containing either a vector of two containg start ands stop codon position, either a 0-elemnent vector
         which indicates no CDS annotation can be found, so non-coding RNA flag
    """
    #Creating empty vector
    CDS_vec=[]  
    for i in 1:length(Transcripts)
        flag=0 #initialization of a flag to point if a transcrit is a non-coding RNA
        for j in 1:length(Transcripts[i]) 
            if count("CDS   ", Transcripts[i][j])!=0 #finding the label CDS (Spaceces to distinguish possible amminoacid combinations)
                push!(CDS_vec, Transcripts[i][j])
                flag+=1
            end
        end
        if flag==0
            push!(CDS_vec, flag) #Int euqal to 0 to signal non-coding RNA
        end
    end
    return CDS_vec 
end

function parse_CDS(CDS_vec)
    """Function to read a .txt line and to convert strings representative of start and stop codon position into a vector of numbers
    The function takes as entry:
        -CDS_vec: a vector containing either a vector of two containg start ands stop codon position, either a 0-elemnent vector
         which indicates no CDS annotation can be found, so non-coding RNA flag
        -CDS_intervals: vector of the same number of human transcipts with either a CDS structure containg start and stop codon potion or empty vector, indicating non-coding RNA
    """
    CDS_intervals=[] #Empty vector
    for i in 1:length(CDS_vec)
        if CDS_vec[i]==0 #when equal to the flag
            push!(CDS_intervals, CDS_vec[i])
        else
            matches = collect(eachmatch(r"\d+", CDS_vec[i]))
            numbers = [parse(Int, m.match) for m in matches]
            push!(CDS_intervals, CDS(numbers[1], numbers[2]))
        end
    end
    return CDS_intervals 
end

function mapping_details(df,positions_A,Code_genes,GENE_to_REF, Annotations_rna,indexes_long,gene_vec)
    GENE_to_SPLICING=Dict()
    for K in keys(positions_A)
        LL=[]
        for j in 1:length(Code_genes[GENE_to_REF[K]]) 
            a=Code_genes[GENE_to_REF[K]][j][1]
            if j!=length(Code_genes[GENE_to_REF[K]]) 
                push!(LL, (isempty(LL)==true ? 0 : LL[end]) +Code_genes[GENE_to_REF[K]][j][2]-a+1) 
            end
        end
        GENE_to_SPLICING[K]=LL 
    end    
    CDS_interval=Annotations_CDS(Annotations_rna);
    Dict_CDS=Dict()
    [Dict_CDS[gene_vec[indexes_long[i]]]=CDS_interval[indexes_long[i]] for i in 1:length(indexes_long)]
    return GENE_to_SPLICING,Dict_CDS
end

function mapping_to_structure(positions_A,positions_cov,Dict_CDS,GENE_to_SPLICING,GENE_to_STRAND,long_transcripts)
    GENE_to_INFO=Dict()
    for a in vcat(keys(positions_A)...)
        ind=sortperm(positions_A[a])
        pos_not=setdiff!(findall('A',long_transcripts[a]), positions_A[a])
        GENE_to_INFO[a]=Gene_info_CLIP(long_transcripts[a], 
        sort(positions_A[a]),
        Dict_CDS[a],
        length(positions_A[a])/count("A",long_transcripts[a])*100, 
        GENE_to_SPLICING[a], 
        GENE_to_STRAND[a],
        pos_not)
    end
    for KK in keys(GENE_to_INFO)
        if count("A", GENE_to_INFO[KK].sequence[GENE_to_INFO[KK].pos_met])!= length(GENE_to_INFO[KK].pos_met) 
            delete!(GENE_to_INFO, KK)
        end
    end
    return GENE_to_INFO
end

function mapping_wrapped(human_transcripts,indexes_long,Annotations_dna,Annotations_rna,gene_vec,Dataset_m6A, long_transcripts)
    Code_genes, GENE_to_REF=exons(human_transcripts,indexes_long,Annotations_dna);
    positions_A,positions_cov,GENE_to_STRAND,GENE_to_CHR=mapping(Dataset_m6A, Code_genes,GENE_to_REF)
    GENE_to_SPLICING,Dict_CDS=mapping_details(Dataset_m6A,positions_A,Code_genes,GENE_to_REF, Annotations_rna,indexes_long,gene_vec)
    GENE_to_INFO=mapping_to_structure(positions_A,positions_cov,Dict_CDS,GENE_to_SPLICING,GENE_to_STRAND,long_transcripts)
    return GENE_to_INFO
end