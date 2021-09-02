#!/mnt/e/Julia-1.6.2
using Base: Float64, Float32
using BioStructures

function get_next_line_frag(chain, frag_num, frag_start::Int64, chain_res_end, frag_matrix, pdb_name, model_num, chain_id)
    truncate_domain = collectatoms(chain, res -> resnumber(res) >= frag_start, calphaselector)
    frag_checkpoints::Int64 = 0
    update_matrix = false
    first_res = true
    residue_num_last::Int64 = 0
    frag_end::Int64 = 0
    for residue in truncate_domain
        b_factor_p::Float32 = tempfactor(residue)
        residue_num_p = resnumber(residue)
        chain_res_end == residue_num_p && begin
            frag_end = residue_num_p
            frag_checkpoints > 0 && ( update_matrix = true )
            break
        end
        first_res == true || begin 
            if residue_num_last != residue_num_p - 1
                frag_end = residue_num_last
                (frag_checkpoints > 0) && ( update_matrix = true )
                println("Note: Gap $residue_num_last - $residue_num_p in $pdb_name-$model_num-$chain_id")
                break
            end
        end
        residue_num_last = residue_num_p
        if b_factor_p >= b_factor_min && b_factor_p <= b_factor_max
            ( frag_checkpoints == 0 ) && ( frag_start = residue_num_p )
            ( frag_checkpoints <= disorder_len_max ) && ( frag_checkpoints = frag_checkpoints + 1 )
        else 
            if frag_checkpoints > 0
                frag_checkpoints = frag_checkpoints - 1
                ( frag_checkpoints == 0 ) && begin
                    frag_end = residue_num_p
                    update_matrix = true
                    break
                end 
            end
        end
        first_res = false
    end
    ( update_matrix == true ) && begin
        frag_len = frag_end - frag_start + 1
        frag_len >= domain_min_len && ( push!(frag_matrix,"$frag_num,$frag_start,$frag_end") )
    end
    return frag_matrix, frag_end
end

function get_terminal_resnum(chain::Chain)
    C_end::Int64 = 1
    N_start::Int64 = 9999999999 
    for residue in chain
        residue_num_p::Int64 = resnumber(residue) 
        ( C_end <= residue_num_p ) && ( C_end = residue_num_p )
        ( N_start >= residue_num_p ) && ( N_start = residue_num_p )
    end
    return N_start::Int64,C_end::Int64
end

function segment_domain(chain, frag_matrix, pdb_name, model_num, chain_id,index_dir,pdb,output_dir)
    res_domain = length(frag_matrix)
    domain_ID = 0
    domain_info = open("$index_dir/$pdb_name-$model_num-$chain_id-domain.info","w")
    while res_domain > 0
        domain = ""
        domain_ID+=1
        got_domain = []
        frag_id = 0
        for frag in frag_matrix
            frag_id+=1
            frag_start = parse(Int64, split(frag, ",")[2])
            frag_end = parse(Int64, split(frag, ",")[3])
            frag_str = collectresidues(chain, res -> frag_start <= resnumber(res) <= frag_end )
            ( frag_id == 1 ) && begin 
                domain = "$frag_start-$frag_end"
                push!(got_domain, frag)
            end
            ( frag_id != 1 ) && begin
                for frag_2 in split(domain,",")
                    frag_2_str = collectresidues(chain, res -> parse(Int64,split(frag_2,"-")[1]) <= resnumber(res) <= parse(Int64,split(frag_2,"-")[2]))
                    Adjacent_residues = 0
                    for residue in frag_2_str
                        frag_dist::Float64 = BioStructures.distance(frag_str, residue, allselector)
                        ( frag_dist <= min_dist_interdomain ) && ( Adjacent_residues+=1 )
                    end
                    ( min_Adjacent_residues <= Adjacent_residues ) && begin
                        domain = "$(domain),$frag_start-$frag_end"
                        push!(got_domain, frag)
                        break
                    end
                end
            end
        end
        for frag in got_domain
            filter!(x->x!=frag, frag_matrix)
        end
        res_domain = length(frag_matrix)
        write(domain_info, "$domain_ID:$domain\n")
    end
    close(domain_info)
    for domain_ddd in readlines("$index_dir/$pdb_name-$model_num-$chain_id-domain.info")
        domain_id = split(domain_ddd,":")[1]
        domain_vec = split(domain_ddd,":")[2]
        com_id = 0
        domain_selection = []
        for domain_range in split(domain_vec,",")
            frag_N = parse(Int64,split(domain_range,"-")[1])
            frag_C = parse(Int64,split(domain_range,"-")[2])
            for resi in frag_N:frag_C
                push!(domain_selection, resi)
            end
        end
        domain_str = collectresidues(chain, res -> resnumber(res) in domain_selection , allselector)
        writepdb("$output_dir/$pdb_name-$model_num-$chain_id-$domain_id.pdb", domain_str)
    end    
end

function main(pdb::String, index_dir::String, output_dir::String)
    Structure = read(pdb, PDB)
    pdb_name = chop(basename(pdb), tail = 4)
    for model in Structure
        model_num = modelnumber(model)
        for chain in model
            chain_id = chainid(chain)
            chain_res_start::Int64, chain_res_end::Int64= get_terminal_resnum(chain)
            frag_start::Int64 = chain_res_start
            frag_end::Int64 = chain_res_start
            frag_num = 0
            frag_matrix = []
            while frag_end < chain_res_end
                frag_num+=1
                frag_start = frag_end
                #将 frag_num, frag_start, frag_end 存入 frag_matrix
                frag_matrix, frag_end= get_next_line_frag(chain, frag_num, frag_start, chain_res_end, frag_matrix, pdb_name, model_num, chain_id)
            end
            segment_domain(chain, frag_matrix, pdb_name, model_num, chain_id,index_dir,pdb,output_dir)
        end
    end
end

println("Name of Input Dir for protein pdb files:")
const input_dir=readline(stdin)
base_input = basename(input_dir)
index_dir="$(base_input)-INDEX"
mkdir(index_dir)
println("Output Dir for domain Structure:")
const output_dir=readline(stdin)
global b_factor_max = 100 # max b_factor(pLDDT) to keep a residue.
global b_factor_min = 50 # min b_factor(pLDDT) to keep a residue.
global disorder_len_max = 5 # the max disorder tail in a ordered domain.
global domain_min_len = 10   # this option is means of the min length of fargs. The domain is consist of frags incessant in sequence.
global min_dist_interdomain = 3.50  # 3.00 is a very strict cutoff.
global min_Adjacent_residues = 5 # How many residues are close together between frags will be identified as a domain.

Threads.@threads for pdb_file in readdir(input_dir)
    pdb = "$(input_dir)/$(pdb_file)"
    endswith(pdb, ".pdb") && begin 
        try
            main(pdb, index_dir, output_dir)
        catch e1
            println("$e1 Error: $pdb")
        end
    end
end
