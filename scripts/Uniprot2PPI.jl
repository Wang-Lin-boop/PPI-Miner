using BioStructures
using CSV, Tables, DataFrames
using Parsers

function check_protein_proximity(structure)
    receptor = collectatoms(structure, ch -> chainid(ch) == "A" )
    ligand = collectatoms(structure, ch -> chainid(ch) == "B" )
    adjacent_atoms = 0 
    proximity_score = 0 
    for receptor_atom in collectatoms(receptor)
        atom_dist::Float64 = distance(ligand, receptor_atom)
        ( atom_dist < min_dist ) && ( adjacent_atoms+=1 )
        for ligand_atom in collectatoms(ligand)
            atom_pair_dist::Float64 = distance(ligand_atom, receptor_atom)
            ( atom_pair_dist < min_dist ) && ( proximity_score += epsilon * ((sigma/atom_pair_dist)^12 - (sigma/atom_pair_dist)^6 ) )
        end
    end
    return adjacent_atoms, proximity_score
end

function write_struct(structure, output_dir, struct_name)
    try
        writepdb("$output_dir/$struct_name.pdb", structure)
    catch e1
        write(error,"PDB writing Error: $output_dir/$struct_name !\n")
        writemmcif("$output_dir/$struct_name.cif", structure)
    end
end

function change_chain_names(ppi_complex, source_chainid, target_chainid, receptor_selection, ligand_selection)
    pseudo_chainname_list = ["X", "Y", "Z", "x", "z", "y", "1", "2", "0", "9", "99"]
    chainnames_in_complex = chainids(collectmodels(ppi_complex)[1])
    candiate_chainnames = filter(x -> !(x in chainnames_in_complex), pseudo_chainname_list)
    for chain in ppi_complex
        for res in chain
            ( chainid(res) == source_chainid && resnumber(res) in receptor_selection ) && chainid!(chain, candiate_chainnames[1])
            ( chainid(res) == target_chainid && resnumber(res) in ligand_selection ) && chainid!(chain, candiate_chainnames[2])
        end
    end
    for chain in ppi_complex
        ( chainid(chain) == "A" ) && chainid!(chain, candiate_chainnames[3])
        ( chainid(chain) == "B" ) && chainid!(chain, candiate_chainnames[4])
    end
    if "A" ∉ chainids(collectmodels(ppi_complex)[1])
        for chain in ppi_complex
            ( chainid(chain) == candiate_chainnames[1] ) && chainid!(chain, "A")
        end
        final_receptor_chainname = "A"
    else
        final_receptor_chainname = candiate_chainnames[1]
    end
    if "B" ∉ chainids(collectmodels(ppi_complex)[1])
        for chain in ppi_complex
            ( chainid(chain) == candiate_chainnames[2] ) && chainid!(chain, "B")
        end
        final_ligand_chainname = "B"
    else
        final_ligand_chainname = candiate_chainnames[2]
    end
    return ppi_complex, final_receptor_chainname, final_ligand_chainname
end

function parse_indiv_domain(Structure, chainid, receptor_selection, ligand_selection)
    function merge_frags_to_domain(target_chain, frag_dict, bond_cutoff)
        domain_id = 0
        domain_dict = Dict()
        for (frag_id,residue_ids_list) in frag_dict
            frag = collectresidues(target_chain, res -> resnumber(res) in residue_ids_list)
            if domain_id ==0
                domain_id += 1
                domain_dict[domain_id] = residue_ids_list
            else
                min_doamin_dist = 9999999
                adjacent_doamin_id = ""
                for (domain_no, domain_res_list) in domain_dict
                    domain = collectresidues(target_chain, res -> resnumber(res) in domain_res_list)
                    dist_domain = distance(frag, domain)
                    if dist_domain <= bond_cutoff && dist_domain <= min_doamin_dist
                        min_doamin_dist = dist_domain
                        adjacent_doamin_id = domain_no
                    end
                end
                if min_doamin_dist <= bond_cutoff
                    domain_residue_list = vcat(domain_dict[adjacent_doamin_id], residue_ids_list)
                    domain_dict[adjacent_doamin_id] = domain_residue_list
                else
                    domain_id += 1
                    domain_dict[domain_id] = residue_ids_list
                end
            end
        end
        return domain_dict
    end
    function update_selection(domain_dict, ppi_residue_list, receptor_selection, ligand_selection)
        max_ident_residue_with_receptor = 0
        max_ident_residue_with_ligand = 0
        candiate_receptor_domainid = 1
        candiate_receptor_list = []
        candiate_ligand_domainid = 2
        candiate_ligand_list = []
        for (domain_id, domain_residue_list) in domain_dict
            ident_residue_with_receptor = length(findall(x->x in domain_residue_list, receptor_selection))
            ident_residue_with_ligand = length(findall(x->x in domain_residue_list, ligand_selection))
            if ident_residue_with_receptor > max_ident_residue_with_receptor
                candiate_receptor_domainid = domain_id
                candiate_receptor_list = domain_residue_list
                max_ident_residue_with_receptor = ident_residue_with_receptor
            end
            if ident_residue_with_ligand > max_ident_residue_with_ligand
                candiate_ligand_domainid = domain_id
                candiate_ligand_list = domain_residue_list
                max_ident_residue_with_ligand = ident_residue_with_ligand
            end
        end
        if max_ident_residue_with_receptor > 0 || max_ident_residue_with_ligand > 0
            if candiate_receptor_domainid != candiate_ligand_domainid
                receptor_selection = candiate_receptor_list
                ligand_selection = candiate_ligand_list
            else
                if max_ident_residue_with_receptor >= max_ident_residue_with_ligand
                    receptor_selection = candiate_receptor_list
                    ligand_selection = filter(x -> !(x in receptor_selection), ppi_residue_list)
                else
                    ligand_selection = candiate_ligand_list
                    receptor_selection = filter(x -> !(x in ligand_selection), ppi_residue_list)
                end
            end
        else
            receptor_selection = domain_dict[1]
            ligand_selection = domain_dict[2]
        end
        return receptor_selection, ligand_selection
    end
    ppi_residue_list = []
    target_chain = collectchains(Structure[1], res -> chainid(res) == chainid, heavyatomselector, notwaterselector)
    for residue in target_chain
        push!(ppi_residue_list, resnumber(residue))
    end
    frag_dict = Dict()
    residue_ids_list = []
    frag_id = 0
    for resid in sort!(ppi_residue_list)
        if frag_id == 0
            push!(residue_ids_list, resid)
            frag_id += 1
        end
        if resid != last(residue_ids_list)+1 && frag_id >= 1
            frag_dict[frag_id] = residue_ids_list
            residue_ids_list = []
            frag_id += 1
        else
            push!(residue_ids_list, resid)
        end
    end
    domain_dict = merge_frags_to_domain(target_chain, frag_dict, covalent_bond_cutoff)
    if length(domain_dict) >= 2
        receptor_selection, ligand_selection = update_selection(domain_dict, ppi_residue_list, receptor_selection, ligand_selection)
        ppi_complex = collectchains(collectresidues(target_chain, res in [receptor_selection, ligand_selection], heavyatomselector, notwaterselector))
    elseif length(domain_dict) == 1 && length(frag_dict) >= 2
        domain_dict = merge_frags_to_domain(target_chain, frag_dict, 0)
        receptor_selection, ligand_selection = update_selection(domain_dict, ppi_residue_list, receptor_selection, ligand_selection)
        ppi_complex = collectchains(collectresidues(target_chain, res in [receptor_selection, ligand_selection], heavyatomselector, notwaterselector))
    else
        ppi_complex = collectchains(collectresidues(target_chain, res in [receptor_selection, ligand_selection], heavyatomselector, notwaterselector))
    end
    return ppi_complex, chainid, chainid, receptor_selection, ligand_selection
end

function extract_ppi_complex(Structure::ProteinStructure, source_chainid::String, source_res_beg::Int64, source_res_end::Int64, target_chainid::String, target_res_beg::Int64, target_res_end::Int64, PDBID::String)
    receptor_selection = []
    ligand_selection = []
    for resi in source_res_beg:source_res_end
        push!(receptor_selection, resi)
    end
    for resi in target_res_beg:target_res_end
        push!(ligand_selection, resi)
    end
    selector1(res) = ( chainid(res) == source_chainid && resnumber(res) in receptor_selection )
    selector2(res) = ( chainid(res) == target_chainid && resnumber(res) in ligand_selection )
    combinedselector(res) = selector1(res) || selector2(res)
    if source_chainid == target_chainid
        receptor = collectchains(collectatoms(Structure[1], at -> selector1(at), heavyatomselector, notwaterselector))
        ligand = collectchains(collectatoms(Structure[1], at -> selector2(at), heavyatomselector, notwaterselector))
        if distance(receptor, ligand) >= covalent_bond_cutoff
            ppi_complex = collectchains(collectatoms(Structure[1], at -> combinedselector(at), heavyatomselector, notwaterselector))
        else
            write(error,"Error: the residue range annotation of $PDBID - $source_chainid maybe wrong, please check it, we will judge it by inter-domain distance!\n")
            ppi_complex, source_chainid, target_chainid, receptor_selection, ligand_selection = parse_indiv_domain(Structure, source_chainid, receptor_selection, ligand_selection)
        end
    else
        ppi_complex = collectchains(collectatoms(Structure[1], at -> combinedselector(at), heavyatomselector, notwaterselector))
    end
    return change_chain_names(ppi_complex, source_chainid, target_chainid, receptor_selection, ligand_selection)
end

function parse_pdb(PDBID::String, source_chainid::String, source_res_beg::Int64, source_res_end::Int64, target_chainid::String, target_res_beg::Int64, target_res_end::Int64)
    try
        Structure = retrievepdb(PDBID, dir=PDB_lib)
        ppi_complex, receptor_chainname, ligand_chainname = extract_ppi_complex(Structure, source_chainid, source_res_beg, source_res_end, target_chainid, target_res_beg, target_res_end, PDBID)
        return ppi_complex, receptor_chainname, ligand_chainname
    catch e1
        try    
            file = downloadpdb(PDBID,dir=PDB_lib,format=MMCIF) 
            Structure = read(file, MMCIF)
            ppi_complex, receptor_chainname, ligand_chainname = extract_ppi_complex(Structure, source_chainid, source_res_beg, source_res_end, target_chainid, target_res_beg, target_res_end, PDBID)
            return ppi_complex, receptor_chainname, ligand_chainname
        catch e2
            file = downloadpdb(PDBID,dir=PDB_lib,format=MMTF)
            Structure = read(file, MMTF)
            ppi_complex, receptor_chainname, ligand_chainname = extract_ppi_complex(Structure, source_chainid, source_res_beg, source_res_end, target_chainid, target_res_beg, target_res_end, PDBID)
            return ppi_complex, receptor_chainname, ligand_chainname
        end
    end
end

function process_uniprot(uniprot_dir)
    input_info = input_dir * "/" * uniprot_dir * "/" * uniprot_dir * "_info.csv"
    uniprot_info_filename = input_dir * "/" * uniprot_dir * "/" * uniprot_dir * ".info"
    uniprot_info = open(uniprot_info_filename,"w")
    output_dir = input_dir * "/" * uniprot_dir
    struct_table = DataFrame(CSV.File(input_info; header=["PDBID", "source_uniprot", "source_chainid", "source_res_beg", "source_res_end", "source_res_beg_diff", "source_res_end_diff", "target_uniprot", "target_chainid", "target_res_beg", "target_res_end", "target_res_beg_diff", "target_res_end_diff"]))
    sequence_regions = Dict()
    sequence_regions_id = 0
    non_proximity_list = []
    struct_table.structure = missings(Vector{Chain}, nrow(struct_table))
    struct_table.proximity_score = missings(Float64, nrow(struct_table))
    struct_table.adjacent_atoms = missings(Int, nrow(struct_table))
    struct_table.receptor_chainname = missings(String, nrow(struct_table))
    struct_table.ligand_chainname = missings(String, nrow(struct_table))
    struct_table.region_dismatch = missings(String, nrow(struct_table))
    struct_table.source_chain_length = missings(Int, nrow(struct_table))
    struct_table.residue_num_error = missings(String, nrow(struct_table))
    for (row_index, row) in enumerate(eachrow(struct_table))
        PDBID::String = row.PDBID
        source_uniprot = row.source_uniprot
        source_chainid::String = string(row.source_chainid)
        target_uniprot = row.target_uniprot
        target_chainid::String = string(row.target_chainid)
        try
            source_res_beg = convert(Int, row.source_res_beg)
            source_res_end = convert(Int, row.source_res_end)
            struct_table.source_chain_length[row_index] = ( source_res_end - source_res_beg )
            source_res_beg_diff = convert(Int, row.source_res_beg_diff)
            source_res_end_diff = convert(Int, row.source_res_end_diff)
            target_res_beg = convert(Int, row.target_res_beg)
            target_res_end = convert(Int, row.target_res_end)
            target_res_beg_diff = convert(Int, row.target_res_beg_diff)
            target_res_end_diff = convert(Int, row.target_res_end_diff)
            real_begin_res = source_res_beg + source_res_beg_diff
            real_end_res = source_res_end + source_res_end_diff
            struct_table.residue_num_error[row_index] = "no"
        catch e1
            struct_table.residue_num_error[row_index] = "yes"
            write(error,"Error: resdue number parsing error $e1, check $source_uniprot, $source_chainid, $PDBID !\n")
        end
        # check the sequence region match
        struct_table.region_dismatch[row_index] = "no"
        ( source_res_beg_diff != source_res_end_diff ) && begin 
            write(error,"Error: chain length bewteen PDB with Uniprot isn't matched, check $source_uniprot, $source_chainid, $PDBID !\n")
            struct_table.region_dismatch[row_index] = "yes"
        end
        ( target_res_beg_diff != target_res_end_diff ) && begin
            write(error,"Error: chain length bewteen PDB with Uniprot isn't matched, check $target_uniprot, $target_chainid, $PDBID !\n")
            struct_table.region_dismatch[row_index] = "yes"
        end
        # process structure
        struct_table.structure[row_index], struct_table.receptor_chainname[row_index], struct_table.ligand_chainname[row_index] = parse_pdb(PDBID, source_chainid, source_res_beg, source_res_end, target_chainid, target_res_beg, target_res_end)
        struct_table.adjacent_atoms[row_index], struct_table.proximity_score[row_index] = check_protein_proximity(struct_table.structure[row_index]) 
    end
    struct_table = struct_table[struct_table.adjacent_atoms .>=min_adjacent_atoms, :]
    struct_table = struct_table[struct_table.residue_num_error .=="no", :]
    if nrow(struct_table) == 0
        write(error,"Warning: no $uniprot_dir ppi complex close enough! \n")
    else
        struct_table.receptor_region_id = missings(Int, nrow(struct_table))
        for (row_index, row) in enumerate(eachrow(struct_table))
            source_res_beg = convert(Int, row.source_res_beg)
            source_res_end = convert(Int, row.source_res_end)
            source_res_beg_diff = convert(Int, row.source_res_beg_diff)
            source_res_end_diff = convert(Int, row.source_res_end_diff)
            real_begin_res = source_res_beg + source_res_beg_diff
            real_end_res = source_res_end + source_res_end_diff
            # diss sequence region 
            if sequence_regions_id == 0
                sequence_regions_id += 1
                sequence_regions[sequence_regions_id] = [real_begin_res, real_end_res]
                struct_table.receptor_region_id[row_index] = sequence_regions_id
            else
                non_matched_region_num = 0
                for (region_id,sequence_region) in sequence_regions
                    if real_begin_res > sequence_region[2] || real_end_res < sequence_region[1]
                        non_matched_region_num += 1
                    else
                        new_real_begin = min(real_begin_res, real_end_res, sequence_region[1], sequence_region[2])
                        new_real_end = max(real_begin_res, real_end_res, sequence_region[1], sequence_region[2])
                        sequence_region = [new_real_begin, new_real_end]
                        struct_table.receptor_region_id[row_index] = region_id
                    end
                end
                non_matched_region_num == sequence_regions_id && begin 
                    sequence_regions_id += 1
                    sequence_regions[sequence_regions_id] = [real_begin_res, real_end_res]
                    struct_table.receptor_region_id[row_index] = sequence_regions_id
                end
            end
        end
        for (region_id,sequence_region) in sequence_regions
            write(uniprot_info, "$region_id,$sequence_region\n")
        end
        for sequence_region in 1:sequence_regions_id
            struct_subset = struct_table[struct_table.receptor_region_id .== sequence_region, :]
            sort!(struct_subset, rev = true, [:source_chain_length])
            matched_ppi = struct_subset[struct_subset.region_dismatch .== "no", :]
            if nrow(matched_ppi) == 0
                ref_ppi = struct_subset[1, :]
                write(error,"Error: sequence region $sequence_region of $uniprot_dir is dismatch !\n")
                # continue
            else
                ref_ppi = matched_ppi[1, :]
            end
            ref_receptor_chain = ref_ppi.receptor_chainname
            ref_source_res_beg = convert(Int, ref_ppi.source_res_beg)
            ref_source_res_end = convert(Int, ref_ppi.source_res_end)
            ref_source_res_diff = convert(Int, ref_ppi.source_res_beg_diff)
            if nrow(struct_subset) == 1
                row = struct_subset[1, :]
                source_uniprot_id = row.source_uniprot
                region_id = string(row.receptor_region_id)
                target_uniprot_id = row.target_uniprot
                PDBID = row.PDBID
                score = string(round.(row.proximity_score, digits=3))
                struct_name = source_uniprot_id * "_" * region_id * "_" * target_uniprot_id * "_" * PDBID * "_" * score
                write_struct(row.structure, output_dir, struct_name)
            elseif nrow(struct_subset) > 1
                for row in eachrow(struct_subset)
                    source_uniprot_id = row.source_uniprot
                    region_id = string(row.receptor_region_id)
                    target_uniprot_id = row.target_uniprot
                    PDBID = row.PDBID
                    score = string(round.(row.proximity_score, digits=3))
                    receptor_chain = row.receptor_chainname
                    struct_name = source_uniprot_id * "_" * region_id * "_" * target_uniprot_id * "_" * PDBID * "_" * score
                    source_res_diff = convert(Int, row.source_res_beg_diff)
                    source_res_beg = convert(Int, row.source_res_beg)
                    source_res_end = convert(Int, row.source_res_end)
                    region_dismatch = string(row.region_dismatch)
                    if source_res_diff == ref_source_res_diff && ref_receptor_chain == receptor_chain
                        common_res_begin = max(source_res_beg, ref_source_res_beg)
                        common_res_end = min(ref_source_res_end, source_res_end)
                        superimpose!(row.structure, ref_ppi.structure, res -> chainid(res) == ref_receptor_chain && common_res_begin < resnumber(res) < common_res_end )
                    else
                        write(error,"Align Error: different chain-name of reorder-residue-id bewteen reference structure ($ref_receptor_chain) with query structure ($receptor_chain), $output_dir/$struct_name won't align!\n")
                    end
                    write_struct(row.structure, output_dir, struct_name)
                end
            end
        end
        CSV.write(input_dir * "/" * uniprot_dir * "/" * uniprot_dir * "_data.csv", struct_table)
    end
end

print("What's your ppi info path?\n")
const input_dir = readline(stdin)
const PDB_lib = "PDB_Library"
const min_dist = 4.0  # 3.30 is hydrogen bond cutoff, 4.0 is van der Waals bonds cutoff.
const min_adjacent_atoms = 10 # How many residues are close together between frags will be identified as a domain.
const sigma = 2*1.094 # 1.094 * 2, 1.094 is length of C-H bond.
const epsilon = 1
const covalent_bond_cutoff = 2.2
const struct_mode = "chain" # chain or residue
global error = open("Parse_PPI_Dimer_log.error","w")
println("Strating Julia with ",Threads.nthreads()," threads ....")
Threads.@threads for uniprot_dir in readdir(input_dir)
    if isfile(input_dir * "/" * uniprot_dir * "/" * uniprot_dir * "_data.csv")
        continue
    else
        process_uniprot(uniprot_dir)
    end
end


