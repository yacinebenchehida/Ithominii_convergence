using DataFrames, CSV

# Capture command-line arguments
function main()
    args = ARGS
    if length(args) < 2
        println("Usage: julia script.jl species_list path_to_data")
        return
    end
    
    # Parse species list and data path from command-line arguments
    species_list = split(args[1], ",")
    path_2_data = args[2]
    
    # Combine data for all species
    Merged_dt = Combine_data(species_list, path_2_data)
    
    # Determine the majority ancestral state
    Alternative_ancestral_state = Get_majority(Merged_dt)
    
    # Set alternative ancestral alleles where 'nan' exists
    Alternative_ancestral_output = Set_alternative_ancestral_nan(path_2_data, Alternative_ancestral_state)
    
    # Print final output
    println(Alternative_ancestral_state)
end

# Function to combine data from species-specific files
function Combine_data(species_list, path_2_data)
    all_species_data = DataFrame(Chromosome = String[], Position = Int[], 
                                  Allele1_frequency = String[], Allele2_frequency = String[], 
                                  Allele1 = String[], frequency1 = Float64[], 
                                  Allele2 = String[], frequency2 = Float64[], 
                                  Species = String[])  # Initialize with specific columns
    
    for species in species_list
        # Read each species-specific file with .frq extension
        species_files = filter(x -> occursin(r"^$species.*\.frq$", x), readdir(path_2_data))
        
        for file in species_files
            # Read file and assign column names manually
            file_path = joinpath(path_2_data, file)
            frequencies = CSV.read(file_path, DataFrame, header = false)
            rename!(frequencies, [:Chromosome, :Position, :V3, :V4, :Allele1_frequency, :Allele2_frequency])
            
            # Split allele frequency columns (Allele1:frequency1, Allele2:frequency2)
            frequencies[!, :Allele1] = first.(split.(frequencies[!, :Allele1_frequency], ":"))
            frequencies[!, :frequency1] = parse.(Float64, last.(split.(frequencies[!, :Allele1_frequency], ":")))
            frequencies[!, :Allele2] = first.(split.(frequencies[!, :Allele2_frequency], ":"))
            frequencies[!, :frequency2] = parse.(Float64, last.(split.(frequencies[!, :Allele2_frequency], ":")))
            
            # Remove unnecessary columns and add species identifier
            select!(frequencies, Not([:V3, :V4]))
            frequencies[!, :Species] .= species
            
            # Filter out rows where frequency1 is NaN
            filter!(row -> !isnan(row.frequency1), frequencies)
            
            # Append data to all_species_data
            append!(all_species_data, frequencies)  # Combine into DataFrame
        end
    end
    
    return all_species_data  # Return combined DataFrame
end

# Function to get the majority allele by summing frequencies for each SNP
function Get_majority(complete_df::DataFrame)
    SNPs_info = DataFrame(Position = String[], MajorityAllele = String[])
    
    for pos in unique(complete_df.Position)
        current_SNP = complete_df[complete_df.Position .== pos, :]
        allele1_sum = sum(current_SNP.frequency1)
        allele2_sum = sum(current_SNP.frequency2)
        
        majority_allele = allele1_sum > allele2_sum ? unique(current_SNP.Allele1)[1] : unique(current_SNP.Allele2)[1]
        
        push!(SNPs_info, (Position = string(pos), MajorityAllele = majority_allele))
    end
    
    return SNPs_info
end

# Function to set alternative ancestral alleles where 'nan' exists
function Set_alternative_ancestral_nan(path_2_data, Alternative_ancestral_state)
    ancestral_file = joinpath(path_2_data, "ancestral.alleles")
    current_ancestral = CSV.read(ancestral_file, DataFrame, header = false)
    rename!(current_ancestral, [:Position, :AncestralAllele])
    
    for row in eachrow(current_ancestral)
        if row.AncestralAllele == "nan"
            alternative_state = Alternative_ancestral_state[Alternative_ancestral_state.Position .== row.Position, :MajorityAllele]
            if !isempty(alternative_state)
                row.AncestralAllele = alternative_state[1]  # Update with found allele
            end
        end
    end
    
    # Return the updated DataFrame with the new ancestral alleles
    return current_ancestral
end

# Run the main function
main()
