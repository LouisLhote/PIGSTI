#!/usr/bin/env python3
"""
PIGSTI Pathogen Database Creator
Complete solution: Creates conda environment, downloads genomes, builds indices, creates spreadsheet
"""

import os
import sys
import csv
import subprocess
import tempfile
import zipfile
import shutil
from pathlib import Path
import argparse
import time

# Known accession numbers for problematic viruses
VIRUS_ACCESSIONS = {
    "Monkeypox virus": "NC_003310.1",
    "Cowpox virus": "NC_003663.1", 
    "Vaccinia virus": "NC_006998.1",
    "Variola virus": "NC_001611.1",
    "ORF virus": "NC_001899.1",
    "Sheeppox virus": "NC_002731.1",
    "Goatpox virus": "NC_004003.1",
    "Lumpy skin disease virus": "NC_003027.1",
    "Swinepox virus": "NC_003389.1",
    "Myxoma virus": "NC_001132.2",
    "Molluscum contagiosum virus": "NC_001731.1",
    "Tanapox virus": "NC_009888.1",
    "Yaba monkey tumor virus": "NC_005179.1",
    "Fowlpox virus": "NC_002188.1",
    "Pigeonpox virus": "NC_005309.1",
    "Canarypox virus": "NC_005309.1",
    "Marek's disease virus": "NC_002229.3",
    "African swine fever virus": "NC_001659.1",
}

# Comprehensive pathogen list with taxIDs - Veterinary and Zoonotic Pathogens
PATHOGENS = [
    # === Poxviruses (Comprehensive) ===
    # Orthopoxviruses
    ("10255", "Variola virus", "Variola virus", "virus"),
    ("10244", "Monkeypox virus", "Monkeypox virus", "virus"),
    ("10243", "Cowpox virus", "Cowpox virus", "virus"),
    ("10245", "Vaccinia virus", "Vaccinia virus", "virus"),
    ("10279", "Camelpox virus", "Camelpox virus", "virus"),
    
    # Parapoxviruses
    ("10258", "ORF virus", "ORF virus", "virus"),
    ("10266", "Sheeppox virus", "Sheeppox virus 17077-99", "virus"),
    ("186805", "Goatpox virus", "Goatpox virus Pellor", "virus"),
    ("59509", "Lumpy skin disease virus", "Lumpy skin disease virus NI-2490", "virus"),
    
    # Suipoxviruses
    ("10275", "Swinepox virus", "Swinepox virus", "virus"),
    
    # Leporipoxviruses
    ("10273", "Myxoma virus", "Myxoma virus", "virus"),
    
    # Molluscipoxviruses
    ("10280", "Molluscum contagiosum virus", "Molluscum contagiosum virus", "virus"),
    
    # Yatapoxviruses
    ("10277", "Tanapox virus", "Tanapox virus", "virus"),
    ("10278", "Yaba monkey tumor virus", "Yaba monkey tumor virus", "virus"),
    
    # Avipoxviruses (Bird poxviruses)
    ("10263", "Fowlpox virus", "Fowlpox virus", "virus"),
    ("10264", "Pigeonpox virus", "Pigeonpox virus", "virus"),
    ("10260", "Canarypox virus", "Canarypox virus", "virus"),
    ("10262", "Flamingopox virus", "Flamingopox virus", "virus"),
    ("10265", "Juncopox virus", "Juncopox virus", "virus"),
    ("10266", "Mynampox virus", "Mynampox virus", "virus"),
    ("10267", "Penguinpox virus", "Penguinpox virus", "virus"),
    ("10268", "Psittacinepox virus", "Psittacinepox virus", "virus"),
    ("10269", "Quailpox virus", "Quailpox virus", "virus"),
    ("10270", "Sparrowpox virus", "Sparrowpox virus", "virus"),
    ("10271", "Starlingpox virus", "Starlingpox virus", "virus"),
    ("10272", "Turkeypox virus", "Turkeypox virus", "virus"),
    
    # Cervidpoxviruses
    ("305674", "Muledeerpox virus", "Muledeerpox virus", "virus"),
    
    
    # === DNA Viruses (Veterinary) ===
    ("10282", "Bovine herpesvirus 1", "BHV1", "virus"),
    ("10920", "Bovine papillomavirus 1", "BPV1", "virus"),
    ("115514", "Porcine circovirus 2", "PCV2", "virus"),
    ("176531", "Canine parvovirus 2", "Canine parvovirus", "virus"),
    ("10344", "Fowl aviadenovirus A", "FAdV_A", "virus"),
    
    # === Important Veterinary Viruses ===
    ("10390", "Marek's disease virus", "Marek's disease virus", "virus"),
    ("10496", "African swine fever virus", "African swine fever virus", "virus"),
    ("10792", "Canine parvovirus", "Canine parvovirus", "virus"),
    ("11232", "Canine distemper virus", "Canine distemper virus", "virus"),
    ("11292", "Rabies virus", "Rabies virus", "virus"),
    ("11976", "Feline calicivirus", "Feline calicivirus", "virus"),
    ("10320", "Feline herpesvirus 1", "Feline herpesvirus 1", "virus"),
    ("11676", "Feline immunodeficiency virus", "Feline immunodeficiency virus", "virus"),
    ("10318", "Equine herpesvirus 1", "Equine herpesvirus 1", "virus"),
    ("10320", "Equine herpesvirus 4", "Equine herpesvirus 4", "virus"),
    ("11777", "Equine infectious anemia virus", "EIAV", "virus"),
    ("11060", "Bluetongue virus", "Bluetongue virus", "virus"),
    ("11709", "Maedi-Visna virus", "Maedi Visna virus", "virus"),
    ("11657", "Caprine arthritis encephalitis virus", "CAEV", "virus"),
    ("31604", "Peste des petits ruminants virus", "PPRV", "virus"),
    ("11062", "African horse sickness virus", "African horse sickness virus", "virus"),
    
    # === Core Bacteria (Original) ===
    ("1504", "Clostridium septicum", "Clostridium septicum", "genome"),
    ("29459", "Brucella melitensis", "Brucella melitensis", "genome"),
    ("238", "Brucella abortus", "Brucella abortus", "genome"),
    ("236", "Brucella ovis", "Brucella ovis", "genome"),
    ("29461", "Brucella suis", "Brucella suis", "genome"),
    ("632", "Yersinia pestis", "Yersinia pestis", "genome"),
    ("777", "Coxiella burnetii", "Coxiella burnetii", "genome"),
    ("160", "Treponema pallidum", "Treponema pallidum", "genome"),
    ("1773", "Mycobacterium tuberculosis", "Mycobacterium tuberculosis", "genome"),
    ("630", "Yersinia enterocolitica", "Yersinia enterocolitica", "genome"),
    ("2102", "Mycoplasma mycoides", "Mycoplasma mycoides", "genome"),
    ("126793", "Plasmodium vivax", "Plasmodium vivax", "genome"),
    ("562", "Escherichia coli", "Escherichia coli", "genome"),
    ("1280", "Staphylococcus aureus", "Staphylococcus aureus", "genome"),
    ("1313", "Streptococcus pneumoniae", "Streptococcus pneumoniae", "genome"),
    ("1423", "Bacillus anthracis", "Bacillus anthracis", "genome"),
    ("28901", "Salmonella enterica", "Salmonella enterica", "genome"),
    ("290", "Pseudomonas aeruginosa", "Pseudomonas aeruginosa", "genome"),
    ("1763", "Mycobacterium leprae", "Mycobacterium leprae", "genome"),
    ("1764", "Mycobacterium bovis", "Mycobacterium bovis", "genome"),
    ("1765", "Mycobacterium avium", "Mycobacterium avium", "genome"),
    
    # === Veterinary Bacteria ===
    ("338733", "Staphylococcus pseudintermedius", "Staphylococcus pseudintermedius", "genome"),
    ("1639", "Listeria monocytogenes", "Listeria monocytogenes", "genome"),
    ("1305", "Streptococcus suis", "Streptococcus suis", "genome"),
    ("747", "Pasteurella multocida", "Pasteurella multocida", "genome"),
    ("100", "Chlamydia abortus", "Chlamydia abortus", "genome"),
    ("828", "Chlamydia pecorum", "Chlamydia pecorum", "genome"),
    ("835", "Chlamydia psittaci", "Chlamydia psittaci", "genome"),
    ("724909", "Mycoplasma bovis", "Mycoplasma bovis", "genome"),
    ("2191", "Mycoplasma hyopneumoniae", "Mycoplasma hyopneumoniae", "genome"),
    ("127", "Anaplasma phagocytophilum", "Anaplasma phagocytophilum", "genome"),
    ("945", "Anaplasma marginale", "Anaplasma marginale", "genome"),
    ("29654", "Francisella tularensis", "Francisella tularensis", "genome"),
    ("194", "Legionella pneumophila", "Legionella pneumophila", "genome"),
    ("581", "Haemophilus influenzae", "Haemophilus influenzae", "genome"),
    
    # === Tick-borne Bacteria ===
    ("139", "Borrelia burgdorferi", "Borrelia burgdorferi", "genome"),
    ("29465", "Borrelia recurrentis", "Borrelia recurrentis", "genome"),
    ("783", "Rickettsia rickettsii", "Rickettsia rickettsii", "genome"),
    ("781", "Rickettsia conorii", "Rickettsia conorii", "genome"),
    ("782", "Rickettsia prowazekii", "Rickettsia prowazekii", "genome"),
    ("784", "Orientia tsutsugamushi", "Orientia tsutsugamushi", "genome"),
    
    # === Protozoa ===
    ("5833", "Plasmodium falciparum", "Plasmodium falciparum", "genome"),
    ("5855", "Plasmodium vivax", "Plasmodium vivax", "genome"),
    ("5850", "Plasmodium knowlesi", "Plasmodium knowlesi", "genome"),
    ("5865", "Babesia microti", "Babesia microti", "genome"),
    ("5864", "Babesia bovis", "Babesia bovis", "genome"),
    ("5875", "Theileria parva", "Theileria parva", "genome"),
    ("5661", "Leishmania donovani", "Leishmania donovani", "genome"),
    ("5671", "Leishmania infantum", "Leishmania infantum", "genome"),
    ("5691", "Trypanosoma brucei", "Trypanosoma brucei", "genome"),
    ("5693", "Trypanosoma cruzi", "Trypanosoma cruzi", "genome"),
    ("5811", "Toxoplasma gondii", "Toxoplasma gondii", "genome"),
    ("5806", "Cryptosporidium parvum", "Cryptosporidium parvum", "genome"),
    ("5741", "Giardia lamblia", "Giardia lamblia", "genome"),
    ("5759", "Entamoeba histolytica", "Entamoeba histolytica", "genome"),
    
    # === Helminths ===
    ("6210", "Echinococcus granulosus", "Echinococcus granulosus", "genome"),
    ("6211", "Echinococcus multilocularis", "Echinococcus multilocularis", "genome"),
    ("6204", "Taenia solium", "Taenia solium", "genome"),
    ("6208", "Taenia saginata", "Taenia saginata", "genome"),
    ("6334", "Trichinella spiralis", "Trichinella spiralis", "genome"),
    ("6253", "Ascaris suum", "Ascaris suum", "genome"),
    ("29170", "Ancylostoma caninum", "Ancylostoma caninum", "genome"),
    ("6248", "Strongyloides stercoralis", "Strongyloides stercoralis", "genome"),
    ("6183", "Schistosoma mansoni", "Schistosoma mansoni", "genome"),
    ("6185", "Schistosoma haematobium", "Schistosoma haematobium", "genome"),
    ("6182", "Schistosoma japonicum", "Schistosoma japonicum", "genome"),
    ("6279", "Haemonchus contortus", "Haemonchus contortus", "genome"),
    ("45464", "Teladorsagia circumcincta", "Teladorsagia circumcincta", "genome"),
    ("6192", "Fasciola hepatica", "Fasciola hepatica", "genome"),
    
    # === Fungi ===
    ("5037", "Histoplasma capsulatum", "Histoplasma capsulatum", "genome"),
    ("5501", "Coccidioides immitis", "Coccidioides immitis", "genome"),
    ("5207", "Cryptococcus neoformans", "Cryptococcus neoformans", "genome"),
    ("746128", "Aspergillus fumigatus", "Aspergillus fumigatus", "genome"),
    
    # === Companion Animals ===
    ("38323", "Bartonella henselae", "Bartonella henselae", "genome"),
    ("28188", "Capnocytophaga canimorsus", "Capnocytophaga canimorsus", "genome"),
    ("1341", "Streptococcus canis", "Streptococcus canis", "genome"),
    ("32257", "Mycoplasma haemofelis", "Mycoplasma haemofelis", "genome"),
    
    # === Small Ruminants ===
    ("1719", "Corynebacterium pseudotuberculosis", "Corynebacterium pseudotuberculosis", "genome"),
    ("75985", "Mannheimia haemolytica", "Mannheimia haemolytica", "genome"),
    ("75686", "Pasteurella trehalosi", "Pasteurella trehalosi", "genome"),
    ("2110", "Mycoplasma agalactiae", "Mycoplasma agalactiae", "genome"),
    ("45362", "Mycoplasma capricolum subsp. capripneumoniae", "Mycoplasma capripneumoniae", "genome"),
    ("83554", "Chlamydophila abortus", "Chlamydophila abortus", "genome"),
    
    # === Horses ===
    ("13373", "Burkholderia mallei", "Burkholderia mallei", "genome"),
    ("28450", "Burkholderia pseudomallei", "Burkholderia pseudomallei", "genome"),
    ("1336", "Streptococcus equi subsp. equi", "Streptococcus equi equi", "genome"),
    ("1335", "Streptococcus equi subsp. zooepidemicus", "Streptococcus equi zooepidemicus", "genome"),
    ("1644", "Rhodococcus equi", "Rhodococcus equi", "genome"),
]

def run_command(cmd, check=True, capture_output=True, env_name=None):
    """Run a command and return the result"""
    try:
        print(f"Running: {' '.join(cmd)}")
        
        # If environment is specified, activate it first
        if env_name:
            # Use conda run to execute command in specific environment
            full_cmd = ['conda', 'run', '-n', env_name] + cmd
        else:
            full_cmd = cmd
            
        result = subprocess.run(full_cmd, check=check, capture_output=capture_output, text=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {' '.join(cmd)}")
        print(f"Error: {e.stderr}")
        return None

def check_mamba():
    """Check if mamba is available"""
    result = run_command(['mamba', '--version'], check=False)
    if result and result.returncode == 0:
        print("✅ Mamba is available")
        return 'mamba'
    else:
        print("❌ Mamba not found, trying conda...")
        result = run_command(['conda', '--version'], check=False)
        if result and result.returncode == 0:
            print("✅ Conda is available")
            return 'conda'
        else:
            print("❌ Neither mamba nor conda found!")
            return False

def create_environment(env_name="pigsti_pathogens"):
    """Create conda environment with required packages"""
    print(f"\n🔧 Creating conda environment: {env_name}")
    
    # Check if mamba/conda is available
    package_manager = check_mamba()
    if not package_manager:
        return False
    
    # Create environment
    if package_manager == 'mamba':
        cmd = ['mamba', 'create', '-n', env_name, '-y']
        print("🐍 Using mamba for environment creation")
    else:
        cmd = ['conda', 'create', '-n', env_name, '-y']
        print("🐍 Using conda for environment creation")
    
    # Add packages
    packages = [
        'python=3.9',
        'ncbi-datasets-cli',
        'entrez-direct',  # For efetch, esearch, etc.
        'bwa',
        'curl',
        'unzip',
        'samtools',
        'pandas',
        'requests'
    ]
    
    cmd.extend(packages)
    
    result = run_command(cmd, check=False)
    if result and result.returncode == 0:
        print(f"✅ Environment '{env_name}' created successfully")
        return True
    else:
        print(f"❌ Failed to create environment '{env_name}'")
        return False

def activate_environment(env_name="pigsti_pathogens"):
    """Get the command to activate the environment"""
    return f"conda activate {env_name}"

def sanitize_name(name):
    """Sanitize pathogen name for filesystem use"""
    return "".join(c for c in name if c.isalnum() or c in (' ', '-', '_')).rstrip().replace(' ', '_')

def download_genome_efetch(taxid, name, tmpdir, env_name="pigsti_pathogens"):
    """Download genome using NCBI efetch (alternative method)"""
    print(f"📥 Trying efetch method for {name} (taxID: {taxid})...")
    
    fasta_path = os.path.join(tmpdir, f"{taxid}.fa")
    
    # Try efetch to get FASTA directly
    cmd = ["efetch", "-db", "nuccore", "-id", str(taxid), "-format", "fasta", "-mode", "text"]
    
    result = run_command(cmd, check=False, env_name=env_name)
    if result and result.returncode == 0 and result.stdout.strip():
        with open(fasta_path, 'w') as f:
            f.write(result.stdout)
        
        if os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 1000:  # At least 1KB
            print(f"✅ Downloaded via efetch: {name}")
            return fasta_path
    
    print(f"❌ efetch failed for {name}")
    return None

def download_genome_curl(taxid, name, tmpdir, env_name="pigsti_pathogens"):
    """Download genome using curl from NCBI FTP"""
    print(f"📥 Trying curl method for {name} (taxID: {taxid})...")
    
    fasta_path = os.path.join(tmpdir, f"{taxid}.fa")
    
    # Try different NCBI FTP URLs
    urls = [
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={taxid}&rettype=fasta&retmode=text",
        f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id={taxid}",
    ]
    
    for url in urls:
        try:
            cmd = ["curl", "-s", "-L", url, "-o", fasta_path]
            result = run_command(cmd, check=False, env_name=env_name)
            
            if result and result.returncode == 0 and os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 1000:
                print(f"✅ Downloaded via curl: {name}")
                return fasta_path
            else:
                print(f"❌ curl returned code {result.returncode if result else 'None'} for {name}")
        except Exception as e:
            print(f"❌ curl failed for {name}: {e}")
            continue
    
    print(f"❌ All curl methods failed for {name}")
    return None

def download_genome_specific_accession(name, accession, tmpdir, env_name="pigsti_pathogens"):
    """Download genome using specific accession number"""
    print(f"📥 Trying specific accession {accession} for {name}...")
    
    fasta_path = os.path.join(tmpdir, f"{accession}.fa")
    
    # Try efetch with specific accession
    cmd = ["efetch", "-db", "nuccore", "-id", accession, "-format", "fasta", "-mode", "text"]
    result = run_command(cmd, check=False, env_name=env_name)
    
    if result and result.returncode == 0 and result.stdout.strip():
        with open(fasta_path, 'w') as f:
            f.write(result.stdout)
        
        if os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 1000:
            print(f"✅ Downloaded via accession {accession}: {name}")
            return fasta_path
    
    # Try curl with specific accession
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    cmd = ["curl", "-s", "-L", url, "-o", fasta_path]
    result = run_command(cmd, check=False, env_name=env_name)
    
    if result and result.returncode == 0 and os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 1000:
        print(f"✅ Downloaded via curl accession {accession}: {name}")
        return fasta_path
    
    print(f"❌ Failed to download {name} with accession {accession}")
    return None

def download_genome(taxid, name, hops_name, dataset_type, outdir, env_name="pigsti_pathogens"):
    """Download genome using multiple methods"""
    print(f"\n📥 Downloading {name} (taxID: {taxid})...")
    
    # Create temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        zip_path = os.path.join(tmpdir, f"{taxid}_{dataset_type}.zip")
        source_fasta = None
        
        # Method 1: Try NCBI datasets first
        if dataset_type == "virus":
            cmd = ["datasets", "download", "virus", "taxon", str(taxid), "--reference", "--filename", zip_path]
        else:
            cmd = ["datasets", "download", "genome", "taxon", str(taxid), "--reference", "--filename", zip_path]
        
        result = run_command(cmd, check=False, env_name=env_name)
        if result and result.returncode == 0 and os.path.exists(zip_path):
            print(f"✅ Downloaded via datasets: {name}")
            # Extract and find FASTA
            extract_dir = os.path.join(tmpdir, "extracted")
            os.makedirs(extract_dir, exist_ok=True)
            
            try:
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    zip_ref.extractall(extract_dir)
                
                # Find FASTA file
                fasta_files = []
                for root, dirs, files in os.walk(extract_dir):
                    for file in files:
                        if file.endswith(('.fna', '.fa', '.fasta')):
                            fasta_files.append(os.path.join(root, file))
                
                if fasta_files:
                    source_fasta = fasta_files[0]
            except Exception as e:
                print(f"❌ Error extracting {name}: {e}")
        
        # Method 2: If datasets failed, try specific accession if available
        if not source_fasta and name in VIRUS_ACCESSIONS:
            accession = VIRUS_ACCESSIONS[name]
            source_fasta = download_genome_specific_accession(name, accession, tmpdir, env_name)
        
        # Method 3: If specific accession failed, try efetch
        if not source_fasta:
            source_fasta = download_genome_efetch(taxid, name, tmpdir, env_name)
        
        # Method 4: If efetch failed, try curl
        if not source_fasta:
            source_fasta = download_genome_curl(taxid, name, tmpdir, env_name)
        
        # Method 4: If all methods failed, try alternative taxIDs for viruses
        if not source_fasta and dataset_type == "virus":
            print(f"🔄 Trying alternative methods for virus {name}...")
            # For some viruses, try downloading representative sequences
            alt_cmd = ["datasets", "download", "virus", "taxon", str(taxid), "--filename", zip_path]
            result = run_command(alt_cmd, check=False, env_name=env_name)
            if result and result.returncode == 0 and os.path.exists(zip_path):
                try:
                    extract_dir = os.path.join(tmpdir, "extracted")
                    os.makedirs(extract_dir, exist_ok=True)
                    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                        zip_ref.extractall(extract_dir)
                    
                    fasta_files = []
                    for root, dirs, files in os.walk(extract_dir):
                        for file in files:
                            if file.endswith(('.fna', '.fa', '.fasta')):
                                fasta_files.append(os.path.join(root, file))
                    
                    if fasta_files:
                        source_fasta = fasta_files[0]
                        print(f"✅ Downloaded via alternative method: {name}")
                except Exception as e:
                    print(f"❌ Alternative method failed for {name}: {e}")
        
        if not source_fasta:
            print(f"❌ All download methods failed for {name}")
            return None
        
        # Create destination directory and copy file
        safe_name = sanitize_name(name)
        dest_dir = os.path.join(outdir, safe_name)
        os.makedirs(dest_dir, exist_ok=True)
        
        dest_fasta = os.path.join(dest_dir, f"{safe_name}.fa")
        
        try:
            # Copy and clean FASTA file
            with open(source_fasta, 'r') as infile, open(dest_fasta, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        outfile.write(line)
                    else:
                        # Remove carriage returns
                        outfile.write(line.rstrip('\r\n') + '\n')
            
            print(f"✅ Saved: {dest_fasta}")
            
            # Create BWA index
            print(f"🔧 Building BWA index for {name}...")
            result = run_command(['bwa', 'index', dest_fasta], check=False, env_name=env_name)
            if result and result.returncode == 0:
                print(f"✅ BWA index created for {name}")
            else:
                print(f"❌ Failed to build BWA index for {name}")
            
            return dest_fasta
            
        except Exception as e:
            print(f"❌ Error processing FASTA for {name}: {e}")
            return None

def create_pathogen_spreadsheet(pathogens, outdir, output_csv):
    """Create the pathogen spreadsheet CSV"""
    print(f"\n📋 Creating pathogen spreadsheet: {output_csv}")
    
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['taxID', 'Krakenuniq name', 'Hops name', 'bwa index', 'min_escore', 'min_reads'])
        
        for taxid, kraken_name, hops_name, dataset_type in pathogens:
            safe_name = sanitize_name(kraken_name)
            bwa_index = os.path.join(outdir, safe_name, f"{safe_name}.fa")
            
            # Set default values
            min_escore = 5
            min_reads = 50
            
            # Adjust thresholds based on pathogen type
            if dataset_type == "virus":
                min_escore = 0  # Viruses often have lower E-scores
                min_reads = 1
            elif any(keyword in kraken_name.lower() for keyword in ['mycobacterium', 'tuberculosis']):
                min_escore = 5
                min_reads = 100  # Higher threshold for Mycobacterium
            elif any(keyword in kraken_name.lower() for keyword in ['brucella', 'yersinia', 'coxiella']):
                min_escore = 2
                min_reads = 20
            
            writer.writerow([taxid, kraken_name, hops_name, bwa_index, min_escore, min_reads])

def main():
    parser = argparse.ArgumentParser(description='Create PIGSTI pathogen database')
    parser.add_argument('--env-name', default='pigsti_pathogens',
                       help='Name of conda environment to create')
    parser.add_argument('--outdir', default='/raid_md0/Reference_Genomes/PIGSTI_pathogens',
                       help='Output directory for pathogen genomes')
    parser.add_argument('--output-csv', default='config/PIGSTI_pathogens.csv',
                       help='Output CSV file for pathogen spreadsheet')
    parser.add_argument('--skip-env', action='store_true',
                       help='Skip environment creation (assume tools are available)')
    parser.add_argument('--skip-download', action='store_true',
                       help='Skip genome download, only create CSV')
    
    args = parser.parse_args()
    
    print("🚀 PIGSTI Pathogen Database Creator")
    print("=" * 50)
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)
    
    # Create conda environment if requested
    if not args.skip_env:
        if not create_environment(args.env_name):
            print("❌ Failed to create environment. Exiting.")
            sys.exit(1)
        
        print(f"\n📝 To activate the environment, run:")
        print(f"   {activate_environment(args.env_name)}")
        print(f"\n⏳ Waiting 5 seconds for environment to be ready...")
        time.sleep(5)
    
    # Create pathogen spreadsheet
    create_pathogen_spreadsheet(PATHOGENS, args.outdir, args.output_csv)
    
    if not args.skip_download:
        print(f"\n📥 Starting genome downloads...")
        downloaded_count = 0
        total_count = len(PATHOGENS)
        
        for i, (taxid, kraken_name, hops_name, dataset_type) in enumerate(PATHOGENS, 1):
            print(f"\n[{i}/{total_count}] Processing {kraken_name}...")
            
            # Download genome
            result = download_genome(taxid, kraken_name, hops_name, dataset_type, args.outdir, args.env_name)
            if result:
                downloaded_count += 1
            
            # Small delay to avoid overwhelming NCBI servers
            time.sleep(2)
        
        print(f"\n✅ Download Summary:")
        print(f"   Downloaded: {downloaded_count}/{total_count} genomes")
        print(f"   Success rate: {(downloaded_count/total_count)*100:.1f}%")
    
    print(f"\n🎉 Setup Complete!")
    print(f"   Pathogen spreadsheet: {args.output_csv}")
    print(f"   Genome directory: {args.outdir}")
    print(f"\n📋 Next Steps:")
    print(f"   1. Copy {args.output_csv} to config/Pathogen_spreadsheet.csv")
    print(f"   2. Update config/config.yaml with your host references")
    print(f"   3. Run PIGSTI pipeline:")
    print(f"      snakemake --cores 50 --use-conda --rerun-incomplete")
    
    if not args.skip_env:
        print(f"\n🔧 Environment Management:")
        print(f"   Activate: conda activate {args.env_name}")
        print(f"   Deactivate: conda deactivate")

if __name__ == "__main__":
    main()
