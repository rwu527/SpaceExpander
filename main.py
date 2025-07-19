# main.py
import os
import shutil
import logging
import traceback
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors
import pandas as pd
from src.concrete_content import (
    generate_fragments,
    process_fragments,
    nest_fragments
)
from src.claims_generate import create_claims
from src.digital import convert_numbers
from src import global_counters



# Configure logging
logging.basicConfig(
    filename='application.log',
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s: %(message)s'
)


def convert_sdf_to_csv(sdf_file, output_folder):
    # Initialize logger for this function
    logger = logging.getLogger('sdf_converter')
    logger.setLevel(logging.INFO)
    
    # Create handlers
    console_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(os.path.join(output_folder, 'conversion.log'))
    
    # Create formatters and add to handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    
    logger.info(f"Starting conversion of {sdf_file}")
    
    try:
        # Initialize standardization tools
        normalizer = rdMolStandardize.Normalizer()
        fragment_chooser = rdMolStandardize.LargestFragmentChooser()
        uncharger = rdMolStandardize.Uncharger()
        metal_disconnector = rdMolStandardize.MetalDisconnector()
    except Exception as e:
        logger.error(f"Failed to initialize standardization tools: {str(e)}")
        return None
    
    # Read SDF file
    try:
        suppl = Chem.SDMolSupplier(sdf_file)
        mol_count = len(suppl)
        logger.info(f"Found {mol_count} molecules in SDF file")
    except Exception as e:
        logger.error(f"Failed to read SDF file: {str(e)}")
        return None
    
    if mol_count == 0:
        logger.error("SDF file contains no molecules")
        return None
    
    data = []
    processed_count = 0
    error_count = 0
    
    for idx, mol in enumerate(suppl):
        if mol is None:
            logger.warning(f"Skipping invalid molecule at position {idx+1}")
            error_count += 1
            continue

        try:
            # Save original SMILES for comparison
            original_smiles = Chem.MolToSmiles(mol) if mol.GetNumAtoms() > 0 else ""
            
            # Standardization steps
            mol = metal_disconnector.Disconnect(mol)
            mol = normalizer.normalize(mol)
            mol = fragment_chooser.choose(mol)
            mol = uncharger.uncharge(mol)
            
            # Skip empty molecules
            if mol is None or mol.GetNumAtoms() == 0:
                logger.warning(f"Molecule at position {idx+1} became empty after standardization")
                error_count += 1
                continue
                
            # Sanitize and generate SMILES
            Chem.SanitizeMol(mol)
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=False)
            
            # Add molecular weight for debugging
            mol_weight = Descriptors.MolWt(mol)
            
            data.append([smiles, mol_weight])
            processed_count += 1
            
            # Log transformation if significant change
            if original_smiles and smiles != original_smiles:
                logger.info(f"Molecule {idx+1} transformed: {original_smiles} -> {smiles}")
            
        except Exception as e:
            logger.error(f"Failed to process molecule at position {idx+1}: {str(e)}")
            error_count += 1
            continue

    logger.info(f"Processed {processed_count} molecules successfully, {error_count} errors")
    
    if processed_count > 0:
        try:
            csv_file = os.path.join(output_folder, os.path.basename(sdf_file).replace(".sdf", ".csv"))
            df = pd.DataFrame(data, columns=["SMILES", "MolWeight"])
            df.to_csv(csv_file, index=False)
            logger.info(f"Saved {processed_count} SMILES to {csv_file}")
            return csv_file
        except Exception as e:
            logger.error(f"Failed to write CSV file: {str(e)}")
            return None
    else:
        logger.error("Conversion failed: No valid molecules processed")
        return None
    
    

def remove_empty_folders_upwards(folder_path):
    for entry in os.listdir(folder_path):
        entry_path = os.path.join(folder_path, entry)
        if os.path.isdir(entry_path):
            remove_empty_folders_upwards(entry_path)
    if not os.listdir(folder_path) or all(f.endswith('.json') for f in os.listdir(folder_path)):
        shutil.rmtree(folder_path)


def run_pipeline(file_path, output_folder):
    global_counters.global_r_group_counter = 1
    global_counters.global_x_group_counter = 1
    global_counters.global_z_group_counter = 1

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder)

    if file_path.endswith('.sdf'):
        csv_file_path = convert_sdf_to_csv(file_path, output_folder)
        if not csv_file_path:
            raise ValueError("Failed to convert SDF to CSV")
        file_path = csv_file_path
    elif not file_path.endswith('.csv'):
        raise ValueError("Only SDF or CSV files are supported")

    total_molecules, r_group_to_atom_info = generate_fragments(file_path, output_folder)
    process_fragments(output_folder)
    nest_fragments(output_folder, r_group_to_atom_info, total_molecules)
    convert_numbers(output_folder, mode="unicode")
    remove_empty_folders_upwards(output_folder)
    create_claims(output_folder)

    return output_folder


if __name__ == "__main__":
    try:
        input_path = input("请输入输入文件路径（SDF 或 CSV）: ").strip()
        output_path = os.path.join(os.path.dirname(input_path), "output")

        result_folder = run_pipeline(input_path, output_path)

    except Exception as e:
        error_message = traceback.format_exc()
        logging.error(f"An error occurred:\n{error_message}")
        print(f"❌ 错误发生：{e}")



