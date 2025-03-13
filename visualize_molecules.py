from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display, HTML
import pandas as pd
import base64
from io import BytesIO

def count_molecules(file_path="Results/New_Mols.txt"):
    """
    Counts the number of molecules in the SMILES file.
    
    Args:
        file_path (str): Path to the file containing SMILES strings
    
    Returns:
        int: Number of molecules in the file
    """
    with open(file_path, 'r') as file:
        return len(file.readlines())

def visualize_molecules(file_path="Results/New_Mols.txt", n_molecules=None):
    """
    Displays molecules from a SMILES file in a two-column table format.
    
    Args:
        file_path (str): Path to the file containing SMILES strings (one per line)
        n_molecules (int, optional): Number of molecules to display. If None, displays all molecules
    
    Returns:
        None: Displays the table in the Jupyter notebook
    """
    # Count total molecules
    total_molecules = count_molecules(file_path)
    
    # Read SMILES from file
    with open(file_path, 'r') as file:
        smiles_list = file.read().splitlines()
    
    # Limit number of molecules if specified
    if n_molecules is not None:
        smiles_list = smiles_list[:n_molecules]
    
    # Create molecules and images
    molecules_data = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                molecules_data.append({
                    'SMILES': smiles,
                    'Molecule': mol
                })
        except Exception as e:
            print(f"Error processing SMILES {smiles}: {str(e)}")
    
    # Create HTML content with centered text
    html_content = '<div style="margin-bottom: 20px; text-align: center;">'
    html_content += f'<p style="font-size: 18px; font-weight: bold;">Total number of molecules generated: {total_molecules}</p>'
    if n_molecules is not None:
        html_content += f'<p style="font-size: 16px; font-weight: bold;">Displaying first {n_molecules} molecules</p>'
    html_content += '</div>'
    
    # Create DataFrame
    df = pd.DataFrame(molecules_data)
    
    # Create HTML table with styled header and centered text
    html_content += '<table style="border-collapse: collapse; margin: 0 auto;">'
    html_content += '<tr>'
    html_content += '<th style="padding: 10px; border: 1px solid black; background-color: #2196F3; color: white; text-align: center;">Structure</th>'
    html_content += '<th style="padding: 10px; border: 1px solid black; background-color: #2196F3; color: white; text-align: center;">SMILES</th>'
    html_content += '</tr>'
    
    for _, row in df.iterrows():
        html_content += '<tr>'
        # Convert molecule to image and then to base64
        img = Draw.MolToImage(row['Molecule'])
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        
        html_content += f'<td style="padding: 10px; border: 1px solid black; text-align: center;">'
        html_content += f'<img src="data:image/png;base64,{img_str}">'
        html_content += '</td>'
        html_content += f'<td style="padding: 10px; border: 1px solid black; text-align: center;">{row["SMILES"]}</td>'
        html_content += '</tr>'
    
    html_content += '</table>'
    
    # Display the HTML content
    display(HTML(html_content)) 