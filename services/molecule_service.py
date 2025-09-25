# services/molecule_service.py
import os
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from models.molecule import Molecule
from extensions import db  # Import db from extensions

class MoleculeService:
    def __init__(self):
        self.upload_folder = 'data/molecules'
        
    def create_molecule_from_smiles(self, smiles, name):
        """Create a molecule object from SMILES string"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, "Invalid SMILES string"
            
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Calculate molecular formula
        formula = self._calculate_molecular_formula(mol)
        
        # Save to database
        molecule = Molecule(name=name, smiles=smiles, formula=formula)
        db.session.add(molecule)
        db.session.commit()
        
        # Save molecule file
        file_path = os.path.join(self.upload_folder, f"{molecule.id}.pkl")
        os.makedirs(self.upload_folder, exist_ok=True)
        with open(file_path, 'wb') as f:
            pickle.dump(mol, f)
            
        molecule.file_path = file_path
        db.session.commit()
        
        return molecule, None
        
    def get_molecule_by_id(self, molecule_id):
        """Get molecule by ID"""
        return Molecule.query.get(molecule_id)
        
    def get_all_molecules(self):
        """Get all molecules"""
        return Molecule.query.order_by(Molecule.created_at.desc()).all()
        
    def load_molecule_file(self, file_path):
        """Load molecule from pickle file"""
        with open(file_path, 'rb') as f:
            return pickle.load(f)
            
    def _calculate_molecular_formula(self, mol):
        """Calculate molecular formula from RDKit molecule"""
        formula = ""
        atoms = {}
        
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atoms[symbol] = atoms.get(symbol, 0) + 1
            
        # Sort atoms by symbol (C first, then H, then alphabetically)
        sorted_atoms = sorted(atoms.items(), key=lambda x: (0 if x[0] == 'C' else 1 if x[0] == 'H' else 2, x[0]))
        
        for symbol, count in sorted_atoms:
            formula += symbol
            if count > 1:
                formula += str(count)
                
        return formula
    
    
    
