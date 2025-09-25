# app.py
from flask import Flask, render_template, request, jsonify, redirect, url_for, flash, Response
from extensions import db  # Import db from extensions
from config import Config
from services.molecule_service import MoleculeService
from services.visualization_service import VisualizationService
from rdkit import Chem
import os

app = Flask(__name__)
app.config.from_object(Config)

# Initialize the db with the app
db.init_app(app)

# Initialize services
molecule_service = MoleculeService()
visualization_service = VisualizationService()

# ===== ROUTES =====

@app.route('/')
def index():
    """Home page with molecule input form"""
    return render_template('index.html')

@app.route('/visualize', methods=['POST'])
def visualize():
    """Process SMILES input and create molecule"""
    smiles = request.form.get('smiles')
    name = request.form.get('name', 'Unnamed Molecule')
    
    if not smiles:
        flash('SMILES string is required', 'error')
        return redirect(url_for('index'))
    
    molecule, error = molecule_service.create_molecule_from_smiles(smiles, name)
    if error:
        flash(error, 'error')
        return redirect(url_for('index'))
    
    return redirect(url_for('visualize_molecule', molecule_id=molecule.id))

@app.route('/molecule/<int:molecule_id>')
def visualize_molecule(molecule_id):
    """Display a specific molecule in 3D"""
    molecule = molecule_service.get_molecule_by_id(molecule_id)
    if not molecule:
        flash('Molecule not found', 'error')
        return redirect(url_for('index'))
    
    mol = molecule_service.load_molecule_file(molecule.file_path)
    bond_angles = visualization_service.calculate_bond_angles(mol)
    
    return render_template('visualize.html', 
                          molecule=molecule, 
                          bond_angles=bond_angles,
                          atom_colors=visualization_service.get_atom_colors(),
                          bond_colors=visualization_service.get_bond_colors())

@app.route('/gallery')
def gallery():
    """Display gallery of all molecules"""
    molecules = molecule_service.get_all_molecules()
    return render_template('gallery.html', molecules=molecules)

@app.route('/api/molecule/<int:molecule_id>/data')
def get_molecule_data(molecule_id):
    """API endpoint to get molecule data for 3Dmol.js"""
    try:
        molecule = molecule_service.get_molecule_by_id(molecule_id)
        if not molecule:
            app.logger.error(f"Molecule with ID {molecule_id} not found")
            return jsonify({'error': 'Molecule not found'}), 404
        
        app.logger.info(f"Loading molecule {molecule_id} from {molecule.file_path}")
        
        # Check if file exists
        if not os.path.exists(molecule.file_path):
            app.logger.error(f"Molecule file not found: {molecule.file_path}")
            return jsonify({'error': 'Molecule file not found'}), 404
        
        mol = molecule_service.load_molecule_file(molecule.file_path)
        if mol is None:
            app.logger.error(f"Failed to load molecule from {molecule.file_path}")
            return jsonify({'error': 'Failed to load molecule'}), 500
        
        # Convert RDKit molecule to SDF format for 3Dmol.js
        sdf_data = Chem.MolToMolBlock(mol)
        
        app.logger.info(f"Successfully generated SDF data for molecule {molecule_id}")
        
        return jsonify({
            'sdf': sdf_data,
            'name': molecule.name,
            'formula': molecule.formula
        })
    except Exception as e:
        app.logger.error(f"Error processing molecule {molecule_id}: {str(e)}")
        return jsonify({'error': f'Internal server error: {str(e)}'}), 500
    

@app.route('/api/molecule/<int:molecule_id>/save', methods=['POST'])
def save_molecule(molecule_id):
    """Save molecule to user's collection (placeholder for future auth)"""
    # This would be implemented with user authentication
    return jsonify({'success': True, 'message': 'Molecule saved successfully'})

@app.route('/api/molecule/<int:molecule_id>/download')
def download_molecule(molecule_id):
    """Download molecule as SDF file"""
    molecule = molecule_service.get_molecule_by_id(molecule_id)
    if not molecule:
        return jsonify({'error': 'Molecule not found'}), 404
    
    mol = molecule_service.load_molecule_file(molecule.file_path)
    
    # Convert to SDF format
    sdf = Chem.MolToMolBlock(mol)
    
    return Response(
        sdf,
        mimetype='chemical/x-mdl-sdfile',
        headers={'Content-Disposition': f'attachment; filename={molecule.name}.sdf'}
    )

@app.route('/test-molecule')
def test_molecule():
    """Create a test molecule for debugging"""
    molecule, error = molecule_service.create_molecule_from_smiles('c1ccccc1', 'Benzene')
    if error:
        return f"Error: {error}"
    return redirect(url_for('visualize_molecule', molecule_id=molecule.id))
# ===== ERROR HANDLERS =====

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404

@app.errorhandler(500)
def internal_server_error(e):
    return render_template('500.html'), 500

# ===== CLI COMMANDS =====

@app.cli.command()
def init_db():
    """Initialize the database"""
    db.create_all()
    print('Database initialized')

@app.cli.command()
def clear_db():
    """Clear all data from the database"""
    if input('Are you sure? This will delete all data. (y/n): ') == 'y':
        db.drop_all()
        db.create_all()
        print('Database cleared')

# ===== MAIN EXECUTION =====

if __name__ == '__main__':
    with app.app_context():
        db.create_all()
    app.run(debug=True)