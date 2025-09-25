import unittest
from app import app, db
from models.molecule import Molecule
from services.molecule_service import MoleculeService

class FlaskTestCase(unittest.TestCase):
    def setUp(self):
        self.app = app.test_client()
        self.app.testing = True
        with app.app_context():
            db.create_all()
    
    def tearDown(self):
        with app.app_context():
            db.session.remove()
            db.drop_all()
    
    def test_home_page(self):
        response = self.app.get('/')
        self.assertEqual(response.status_code, 200)
        self.assertIn(b'Chemical Structure Visualizer', response.data)
    
    def test_valid_smiles(self):
        response = self.app.post('/visualize', data=dict(
            smiles='CCO',
            name='Ethanol'
        ), follow_redirects=True)
        self.assertEqual(response.status_code, 200)
        self.assertIn(b'Ethanol', response.data)
    
    def test_invalid_smiles(self):
        response = self.app.post('/visualize', data=dict(
            smiles='INVALID',
            name='Invalid'
        ), follow_redirects=True)
        self.assertEqual(response.status_code, 200)
        self.assertIn(b'Invalid SMILES string', response.data)
    
    def test_molecule_service(self):
        service = MoleculeService(db)
        molecule, error = service.create_molecule_from_smiles('CCO', 'Ethanol')
        self.assertIsNotNone(molecule)
        self.assertIsNone(error)
        self.assertEqual(molecule.name, 'Ethanol')
        self.assertEqual(molecule.formula, 'C2H6O')

if __name__ == '__main__':
    unittest.main()