"""Packages pandas/RDKit. Inputs frozen public compound structures and PubChem records.
Outputs chemical standardization audit; retain stereo uncertainty and original identifiers.
No experimental material identity/purity inference. Does not alter statistical ranking.
"""
import json
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from pipeline_utils import ROOT,path,acquire,write_json,run_standard_module
def main():
 out=path('analysis')/'mechanism';raw=path('reference')/'mechanism'
 p=raw/'BRD-K65182930_name_properties.json';cid=json.loads(p.read_text(encoding='utf-8'))['PropertyTable']['Properties'][0]['CID']
 dest=raw/'NVP_AUY922_name_inchi.json';acquire(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChI,InChIKey,Title/JSON',dest,'PubChem name-derived structure for tautomer audit')
 dictionary=pd.read_csv(path('preparation')/'lincs_compound_dictionary.tsv',sep='\t',dtype=str).set_index('pert_id')
 table=pd.read_csv(out/'compound_identity_crosscheck.tsv',sep='\t');computed=[]
 for row in table.itertuples():
  mol=Chem.MolFromSmiles(dictionary.loc[row.pert_id,'canonical_smiles']);key=inchi.MolToInchiKey(mol)
  if key!=row.lincs_inchikey:raise ValueError('LINCS SMILES and key disagree')
  computed.append(key)
 table['smiles_computed_key']=computed
 original=dictionary.loc['BRD-K65182930'];pub=json.loads(dest.read_text(encoding='utf-8'))['PropertyTable']['Properties'][0]
 te=rdMolStandardize.TautomerEnumerator();a=Chem.MolFromSmiles(original.canonical_smiles);b=inchi.MolFromInchi(pub['InChI'])
 ka=inchi.MolToInchiKey(te.Canonicalize(a));kb=inchi.MolToInchiKey(te.Canonicalize(b))
 result={'lincs_key':original.inchi_key,'pubchem_name_key':pub['InChIKey'],'pubchem_name_cid':cid,'rdkit_lincs_canonical_tautomer_key':ka,'rdkit_pubchem_canonical_tautomer_key':kb,'tautomer_equivalent_under_rdkit':ka==kb,'limitation':'Computational standardization, not confirmation of experimental compound batch or purity'}
 write_json(out/'NVP_AUY922_tautomer_audit.json',result)
 table['reviewed_identity_status']=table.name_structure_status
 table.loc[table.pert_id=='BRD-K31342827','reviewed_identity_status']='source_structure_identifies_bisindolylmaleimide_I; generic_name_ambiguous'
 if ka==kb:table.loc[table.pert_id=='BRD-K65182930','reviewed_identity_status']='luminespib_tautomer_equivalent_under_RDKit'
 table.to_csv(out/'reviewed_compound_identity.tsv',sep='\t',index=False)
 return {'source_SMILES_key_consistency':len(table),'NVP_AUY922_tautomer_equivalent':ka==kb,'connectivity_only_stereo_unresolved':int((table.name_structure_status=='connectivity_only_agreement').sum())}
if __name__=='__main__':run_standard_module('26_chemical_standardization',main)
