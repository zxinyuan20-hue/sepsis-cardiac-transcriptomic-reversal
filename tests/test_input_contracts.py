"""Regression checks for real ingestion risks discovered during preparation."""
import gzip,importlib.util,sys,unittest
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
spec=importlib.util.spec_from_file_location("qc",ROOT/"scripts/02_validate_inputs.py")
qc=importlib.util.module_from_spec(spec);spec.loader.exec_module(qc)

class InputContracts(unittest.TestCase):
    def test_units_are_explicit_and_case_insensitive(self):
        out=qc.parse_dose_um(pd.Series(["10.0 um","0.37 uM","-666","10 mg/ml","10","1e-2 uM"]))
        self.assertEqual(out.iloc[0],10)
        self.assertEqual(out.iloc[1],.37)
        self.assertTrue(out.iloc[2:5].isna().all())
        self.assertEqual(out.iloc[5],.01)
    def test_fractional_counts_are_detected_not_rounded(self):
        x=pd.DataFrame({"a":[0.,1.2],"b":[2.,3.]},index=["g1","g2"])
        self.assertEqual(qc.matrix_profile(x)["noninteger_values"],1)
        self.assertEqual(x.loc["g2","a"],1.2)
    def test_duplicate_sample_content_is_reported(self):
        x=pd.DataFrame({"a":[1.,2.],"b":[1.,2.]})
        self.assertEqual(qc.matrix_profile(x)["exact_duplicate_sample_pairs"],[["b","a"]])
    def test_unlabelled_metrics_index_is_not_signature_id(self):
        folder=ROOT/"outputs/audit/test_fixtures";folder.mkdir(parents=True,exist_ok=True)
        p=folder/"metrics_extra_index.tsv.gz"
        with gzip.open(p,"wt") as f:f.write("sig_id\tdistil_cc_q75\n0\tplate:A01\t0.8\n")
        result=qc.read_metrics(p)
        self.assertEqual(result.sig_id.iloc[0],"plate:A01")
        self.assertEqual(float(result.distil_cc_q75.iloc[0]),.8)
    def test_sample_mapping_uses_explicit_replicate_not_order(self):
        meta=pd.DataFrame([{"gsm":"GSM2","title":"Heart_PBS_rep2"},{"gsm":"GSM1","title":"Heart_LPS_rep1"}])
        data=pd.DataFrame({"read_count_WT_LPS_1":[10,11],"read_count_WT_PBS_2":[20,21]},index=["a","b"])
        mp,out=qc.align_mouse(data,meta,"GSE267388")
        self.assertEqual(out.columns.tolist(),["GSM2","GSM1"])
        self.assertEqual(out.GSM2.iloc[0],20)
        self.assertEqual(mp.group.tolist(),["control","LPS"])
    def test_unknown_expression_column_blocks_mapping(self):
        meta=pd.DataFrame([{"gsm":"GSM2","title":"Heart_PBS_rep2"}])
        data=pd.DataFrame({"read_count_WT_PBS_2":[20],"read_count_WT_LPS_1":[10]})
        with self.assertRaises(ValueError):qc.align_mouse(data,meta,"GSE267388")

if __name__=="__main__":unittest.main(verbosity=2)
