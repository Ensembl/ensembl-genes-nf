import json
import unittest
from pathlib import Path


PIPELINE = Path(__file__).parents[1]


class DiamondScopeContractTests(unittest.TestCase):
    def test_scope_parameter_and_explicit_metadata_normalisation(self):
        main = (PIPELINE / "main.nf").read_text()
        config = (PIPELINE / "nextflow.config").read_text()
        schema = json.loads((PIPELINE / "nextflow_schema.json").read_text())

        self.assertIn("diamond_scope", main)
        self.assertIn("tuple(meta + [scope: 'cohort'], backend, bed)", main)
        self.assertIn("tuple(meta + [scope: 'accession'], backend, bed)", main)
        self.assertIn("cohort_diamond_inputs.mix(accession_diamond_inputs)", main)
        self.assertIn("diamond_scope = 'cohort'", config)
        self.assertEqual(schema["properties"]["diamond_scope"]["default"], "cohort")
        self.assertEqual(
            schema["properties"]["diamond_scope"]["enum"],
            ["cohort", "accession", "both"],
        )


if __name__ == "__main__":
    unittest.main()
