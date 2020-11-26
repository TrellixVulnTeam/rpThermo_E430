import unittest
import os
import sys

sys.path.insert(0, '..')

import rpEquilibrator
#WARNING: Need to copy a version of rpSBML locally
import rpSBML

class TestRPequilibrator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.rpsbml = rpSBML.rpSBML('test', path=os.path.join('data', 'rpsbml.xml'))
        self.rpeq = rpEquilibrator.rpEquilibrator(self.rpsbml)

    def test_pathway(self):
        res = self.rpeq.pathway()
        self.assertAlmostEqual(res[0][0], -148.48500016184292)
        self.assertAlmostEqual(res[0][1], 0.0)
        self.assertAlmostEqual(res[1][0], -165.59983769656395)
        self.assertAlmostEqual(res[1][1], 0,0)
        self.assertAlmostEqual(res[2][0], -148.48500016184292)
        self.assertAlmostEqual(res[2][1], 0.0)
        #check that everything has been written to the file
        res_json = self.rpsbml.genJSON()
        self.assertAlmostEqual(res_json['pathway']['brsynth']['dfG_prime_m']['value'], -165.59983769656395)

    def test_makeSpeciesStr(self):
        spe = self.rpsbml.model.getSpecies('MNXM100__64__MNXC3')
        self.assertIsNotNone(spe)
        self.assertTrue(self.rpeq._makeSpeciesStr(spe)=='CHEBI:5332')

    def test_makeReactionStr(self):
        reac = self.rpsbml.model.getReaction('RP1')
        self.assertIsNotNone(reac)
        self.assertTrue(self.rpeq._makeReactionStr(reac)=='1.0 CHEBI:5332 <=> 1.0 XMGQYMWWDOXHJM-UHFFFAOYSA-N + 1.0 CHEBI:8683 ')

    def test_speciesCmpQuery(self):
        reac = self.rpsbml.model.getSpecies('TARGET_0000000001__64__MNXC3')
        self.assertIsNotNone(reac)
        mu, sigma = self.rpeq._speciesCmpQuery(reac)
        self.assertAlmostEqual(mu, 143.49986576845282)
        
    def test_MDF(self):
        rpsbml = rpSBML.rpSBML('test', path=os.path.join('data', 'rpsbml_pathway.xml'))
        rpeq = rpEquilibrator.rpEquilibrator(rpsbml)
        mdf = rpeq.MDF()
        self.assertAlmostEqual(mdf, 205.53445861091302)
