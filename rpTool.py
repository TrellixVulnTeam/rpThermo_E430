import rpEquilibrator
import rpComponentContribution

class rpThermo:
    def __init__(self, rpsbml=None, kegg_dG={}, cc_preprocess={}, ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15):
        self.rpsbml = rpsbml
        self.ph = ph
        self.ionic_strength = ionic_strength
        self.pMg = pMg
        self.temp_k = temp_k
        self.rpequilibrator = rpEquilibrator.rpEquilibrator(self.rpsbml, self.ph, self.ionic_strength, self.pMg, self.temp_k)
        self.rpcomponentcontribution = rpComponentContribution.rpComponentContribution(self.rpsbml, self.ph, self.pMg, self.ionic_strength/100.0, self.temp_k)
        self.rpcomponentcontribution.kegg_dG = kegg_dG 
        self.rpcomponentcontribution.cc_preprocess = cc_preprocess

    def rpSBMLPass(rpsbml):
        self.rpsbml = rpsbml
        self.rpequilibrator.rpsbml = rpsbml
        self.rpcomponentcontribution.rpsbml = rpsbml


    def pathway(self, pathway_id='rp_pathway', write_results=True):
        """
        WARNING: taking the sum of the reaction thermodynamics is perhaps not the best way to do it
        Using the equilibrator-api as well as the legacy component contribution method to calculate the mean thermodynamics of the rp_pathway
        """
        groups = self.rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        if not rp_pathway:
            self.logger.error('Cannot retreive the pathway: '+str(pathway_id))
            reurn False
        pathway_balanced = []
        pathway_reversibility_index = []
        pathway_reversibility_index_error = []
        pathway_standard_dg = []
        pathway_standard_dg_error = []
        pathway_standard_dg_prime = []
        pathway_standard_dg_prime_error = []
        pathway_physiological_dg_prime = []
        pathway_physiological_dg_prime_error = []
        for react in [self.rpsbml.model.getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]:
            res = self.rpequilibrator.reaction(react, write_results)
            if res:
                #WARNING: the uncertainty for the three thermo calculations should be the same
                pathway_balanced.append(res[0])
                pathway_reversibility_index.append(res[1][0])
                pathway_reversibility_index_error.append(res[1][1])
                #ignoring --  need to see if legacy component contribution can return these values
                #pathway_standard_dg.append(res[2][0])
                #pathway_standard_dg_error.append(res[2][1])
                pathway_standard_dg_prime.append(res[3][0])
                pathway_standard_dg_prime_error.append(res[3][1])
                pathway_physiological_dg_prime.append(res[4][0])
                pathway_physiological_dg_prime_error.append(res[4][1])
            else:
                self.logger.debug('Reverting to legacy component contribution')
                standard_dg_prime, physiological_dg_prime, dg_prime_uncert = self.rpcomponentcontribution.reaction(react, write_results)
                res = self.rpcomponentcontribution.reaction(react, write_results)
                if res:
                    pathway_standard_dg_prime.append(res[0])
                    pathway_standard_dg_prime_error.append(res[2])
                    pathway_physiological_dg_prime.append(res[1])
                    pathway_physiological_dg_prime_error.append(res[2])
                    #TODO: need to implement
                    pathway_balanced.append(None)
                    pathway_reversibility_index.append(None)
                    pathway_reversibility_index_error.append(None)
                else:
                    self.logger.error('Cannot calculate the thermodynmics for the reaction: '+str(react))
                    return False
        #WARNING return is ignoring balanced and reversibility index -- need to implement in legacy to return it (however still writing these results to the SBML)
        if write_results:
            self.rpsbml.addUpdateBRSynth(reac, 'dfG_prime_o', np.mean(pathway_standard_dg_prime), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(reac, 'dfG_prime_o_std', np.std(pathway_standard_dg_prime), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(reac, 'dfG_prime_m', np.mean(pathway_physiological_dg_prime), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(reac, 'dfG_prime_m_std', np.std(pathway_physiological_dg_prime), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(reac, 'dfG_uncert', np.mean(pathway_standard_dg_prime_error), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(reac, 'dfG_uncert_std', np.std(pathway_standard_dg_prime_error), 'kj_per_mol')
        return (np.median(pathway_standard_dg_prime), np.std(pathway_standard_dg_prime)), (np.median(pathway_physiological_dg_prime), np.std(pathway_physiological_dg_prime)), (np.median(pathway_standard_dg_prime), np.std(pathway_standard_dg_prime))

