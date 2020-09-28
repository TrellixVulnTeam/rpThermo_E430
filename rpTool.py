import rpEquilibrator
import logging
import numpy as np
import json


class rpThermo:
    """
    TODO: drop the formation energy of chemical species in the SBML file
    """
    def __init__(self, kegg_dG={}, cc_preprocess={}, ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpThermo')
        self.rpsbml = None
        self.ph = ph
        self.ionic_strength = ionic_strength
        self.pMg = pMg
        self.temp_k = temp_k
        self.rpequilibrator = rpEquilibrator.rpEquilibrator(self.rpsbml, self.ph, self.ionic_strength, self.pMg, self.temp_k)

    ## Small contructor function that passes the rpsbml object to the other classes
    #
    #
    def rpSBMLPass(self, rpsbml):
        self.logger.debug('Passing rpsbml: '+str(rpsbml))
        self.rpsbml = rpsbml
        self.rpequilibrator.rpsbml = rpsbml


    def pathway(self, pathway_id='rp_pathway', write_results=True):
        """
        WARNING: taking the sum of the reaction thermodynamics is perhaps not the best way to do it
        Using the equilibrator-api as well as the legacy component contribution method to calculate the mean thermodynamics of the rp_pathway
        """
        groups = self.rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        if not rp_pathway:
            self.logger.error('Cannot retreive the pathway: '+str(pathway_id))
            return False
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
            self.logger.debug('Sending the following reaction to reactionStrQuery: '+str(react))
            res = self.rpequilibrator.reactionStrQuery(react, write_results)
            self.logger.debug('The result is :'+str(res))
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
                self.logger.info('Native equilibrator string query failed')
                self.logger.info('Trying equilibrator_api component contribution')
                self.logger.debug('Trying to calculate using CC: '+str(react))
                res = self.rpequilibrator.reactionCmpQuery(react, write_results)
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
                    if write_results:
                        self.logger.warning('Adding 0 values everything to 0')
                        pathway_standard_dg_prime.append(0.0)
                        pathway_standard_dg_prime_error.append(0.0)
                        pathway_physiological_dg_prime.append(0.0)
                        pathway_physiological_dg_prime_error.append(0.0)
                        #TODO: need to implement
                        pathway_balanced.append(None)
                        pathway_reversibility_index.append(None)
                        pathway_reversibility_index_error.append(None)
                    return False
        #WARNING return is ignoring balanced and reversibility index -- need to implement in legacy to return it (however still writing these results to the SBML)
        if write_results:
            self.rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_prime_o', np.sum(pathway_standard_dg_prime), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_prime_o_std', np.std(pathway_standard_dg_prime), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_prime_m', np.sum(pathway_physiological_dg_prime), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_prime_m_std', np.std(pathway_physiological_dg_prime), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_uncert', np.mean(pathway_standard_dg_prime_error), 'kj_per_mol')
            self.rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_uncert_std', np.std(pathway_standard_dg_prime_error), 'kj_per_mol')
        return (np.sum(pathway_standard_dg_prime), np.std(pathway_standard_dg_prime)), (np.sum(pathway_physiological_dg_prime), np.std(pathway_physiological_dg_prime)), (np.sum(pathway_standard_dg_prime), np.std(pathway_standard_dg_prime))

