from equilibrator_api import ComponentContribution, Q_
import equilibrator_cache
import logging


class rpEquilibrator:
    def __init__(self, rpsbml=None, ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15):
        self.logger = logging.getLogger(__name__)
        self.logger.debug('Started instance of rpEquilibrator')
        self.cc = ComponentContribution()
        self.cc.p_h = Q_(ph)
        self.cc.ionic_strength = Q_(str(ionic_strength)+' mM')
        self.cc.p_mg = Q_(pMg)
        self.cc.temperature = Q_(str(temp_k)+' K')
        self.rpsbml = rpsbml
    

    def _makeSpeciesStr(self, libsbml_species):
        """
        example: {'inchikey': ['GPRLSGONYQIRFK-UHFFFAOYSA-N'], 'seed': ['cpd00067'], 'sabiork': ['39'], 'reactome': ['R-ALL-74722', 'R-ALL-70106', 'R-ALL-5668577', 'R-ALL-428548', 'R-ALL-428040', 'R-ALL-427899', 'R-ALL-425999', 'R-ALL-425978', 'R-ALL-425969', 'R-ALL-374900', 'R-ALL-372511', 'R-ALL-351626', 'R-ALL-2872447', 'R-ALL-2000349', 'R-ALL-194688', 'R-ALL-193465', 'R-ALL-163953', 'R-ALL-156540', 'R-ALL-1470067', 'R-ALL-113529', 'R-ALL-1132304'], 'metacyc': ['PROTON'], 'hmdb': ['HMDB59597'], 'chebi': ['5584', '13357', '10744', '15378'], 'bigg': ['M_h', 'h'], 'metanetx': ['MNXM89553', 'MNXM145872', 'MNXM1', 'MNXM01']}
        Take a libsbml species object, parse the MIRIAM or the brsynth (if present) to return 
        the equilibrator appropriate string. The order of preference is the following:
        -KEGG
        -CHEBI
        #-bigg
        #-MNX
        -inchikey
        
        TODO: metanetx.chemical:MNXM7 + bigg.metabolite:pi
        """
        annot = libsbml_species.getAnnotation()
        if not annot:
            self.logger.error('Cannot retreive the annotation')
            return False
        miriam_dict = self.rpsbml.readMIRIAMAnnotation(annot)
        if not miriam_dict:
            self.logger.error('The object annotation does not have any MIRIAM entries')
            return False
        if 'kegg' in miriam_dict:
            if miriam_dict['kegg']:
                try:
                    #take the lowest value
                    int_list = [int(i.replace('C', '')) for i in miriam_dict['kegg']]
                    return 'KEGG:'+str(miriam_dict['kegg'][int_list.index(min(int_list))])
                except ValueError:
                    self.logger.warning('There is a non int value in: '+str(miriam_dict['kegg']))
        elif 'chebi' in miriam_dict:
            if miriam_dict['chebi']:
                try:
                    #take the lowest value
                    int_list = [int(i) for i in miriam_dict['chebi']]
                    return 'CHEBI:'+str(miriam_dict['chebi'][int_list.index(min(int_list))])
                except ValueError:
                    self.logger.warning('There is a non int value in: '+str(miriam_dict['chebi']))
        elif 'inchikey' in miriam_dict:
            if miriam_dict['inchikey']:
                if len(miriam_dict['inchikey'])==1:
                    return miriam_dict['inchikey'][0]
                else:
                    self.logger.warning('There are multiple values of inchikey: '+str(miriam_dict['inchikey']))
                    self.logger.warning('Taking the first one')
                    return miriam_dict['inchikey'][0]
        else:
            self.logger.warning('Could not extract string input for '+str(miriam_dict))
            return False
        self.logger.warning('Got to the end without return... should not happen')
        return False


    def species(self, libsbml_species, write_results=False):
        """
        Return the formation energy of a chemical species
        """
        return False


    def reaction(self, libsbml_reaction, write_results=False):
        """
        Build the string reaction from a libSBML reaction object to send to equilibrator and return the reversibility index, the different 
        """
        reac_str = ''
        for rea in libsbml_reaction.getListOfReactants():
            rea_str = self._makeSpeciesStr(self.rpsbml.model.getSpecies(rea.getSpecies()))
            if rea_str:
                reac_str += str(rea.getStoichiometry())+' '+str(rea_str)+' + '
            else:
                return False
        reac_str = reac_str[:-2]
        reac_str += '<=> ' #TODO: need to find a way to determine the reversibility of the reaction
        for pro in libsbml_reaction.getListOfProducts():
            pro_str = self._makeSpeciesStr(self.rpsbml.model.getSpecies(pro.getSpecies()))
            if pro_str:
                reac_str += str(pro.getStoichiometry())+' '+str(pro_str)+' + '
            else:
                return False
        reac_str = reac_str[:-2]
        self.logger.debug('reac_str: '+str(reac_str))
        try:
            rxn = self.cc.parse_reaction_formula(reac_str)
            standard_dg = self.cc.standard_dg(rxn)
            standard_dg_prime = self.cc.standard_dg_prime(rxn)
            physiological_dg_prime = self.cc.physiological_dg_prime(rxn)
            ln_reversibility_index = self.cc.ln_reversibility_index(rxn)
            self.logger.debug(rxn.is_balanced())
            self.logger.debug('ln_reversibility_index: '+str(ln_reversibility_index.value.m))
            self.logger.debug('standard_dg.value.m: '+str(standard_dg.value.m))
            self.logger.debug('standard_dg.error.m: '+str(standard_dg.error.m))
            self.logger.debug('standard_dg_prime.value.m: '+str(standard_dg_prime.value.m))
            self.logger.debug('standard_dg_prime.error.m: '+str(standard_dg_prime.error.m))
            self.logger.debug('physiological_dg_prime.value.m: '+str(physiological_dg_prime.value.m))
            self.logger.debug('physiological_dg_prime.error.m: '+str(physiological_dg_prime.error.m))
            if write_results:
                self.rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_o', standard_dg_prime.value.m, 'kj_per_mol')
                self.rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_m', physiological_dg_prime.value.m, 'kj_per_mol')
                self.rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_uncert', standard_dg.error.m, 'kj_per_mol')
                self.rpsbml.addUpdateBRSynth(libsbml_reaction, 'reversibility_index', ln_reversibility_index.value.m)
                self.rpsbml.addUpdateBRSynth(libsbml_reaction, 'balanced', rxn.is_balanced())
            return (rxn.is_balanced(),
                   (float(ln_reversibility_index.value.m), float(ln_reversibility_index.error.m)),
                   (float(standard_dg.value.m), float(standard_dg.error.m)),
                   (float(standard_dg_prime.value.m), float(standard_dg_prime.error.m)), 
                   (float(physiological_dg_prime.value.m), float(physiological_dg_prime.error.m)))
        except equilibrator_cache.exceptions.ParseException:
            self.logger.warning('One of the reaction species cannot be parsed by equilibrator: '+str(reac_str))
            return False
        except equilibrator_cache.exceptions.MissingDissociationConstantsException:
            self.logger.warning('Some of the species have not been pre-caclulated using ChemAxon')
            return False

    def MDF():
        pass

#used to initialise and download the data for equilibrator
if __name__ == "__main__":
    from equilibrator_api import ComponentContribution, Q_
    cc = ComponentContribution()


#######################
#m = component_contribution.molecule.Molecule()
#s = m.FromInChI('InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+')
