from equilibrator_api import ComponentContribution, Q_
import equilibrator_cache
import logging
import json
import tempfile


class rpEquilibrator:
    """
    Collection of functions to intereact between rpSBML files and equilibrator. Includes a function to convert an rpSBML file to a SBtab format for MDF analysis
    """
    def __init__(self, rpsbml=None, ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
        self.logger = logging.getLogger(__name__)
        self.logger.debug('Started instance of rpEquilibrator')
        self.cc = ComponentContribution()
        self.cc.p_h = Q_(ph)
        self.cc.ionic_strength = Q_(str(ionic_strength)+' mM')
        self.cc.p_mg = Q_(pMg)
        self.cc.temperature = Q_(str(temp_k)+' K')
        self.stdev_factor = stdev_factor
        self.ph = ph
        self.ionic_strength = ionic_strength
        self.pMg = pMg
        self.temp_k = temp_k
        self.mnx_default_conc = json.load(open('data/mnx_default_conc.json', 'r'))
        self.rpsbml = rpsbml
    

    ##################################################################################
    ############################### PRIVATE ##########################################
    ##################################################################################


    def _makeSpeciesStr(self, libsbml_species, ret_type='xref'):
        """
        example: {'inchikey': ['GPRLSGONYQIRFK-UHFFFAOYSA-N'], 'seed': ['cpd00067'], 'sabiork': ['39'], 'reactome': ['R-ALL-74722', 'R-ALL-70106', 'R-ALL-5668577', 'R-ALL-428548', 'R-ALL-428040', 'R-ALL-427899', 'R-ALL-425999', 'R-ALL-425978', 'R-ALL-425969', 'R-ALL-374900', 'R-ALL-372511', 'R-ALL-351626', 'R-ALL-2872447', 'R-ALL-2000349', 'R-ALL-194688', 'R-ALL-193465', 'R-ALL-163953', 'R-ALL-156540', 'R-ALL-1470067', 'R-ALL-113529', 'R-ALL-1132304'], 'metacyc': ['PROTON'], 'hmdb': ['HMDB59597'], 'chebi': ['5584', '13357', '10744', '15378'], 'bigg': ['M_h', 'h'], 'metanetx': ['MNXM89553', 'MNXM145872', 'MNXM1', 'MNXM01']}
        Take a libsbml species object, parse the MIRIAM or the brsynth (if present) to return 
        the equilibrator appropriate string. The order of preference is the following:
        -KEGG
        -CHEBI
        #-bigg
        #-MNX
        -inchikey
        ret_type -> valid options (xref, id, name)
        TODO: metanetx.chemical:MNXM7 + bigg.metabolite:pi
        """
        if ret_type=='name':
            return libsbml_species.getName()
        elif ret_type=='id':
            return libsbml_species.getId()
        elif ret_type=='xref':
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
        else:
            self.logger.warning('Cannot determine ret_type: '+str(ret_type))


    def _makeReactionStr(self, libsbml_reaction, ret_type='xref', ret_stoichio=True):
        reac_str = ''
        for rea in libsbml_reaction.getListOfReactants():
            rea_str = self._makeSpeciesStr(self.rpsbml.model.getSpecies(rea.getSpecies()), ret_type)
            if rea_str:
                if ret_stoichio:
                    reac_str += str(rea.getStoichiometry())+' '+str(rea_str)+' + '
                else:
                    reac_str += str(rea_str)+' + '
            else:
                return False
        reac_str = reac_str[:-2]
        reac_str += '<=> ' #TODO: need to find a way to determine the reversibility of the reaction
        for pro in libsbml_reaction.getListOfProducts():
            pro_str = self._makeSpeciesStr(self.rpsbml.model.getSpecies(pro.getSpecies()), ret_type)
            if pro_str:
                if ret_stoichio:
                    reac_str += str(pro.getStoichiometry())+' '+str(pro_str)+' + '
                else:
                    reac_str += str(pro_str)+' + '
            else:
                return False
        reac_str = reac_str[:-2]
        self.logger.debug('reac_str: '+str(reac_str))
        return reac_str


    ################################################################################
    ########################### PUBLIC FUNCTIONS ###################################
    ################################################################################


    def species(self, libsbml_species, write_results=False):
        """
        Return the formation energy of a chemical species
        """
        return False


    def reactionThermo(self, libsbml_reaction, write_results=False):
        """
        Build the string reaction from a libSBML reaction object to send to equilibrator and return the different thermodynamics analysis available
        """
        #TODO: when an inchikey is passed, (and you don't have any other xref) and equilibrator finds the correct species then update the MIRIAM annotations
        try:
            rxn = self.cc.parse_reaction_formula(self._makeReactionStr(libsbml_reaction))
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


    def toNetworkSBtab(self, output, pathway_id='rp_pathway'):
        """
        Convert an SBML pathway to a simple network for input to equilibrator-pathway for MDF
        """
        groups = self.rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        all_reac_ids = []
        if not rp_pathway:
            self.logger.error('Cannot retreive the pathway: '+str(pathway_id))
            return False
        with open(output, 'w') as fo:
            fo.write("!!!SBtab DocumentName='E. coli central carbon metabolism - balanced parameters' SBtabVersion='1.0'\t\t\t\n")
            fo.write("!!SBtab TableID='Configuration' TableType='Config'\t\t\t\n")
            fo.write("!Option\t!Value\t!Comment\t\n")
            fo.write("algorithm\tMDF\tECM, or MDF\t\n")
            fo.write("p_h\t"+str(self.ph)+"\t\t\n")
            fo.write("ionic_strength\t"+str(self.ionic_strength)+" mM\t\t\n")
            fo.write("p_mg\t"+str(self.pMg)+"\t\t\n")
            fo.write("stdev_factor    "+str(self.stdev_factor)+"\n")
            fo.write("\t\t\t\n")
            fo.write("!!SBtab TableID='Reaction' TableType='Reaction'\t\t\t\n")
            fo.write("!ID\t!ReactionFormula\t\t\n")
            for react in [self.rpsbml.model.getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]:
                all_reac_ids.append(react.getId())
                react_str = self._makeReactionStr(react, 'id', True)
                if react_str:
                    fo.write(str(react.getId())+"\t"+str(react_str)+"\n")
                else:
                    self.logger.error('Cannot build the reaction: '+str(rect))
                    return False
            fo.write("\t\t\t\n")
            fo.write("\t\t\t\n")
            fo.write("!!SBtab TableID='Compound' TableType='Compound'\t\t\t\n")
            fo.write("!ID\t!Identifiers\t\t\n")
            rp_species = self.rpsbml.readUniqueRPspecies(pathway_id)
            for spe_id in rp_species:
                spe = self.rpsbml.model.getSpecies(spe_id)
                miriam_dict = self.rpsbml.readMIRIAMAnnotation(spe.getAnnotation())
                if not miriam_dict:
                    self.logger.error('The object annotation does not have any MIRIAM entries')
                    return False
                iden_str = None
                if 'kegg' in miriam_dict:
                    if miriam_dict['kegg']:
                        try:
                            #take the lowest value
                            int_list = [int(i.replace('C', '')) for i in miriam_dict['kegg']]
                            iden_str = 'KEGG:'+str(miriam_dict['kegg'][int_list.index(min(int_list))])
                        except ValueError:
                            self.logger.warning('There is a non int value in: '+str(miriam_dict['kegg']))
                if 'chebi' in miriam_dict and not iden_str:
                    if miriam_dict['chebi']:
                        try:
                            #take the lowest value
                            int_list = [int(i) for i in miriam_dict['chebi']]
                            iden_str = 'CHEBI:'+str(miriam_dict['chebi'][int_list.index(min(int_list))])
                        except ValueError:
                            self.logger.warning('There is a non int value in: '+str(miriam_dict['chebi']))
                if 'inchikey' in miriam_dict and not iden_str:
                    if miriam_dict['inchikey']:
                        if len(miriam_dict['inchikey'])==1:
                            iden_str = miriam_dict['inchikey'][0]
                        else:
                            self.logger.warning('There are multiple values of inchikey: '+str(miriam_dict['inchikey']))
                            self.logger.warning('Taking the first one')
                            iden_str = miriam_dict['inchikey'][0]
                if not iden_str:
                    self.logger.warning('Could not extract string input for '+str(miriam_dict))
                    return False
                fo.write(str(spe_id)+"\t"+str(iden_str)+"\t\t\n")
            fo.write("\t\t\t\n")
            #TODO: perhaps find a better way than just setting this to 1
            fo.write("!!SBtab TableID='Flux' TableType='Quantity' Unit='mM/s'\t\t\t\n")
            fo.write("!QuantityType\t!Reaction\t!Value\t\n")
            for rea_id in all_reac_ids:
                fo.write("rate of reaction\t"+str(rea_id)+"\t1\t\n")
            fo.write("\t\t\t\n")
            fo.write("!!SBtab TableID='ConcentrationConstraint' TableType='Quantity' Unit='mM'\t\t\t\n")
            fo.write("!QuantityType\t!Compound\t!Min\t!Max\n")
            for spe_id in rp_species: 
                self.logger.debug('========= '+str(spe_id)+' ========')
                is_found = False
                spe = self.rpsbml.model.getSpecies(spe_id)
                miriam_dict = self.rpsbml.readMIRIAMAnnotation(spe.getAnnotation())
                self.logger.debug(miriam_dict)
                if not miriam_dict:
                    self.logger.warning('The object annotation does not have any MIRIAM entries')
                    continue
                if 'metanetx' in miriam_dict:
                    self.logger.debug(miriam_dict['metanetx'])
                    for mnx in miriam_dict['metanetx']:
                        if mnx in list(self.mnx_default_conc.keys()) and not is_found:
                            self.logger.debug('Found default concentration range for '+str(spe.getId())+' ('+str(mnx)+'): '+str(self.mnx_default_conc[mnx]))
                            fo.write("concentration\t"+spe.getId()+"\t"+str(self.mnx_default_conc[mnx]['c_min'])+"\t"+str(self.mnx_default_conc[mnx]['c_max'])+"\n")
                            is_found = True
                if not is_found:
                    self.logger.debug('Using default values for '+str(spe.getId()))
                    fo.write("concentration\t"+spe.getId()+"\t0.001\t10\n")
        return True


    def MDF(self, pathway_id='rp_pathway', write_results=True):
        """
        Perform MDF analysis on the retropath pathways
        """
        from equilibrator_pathway import Pathway
        to_ret_mdf = None
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            path_sbtab = os.path.join(tmpOutputFolder, 'tmp_sbtab.tsv')
            self.toNetworkSBtab(path_sbtab, pathway_id)
            pp = Pathway.from_sbtab(path_sbtab, comp_contrib=self.cc)
            pp.update_standard_dgs()
            mdf_sol = pp.calc_mdf()
            to_ret_mdf = float(mdf_sol.mdf)
            self.rpsbml.addUpdateBRSynth(rp_pathway, 'MDF', float(mdf_sol.mdf), 'kj_per_mol')
        return to_ret_mdf


#used to initialise and download the data for equilibrator
if __name__ == "__main__":
    from equilibrator_api import ComponentContribution, Q_
    cc = ComponentContribution()


#######################
#m = component_contribution.molecule.Molecule()
#s = m.FromInChI('InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+')
