from equilibrator_api import ComponentContribution, Q_

class rpEquilibrator:
    def __init__(self, ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15):
        self.cc = ComponentContribution()
        self.cc.p_h = Q_(ph)
        self.cc.ionic_strength = Q_(str(ionic_strength)+' mM')
        self.cc.p_mg = Q_(pMg)
        self.cc.temperature = Q_(str(temp_k)+' K')
    

    def _makeSpeciesStr():
        """
        Take a libsbml species object, parse the MIRIAM or the brsynth (if present) to return 
        the equilibrator appropriate string. The order of preference is the following:
        -KEGG
        -CHEBI
        -MNX
        -
        """

    def species(self, rpsbml_species, write_results=False):
        """
        TODO
        """
        pass

    def reaction(self, rpsbml_reaction, write_results=False):
        reac_str = ''
        for rea in rpsbml_reaction.getListOfReactants():
        for pro in rpsbml_reaction.getListOfProducts():

        #"kegg:C00002 + kegg:C00001 = kegg:C00008 + kegg:C00009"
        reac_string = ''
        rxn = self.cc.parse_reaction_formula(reac_string)
        ### need to test if all the species are recognised and if not you need to bail

    def pathway(self, rpsbml, pathway_id='rp_pathway', write_results=True):
        """
        WARNING: taking the sum of the reaction thermodynamics is perhaps not the best way to do it
        """
        groups = rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        for react in [rpsbml.model.getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]:
            self.reaction(react, write_results)

    def MDF():


#######################
#m = component_contribution.molecule.Molecule()
#s = m.FromInChI('InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)/b3-1+,4-2+')
