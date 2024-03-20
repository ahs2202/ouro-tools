from . import biobookshelf as bk

"""
Implementing Ontology functions
"""
import owlready2 as ol

class OntologyTerms :
    def __init__( self, path_file_owl : str, name_prefix : str, name_root_term : str ) :
        """ 
        load an ontology file, given as 'path_file_owl'
        
        name_prefix : str # prefix of the name
        name_root_term : str # name of the root term
        # 2024-02-29 22:50:55 
        """
        self.path_file_owl = path_file_owl
        self.name_prefix = name_prefix
        self.name_root_term = name_root_term
        self.onto = ol.get_ontology(f"file://{path_file_owl}").load( )
        self._set_terms = set( self.onto.classes( ) )
        self._root_term = self[ self.name_root_term ]
    def __repr__( self ) :
        return f"<{len( self._set_terms )} ontology terms stored at {self.path_file_owl}>"
    def __contains__( self, term ) :
        """
        # 2024-03-01 21:35:55 
        """
        return self[ term ] in self._set_terms
    def __getitem__( self, ontology_id : str ) :
        """
        get ontology term using an ID
        # 2024-02-29 22:51:49 
        """
        # handle 'ontology_id' that is not a string
        if not isinstance( ontology_id, str ) :
            ontology_id
        l = self.onto.search( iri = f"*{ontology_id}*")
        if len( l ) == 0 :
            return None 
        elif len( l ) == 1 :
            return l[ 0 ]
        else :
            return l # if more then one terms are matched, return more than one elements
    def __iter__( self ) :
        """
        return an iterater returning each class
        # 2024-03-01 00:54:16 
        """
        return self.onto.classes( )
    def get_ancestor_chain( self, term ) :
        """
        get a chain of ancestors (excluding restriction objects), excluding self, from the most distant ancestor (owl.Thing) to the closest ancestor.
        Note)
        This function utilizes a recursive algorithm to explore the tree structure.
        # 2024-03-01 21:45:16 
        """
        # get the ontology term
        term = self[ term ]
        # initialize the 'l_ancestor'
        def get_superclasses( term ) :
            """
            get filtered super classes of a term
            # 2024-03-01 23:15:04 
            """
            return list( e for e in term.is_a if hasattr( e, 'name' ) and e.name[ : len( self.name_prefix ) ] == self.name_prefix )
        def get_ancestor_chains( term ) :
            l_ancestor_chain = [ ]
            l_term_super = get_superclasses( term )
            # termination condition
            if len( l_term_super ) == 0 :
                if term == self._root_term :
                    return [ [ ] ]
                else : # if the chain terminate with a term that is not a root term, add the root term 
                    return [ [ self._root_term ] ]
            # recursive condition
            for e in l_term_super :
                for l_ancestor in get_ancestor_chains( e ) :
                    l_ancestor_chain.append( [ e ] + l_ancestor )
            return l_ancestor_chain
        # reverse the order (from the most distant ancestor (the root term) to the closest ancestor)
        l_ancestor_chain_reverse_order = get_ancestor_chains( term )
        l_ancestor_chain = [ ] # l_ancestor_chain
        for ancestor_chain_reverse_order in l_ancestor_chain_reverse_order :
            l_ancestor_chain.append( ancestor_chain_reverse_order[ : : -1 ] )
        return l_ancestor_chain
    def get_longest_shared_ancestor_chains( self, term_1, term_2 ) -> set :
        """
        return the ancestor chains to the most closest shared ancestors between the term1 and term2
        
        # 2024-02-29 22:58:13 
        """
        # retrieve ancestor chains
        l_ancestor_chain_1 = self.get_ancestor_chain( term_1 )
        l_ancestor_chain_2 = self.get_ancestor_chain( term_2 )
        
        # collect the ancestor chains to the most closest shared ancestors
        set_ancestor_chain_to_most_closest_shared_ancestor = set( ) # initialize 'set_ancestor_chain_to_most_closest_shared_ancestor'
        for ancestor_chain_1 in l_ancestor_chain_1 : # iterate over chain list #1
            # find the chain in the chain list # 2 that contains the longest shared chain with the chain in the chain list #1
            l_index_most_closest_shared_ancestor = [ ]
            for ancestor_chain_2 in l_ancestor_chain_2 : # iterate over chain list #2
                index_most_closest_shared_ancestor = 0 # initialize the index that indicate the location of the most closest shared ancestor between the two chains # initialize with the index of the root term 
                for ancestor_1, ancestor_2 in zip( ancestor_chain_1[ 1 : ], ancestor_chain_2[ 1 : ] ) : # retrieve ancester from chain #1 and chain #2 (from the most distant ancestor (excluding the root term) to the closest ancestor) 
                    if ancestor_1 != ancestor_2 : # if the ancestors diverged between chain_1 and chain_2
                        break
                    index_most_closest_shared_ancestor += 1 # increase the pointer (take into account the current shared ancestor)
                l_index_most_closest_shared_ancestor.append( index_most_closest_shared_ancestor )
            index_chain_2_with_most_closest_shared_ancestor = np.argmax( l_index_most_closest_shared_ancestor ) # define 'most closest shared ancestor' as the ancestor that has the the largest number of ancestors between itself and the root term.
            ancestor_chain_to_most_closest_shared_ancestor = tuple( l_ancestor_chain_2[ index_chain_2_with_most_closest_shared_ancestor ][ : l_index_most_closest_shared_ancestor[ index_chain_2_with_most_closest_shared_ancestor ] + 1 ] ) # including the most_closest_shared_ancestor in the chain
            set_ancestor_chain_to_most_closest_shared_ancestor.add( ancestor_chain_to_most_closest_shared_ancestor )
        return set_ancestor_chain_to_most_closest_shared_ancestor
    def get_properties( self, term ) :
        """
        return the properties of the termk
        # 2024-03-01 12:51:18 
        """
        # get ontology terms
        term = self[ term ]
        # retrieve properties    
        l_label, l_comment, l_broadsynonym, l_exactsynonym = list( set( term.label ) ), list( set( term.comment ) ), list( set( term.hasBroadSynonym ) ), list( set( term.hasExactSynonym ) )
        def _parse_property( l ) :
            return l[ 0 ] if len( l ) > 0 else None
        dict_property = { 'label' : _parse_property( l_label ), 'comment' : _parse_property( l_comment ), 'broad_synonym' : l_broadsynonym, 'exact_synonym' : l_exactsynonym }
        return dict_property