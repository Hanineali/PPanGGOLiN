#!/usr/bin/env python3

# default libraries
from collections import defaultdict
from typing import Dict, Generator, List, Tuple, Union, Optional

#from ppanggolin.formats.writeAnnotations import intergenic_sequences_desc
from ppanggolin.genome import Gene, Organism, RNA, Intergenic


class Edge:
    """
    The Edge class represents an edge between two gene families in the pangenome graph. It is associated with all the
    organisms in which the neighborship is found, and all the involved genes as well.

    Methods:
        - get_org_dict: Returns a dictionary with organisms as keys and an iterable of the pairs in genes as values.
        - gene_pairs: Returns a list of all the gene pairs in the Edge.
        - add_genes: Adds genes to the edge. They are supposed to be in the same organism.

    Fields:
        - source: A GeneFamily object representing the source gene family of the edge.
        - target: A GeneFamily object representing the target gene family of the edge.
        - organisms: A defaultdict object representing the organisms in which the edge is found and the pairs of genes involved.
    """

    def __init__(self, source_feature: Union[Gene,RNA], target_feature: Union[Gene,RNA]):
        """Constructor method

        :param source_feature: First feature to initialize the edge
        :param target_feature: Second feature to initialize the edge
        """
        assert isinstance(source_feature, (Gene, RNA)) and isinstance(target_feature, (Gene, RNA)), \
            "Both source and target must be Gene or RNA objects."

        # TODO try to change for gene family ?
        if source_feature.family is None:
            raise AttributeError(
                f"You cannot create a graph without gene families. "
                f"gene {source_feature.ID} did not have a gene family."
            )
        if target_feature.family is None:
            raise AttributeError(
                f"You cannot create a graph without gene families. "
                f"gene {target_feature.ID} did not have a gene family."
            )
        self.source = source_feature.family
        self.target = target_feature.family
        src, tgt = sorted((self.source.name, self.target.name))
        self.name = f"{src} | {tgt}"
        #self.intergenics=[]
        self.source.set_edge(self.target, self)
        self.target.set_edge(self.source, self)
        self._organisms: Dict[Organism, Dict[str, List]] = defaultdict(lambda : {"pairs": [], "intergenic": []})
        self.add_features(source_feature, target_feature)

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """Get all the organisms belonging to the edge

        :return: Generator with organisms as the key and an iterable of the gene pairs as value
        """
        yield from self._organisms.keys()

    @property
    def number_of_organisms(self) -> int:
        """Get the number of organisms in the edge

        :return: Number of organisms
        """
        return len(self._organisms)

    def get_organism_feature_pairs(self, organism: Organism) -> List[Tuple[Gene, Gene]]:
        """Get the gene pair corresponding to the given organism

        :param organism: Wanted organism

        :return: Pair of genes in the edge corresponding to the given organism
        """
        return self._organisms[organism]["pairs"]

    def get_organisms_dict(self) -> Dict[Organism, List[Tuple[Gene, Gene]]]:
        """Get all the organisms with their corresponding pair of genes in the edge

        :return: Dictionary with the organism as the key and list of gene pairs as value
        """
        return self._organisms

    @property
    def feature_pairs(self) -> List[Tuple[Union[Gene, RNA], Union[Gene, RNA]]]:
        """Get the list of all the feature pairs in the Edge

        :return: A list of all the feature pairs in the Edge
        """
        return [
            feature_pair
            for org_data in self._organisms.values()
            for feature_pair in org_data["pairs"]
        ]

    def add_features(self, source_feature, target_feature):
        """
        Adds features (genes/rnas) to the edge.
        They are supposed to be in the same organism.

        :param source_feature: Feature corresponding to the source of the edge
        :param target_feature: Feature corresponding to the target of the edge

        :raises TypeError: If the features are not with Gene/RNA type
        :raises ValueError: If features are not associated with an organism
        :raises Exception: If the features are not in the same organism.
        """
        if not isinstance(source_feature, (Gene,RNA)) or not isinstance(target_feature, (Gene,RNA)):
            raise TypeError(
                f"Genes/RNAs are expected to be added to edge. "
                f"Given type for source: {type(source_feature)} and target: {type(target_feature)}"
            )
        if source_feature.organism is None or target_feature.organism is None:
            raise ValueError(
                "Genes/RNAs are not associated to genome. It's needed to create add genes/rnas to edge"
            )
        if source_feature.organism != target_feature.organism:
            raise Exception(
                f"You tried to create an edge between two genes/rnas that are not even in the same genome ! "
                f"(features are '{source_feature.ID}' and '{target_feature.ID}')"
            )

        self._organisms[source_feature.organism]["pairs"].append((source_feature, target_feature))

    @property
    def intergenic_chain(self) -> List[Tuple[Union[Gene, RNA], Union[Gene, RNA]]]:
        """Get the list of all the intergenic chains in the Edge

        :return: A list of all the intergenic's chain in the Edge
        """
        return [
            intergenic_chain
            for org_data in self._organisms.values()
            for intergenic_chain in org_data["intergenic"]
        ]

    def add_intergenic_chain(self, intergenic_chain: Tuple[Union[Intergenic, Tuple[Union[Intergenic, Gene, RNA], ...]]]):
        """
        Adds the corresponding intergenic (if no removed family) / intergenic_feature chain (if removed family).
        :param intergenic_chain: A list containing Intergenic, or a tuple of (Intergenic, Gene/RNA, Intergenic), etc.
        """
        if not intergenic_chain:
            return

        #  Identify the organism from the first valid item in the chain
        #  It might be a tuple or a single element. We just need to find an object that has an organism
        org = None

        # function to extract an organism from a chain item
        def get_organism_from_item(item) -> Optional[Organism]:
            """Return the organism from item if it exists, or from sub-items if item is a tuple."""
            if isinstance(item, tuple):
                for sub in item:
                    # Check if sub has an organism attribute
                    if hasattr(sub, 'organism') and sub.organism is not None:
                        return sub.organism
            else:
                if hasattr(item, 'organism') and item.organism is not None:
                    return item.organism
            return None

        for element in intergenic_chain:
            org = get_organism_from_item(element)
            if org is not None:
                break

        if org is None:
            return

        if org not in self._organisms:
            self._organisms[org] = {"pairs": [], "intergenic": []}

        self._organisms[org]['intergenic'].append(intergenic_chain)


