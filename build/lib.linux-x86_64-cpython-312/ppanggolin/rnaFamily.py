#!/usr/bin/env python3

# default libraries
from __future__ import annotations
from collections import defaultdict
import logging

# installed libraries
from typing import Dict, Generator, Set
import gmpy2

# local libraries
from ppanggolin.edge import Edge
from ppanggolin.genome import Gene, Organism, RNA
from ppanggolin.metadata import MetaFeatures


class rnaFamily(MetaFeatures):
    """
    Represents a single RNA family in the pangenome.

    Methods:
        - neighbors: Returns all the RNAFamilies linked with an edge.
        - edges: Returns all Edges linked to this RNA family.
        - rnas: Returns all the RNA genes associated with the family.
        - organisms: Returns all the Organisms that have this RNA family.
        - set_edge: Sets an edge between the current family and a target family.
        - add_rna: Adds an RNA gene to the RNA family and sets the gene's family accordingly.
        - get_org_dict: Returns a dictionary of organisms as keys and sets of RNA genes as values.

    Fields:
        - name: The name of the RNA family.
        - ID: The internal identifier of the RNA family.
        - removed: Boolean indicating whether the family has been removed from the main graph.
        - _representative: The representative RNA sequence.
        - _rnas_getter: Dictionary storing RNA genes in the family.

    """

    def __init__(self, family_id: int, name: str):
        # TODO edges as genes in contig to get and set
        """Constructor method
        :param family_id: The internal identifier to give to the gene family
        :type family_id: any
        :param name: The name of the gene family (to be printed in output files)
        :type name: str
        """
        assert isinstance(family_id, int), "rnaFamily object id should be an integer"
        assert isinstance(name, str), "rnaFamily object name should be a string"
        assert name != "", "rnaFamily object cannot be created with an empty name"

        super().__init__()
        self.name = str(name)
        self.ID = family_id
        self._representative = None
        self._edges_getter = {}
        self._rnaPerOrg = defaultdict(set)
        self._rnas_getter = {}
        self.removed = False  # for the repeated family not added in the main graph
        self.sequence = ""
        self._family_type = None  # Type of RNA family (e.g., "tRNA", "tmRNA", "ncRNA")

    def __repr__(self) -> str:
        """Family representation"""
        return f"{self.ID}: {self.name}"

    def __len__(self) -> int:
        """Get the number of rnas in the family

        :return: The length of the list of rnas
        """
        return len(self._rnas_getter)

    def __setitem__(self, identifier: str, rna: RNA):
        """Set RNA gene to RNA Family

        :param identifier: ID of the RNA gene
        :param rna: RNA object to add
        """
        if not isinstance(rna, RNA):
            raise TypeError(f"'RNA' type was expected but got '{type(rna)}'")
        if not isinstance(identifier, str):
            raise TypeError(f"RNA ID should be a string, got '{type(identifier)}'")
        if identifier in self._rnas_getter:
            raise KeyError(f"RNA with ID {identifier} already exists in the RNA family")

        self._rnas_getter[identifier] = rna


    # retrieve gene by start position
    def __getitem__(self, identifier: str) -> RNA:
        """Get the RNA gene for the given name

        :param identifier: ID of the RNA gene

        :return: RNA gene object
        """
        if not isinstance(identifier, str):
            raise TypeError(f"RNA ID should be a string, got '{type(identifier)}'")
        try:
            return self._rnas_getter[identifier]
        except KeyError:
            raise KeyError(f"RNA with ID {identifier} does not exist in the family")

    def __delitem__(self, identifier: str):
        """Remove the gene for the given name in the gene family

        :param position: ID of the gene in the family

        :raises TypeError: If the identifier is not instance string
        :raises KeyError: Gene with the given identifier does not exist in the contig
        """
        if not isinstance(identifier, str):
            raise TypeError(
                f"Gene ID should be a string. You provided a '{type(identifier)}' type object"
            )
        try:
            del self._genes_getter[identifier]
        except KeyError:
            raise KeyError(
                f"Gene with the name: {identifier} does not exist in the family"
            )

    def add(self, rna: RNA):
        """Add an RNA gene to the RNA family, and set the RNA's family accordingly.

        :param rna: The RNA to add
        """
        if not isinstance(rna, RNA):
            raise TypeError(f"Expected 'RNA' type, but got '{type(rna)}'")

        self[rna.ID] = rna
        rna.family = self

        if rna.organism is not None and rna.organism in self._rnaPerOrg:
            self._rnaPerOrg[rna.organism].add(rna)

    def get_org_dict(self) -> Dict[Organism, Set[RNA]]:
        """Returns organisms and their associated RNA genes in the family

        :return: A dictionary {Organism: Set of RNA genes}
        """
        if len(self._rnaPerOrg) == 0:
            for rna in self._rnas_getter.values():
                if rna.organism is None:
                    raise AttributeError(f"RNA {rna.ID} is not linked to an organism")
                self._rnaPerOrg[rna.organism].add(rna)

        return self._rnaPerOrg

    def remove(self, identifier: str):
        """Remove an RNA gene by its name (ID)

        :param identifier: Name (ID) of the RNA gene to remove

        :raises TypeError: If the identifier is not a string
        :raises KeyError: If the RNA gene does not exist in the family
        """
        if not isinstance(identifier, str):
            raise TypeError(f"RNA ID should be a string. You provided a '{type(identifier)}' type object")

        try:
            rna_to_remove = self._rnas_getter[identifier]
            del self._rnas_getter[identifier]  # Remove RNA from the family

            # Also remove from the organism dictionary
            if rna_to_remove.organism and rna_to_remove.organism in self._rnaPerOrg:
                self._rnaPerOrg[rna_to_remove.organism].discard(rna_to_remove)

        except KeyError:
            raise KeyError(f"RNA with ID '{identifier}' does not exist in the RNA family")

    @property
    def representative(self) -> RNA:
        """Get the representative RNA of the family"""
        if self._representative is None:
            raise Exception("Representative RNA has not been set")
        return self._representative

    @representative.setter
    def representative(self, rna: RNA) -> None:
        """Set the representative RNA of the family"""
        if not isinstance(rna, RNA):
            raise TypeError(f"Representative RNA should be of type 'RNA', got '{type(rna)}'")
        self._representative = rna

    @property
    def rnas(self):
        """Return all the RNA genes belonging to the family"""
        yield from self._rnas_getter.values()

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """Returns all the Organisms that have this RNA family"""
        if len(self._rnaPerOrg) == 0:
            _ = self.get_org_dict()
        yield from self._rnaPerOrg.keys()

    def set_edge(self, target: rnaFamily, edge):
        """Set an edge between RNA families

        :param target: Neighbor RNA family
        :param edge: Edge connecting families
        """
        self._edges_getter[target] = edge

    def get_edge(self, target: rnaFamily):
        """Get the edge linking this RNA family to another"""
        return self._edges_getter[target]

    def add_sequence(self, seq: str):
        """Assigns a protein sequence to the gene family.

        :param seq: The sequence to add to the gene family
        """
        assert isinstance(seq, str), "Sequence must be a string"

        self.sequence = seq






