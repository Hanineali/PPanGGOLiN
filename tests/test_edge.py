#! /usr/bin/env python3

import pytest
from typing import Generator, Tuple, Union

from ppanggolin.genome import Gene, Organism, RNA, Intergenic
from ppanggolin.edge import Edge
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.rnaFamily import rnaFamily


class TestEdge:
    @pytest.fixture
    def organism(self) -> Generator[Organism, None, None]:
        """Generate a basic organism object"""
        yield Organism("organism")

    @pytest.fixture
    def families_pair(self) -> Generator[Tuple[GeneFamily, rnaFamily], None, None]:
        """Generate a pair of one GeneFamily and one rnaFamily."""
        yield GeneFamily(1, "family1"), rnaFamily(2, "family2")

    @pytest.fixture
    def feature_pairs(
        self, organism, families_pair
    ) -> Generator[Tuple[Union[Gene,RNA], Union[Gene,RNA]], None, None]:
        """Generate genes_pair"""
        gene1, rna2 = Gene("gene1"), RNA("rna2")
        gene1.fill_parents(organism, None)
        rna2.fill_parents(organism, None)
        gene1.family, rna2.family = GeneFamily(1, "family1"), rnaFamily(2, "family2")
        yield gene1, rna2

    @pytest.fixture
    def edge(self, feature_pairs):
        """Generate a basic edge"""
        edge = Edge(*feature_pairs)
        yield edge

    def test_constructor(self, feature_pairs, organism, families_pair):
        """Tests that an Edge object can be created with two genes belonging to different families"""
        gene1, rna2 = feature_pairs
        edge = Edge(gene1, rna2)
        assert edge.source == gene1.family
        assert edge.target == rna2.family
        assert edge.source._edges_getter[edge.target] == edge
        assert edge.target._edges_getter[edge.source] == edge
        assert set(edge._organisms.keys()) == {organism}
        assert edge._organisms[organism]["pairs"] == [(gene1, rna2)]
        assert edge._organisms[organism]["intergenic"] == []

    def test_constructor_attribute_error(self):
        """
        Tests that an AttributeError is raised when creating an Edge object
        with a gene that does not belong to any family
        """
        gene1 = Gene("gene1")
        gene1.family = GeneFamily(0, "test")
        gene2 = Gene("gene2")
        with pytest.raises(AttributeError):
            # Test target attribute error
            Edge(gene1, gene2)
        with pytest.raises(AttributeError):
            # Test source attribute error
            Edge(gene2, gene1)

    def test_feature_pairs(self, edge, feature_pairs):
        """Tests that gene pairs' generator return what's expected"""
        assert set(edge.feature_pairs) == {feature_pairs}

    def test_get_organisms(self, edge, organism):
        """Tests that organism generator return what's expected"""
        assert set(edge.organisms) == {organism}

    def test_get_number_of_organisms(self, edge):
        """Tests that the good number of organism is returned"""
        assert isinstance(edge.number_of_organisms, int)
        assert edge.number_of_organisms == 1

    def test_get_organisms_dict(self, edge, organism, feature_pairs):
        """Tests that get_organisms_dict returns the expected nested structure"""
        org_dict = edge.get_organisms_dict()
        assert set(org_dict.keys()) == {organism}
        assert isinstance(org_dict[organism], dict)
        assert org_dict[organism]["pairs"] == [feature_pairs]
        assert org_dict[organism]["intergenic"] == []

    def test_get_organism_features_pairs(self, edge, organism, feature_pairs):
        """Tests that the gene pairs corresponding to the organism is returned"""
        assert edge.get_organism_feature_pairs(organism) == [feature_pairs]

    def test_edge_add_genes_same_organism(self, edge, feature_pairs, organism):
        """Tests that genes can be added to the edge that are on the same organism"""
        gene1, gene2, gene3, gene4 = *feature_pairs, Gene("gene3"), RNA("gene4")
        gene3.fill_parents(organism, None)
        gene4.fill_parents(organism, None)
        edge.add_features(gene3, gene4)
        assert edge.get_organism_feature_pairs(organism) == [
            (gene1, gene2),
            (gene3, gene4),
        ]

    def test_edge_add_genes_different_organisms(self, edge, organism):
        """Tests that an Exception is raised when adding genes to the edge that are not on the same organism"""
        gene1, gene2 = Gene("gene3"), Gene("gene4")
        gene1.fill_parents(organism, None)
        org = Organism("org")
        gene2.fill_parents(org, None)
        with pytest.raises(Exception):
            edge.add_features(gene1, gene2)

    def test_edge_add_genes_one_none_gene(self, edge, organism):
        """Tests that a TypeError is raised when adding genes to the edge where one gene is None"""
        gene1 = Gene("gene1")
        gene1.fill_parents(organism)
        with pytest.raises(TypeError):
            edge.add_features(gene1, None)
        with pytest.raises(TypeError):
            edge.add_features(None, gene1)

    def test_edge_add_genes_without_organisms(self, edge, organism):
        """Tests that a ValueError is raised when adding genes not filled with organism"""
        gene1, gene2 = Gene("gene1"), Gene("gene2")
        gene1.fill_parents(organism, None)
        with pytest.raises(ValueError):
            edge.add_features(gene1, gene2)
        with pytest.raises(ValueError):
            edge.add_features(gene2, gene1)
