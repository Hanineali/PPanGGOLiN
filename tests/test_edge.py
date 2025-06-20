#! /usr/bin/env python3

import pytest
from typing import Generator, Tuple, Union

from ppanggolin.genome import Gene, Organism, Intergenic, RNA
from ppanggolin.edge import Edge
from ppanggolin.geneFamily import GeneFamily
from ppanggolin.rnaFamily import rnaFamily


class TestEdge:
    @pytest.fixture
    def organism(self) -> Generator[Organism, None, None]:
        yield Organism("organism")

    @pytest.fixture
    def families_pair(self) -> Generator[Tuple[Union[GeneFamily,rnaFamily], Union[GeneFamily,rnaFamily]], None, None]:
        yield GeneFamily(1, "family1"), rnaFamily(2, "rnafam2")

    @pytest.fixture
    def features_pair(self, organism, families_pair) -> Generator[Tuple[Union[Gene,RNA], Union[Gene,RNA]], None, None]:
        gene1, gene2 = Gene("gene1"), RNA("gene2")
        gene1.fill_parents(organism, None)
        gene2.fill_parents(organism, None)
        gene1.family, gene2.family = families_pair
        yield gene1, gene2

    @pytest.fixture
    def edge(self, features_pair):
        edge = Edge(*features_pair)
        yield edge

    def test_constructor(self, features_pair, organism):
        gene1, gene2 = features_pair
        edge = Edge(gene1, gene2)
        assert edge.source == gene1.family
        assert edge.target == rna2.family
        assert edge.source._edges_getter[edge.target] == edge
        assert edge.target._edges_getter[edge.source] == edge
        assert list(edge.organisms) == [organism]

    def test_constructor_attribute_error(self):
        gene1 = Gene("gene1")
        gene1.family = GeneFamily(0, "test")
        gene2 = Gene("gene2")
        with pytest.raises(AttributeError):
            Edge(gene1, gene2)
        with pytest.raises(AttributeError):
            Edge(gene2, gene1)

    def test_feature_pairs(self, edge, features_pair):
        assert set(edge.feature_pairs) == {features_pair}

    def test_get_organisms(self, edge, organism):
        assert set(edge.organisms) == {organism}

    def test_number_of_organisms(self, edge):
        assert isinstance(edge.number_of_organisms, int)
        assert edge.number_of_organisms == 1

    def test_get_organisms_dict(self, edge, organism, features_pair):
        expected = {organism: {"pairs": [features_pair], "intergenic": []}}
        assert edge.get_organisms_dict() == expected

    def test_get_organism_feature_pairs(self, edge, organism, features_pair):
        assert edge.get_organism_feature_pairs(organism) == [features_pair]

    def test_add_features_same_organism(self, edge, features_pair, organism):
        gene3, gene4 = Gene("gene3"), Gene("gene4")
        gene3.fill_parents(organism, None)
        gene4.fill_parents(organism, None)
        gene3.family = GeneFamily(3, "family3")
        gene4.family = GeneFamily(4, "family4")
        edge.add_features(gene3, gene4)
        assert edge.get_organism_feature_pairs(organism) == [features_pair, (gene3, gene4)]

    def test_add_features_different_organisms(self, edge, organism):
        gene1, gene2 = Gene("gene3"), Gene("gene4")
        gene1.fill_parents(organism, None)
        org2 = Organism("org2")
        gene2.fill_parents(org2, None)
        gene1.family = GeneFamily(3, "family3")
        gene2.family = GeneFamily(4, "family4")
        with pytest.raises(Exception):
            edge.add_features(gene1, gene2)

    def test_add_features_missing_organism(self, edge):
        gene1, gene2 = Gene("gene1"), Gene("gene2")
        gene1.family = GeneFamily(1, "f1")
        gene2.family = GeneFamily(2, "f2")
        with pytest.raises(ValueError):
            edge.add_features(gene1, gene2)

    def test_add_features_wrong_type(self, edge, organism):
        gene = Gene("gene")
        gene.fill_parents(organism, None)
        gene.family = GeneFamily(1, "f")
        with pytest.raises(TypeError):
            edge.add_features(gene, None)
        with pytest.raises(TypeError):
            edge.add_features(None, gene)

    def test_add_intergenic_chain(self, edge, organism):
        inter1 = Intergenic("int1")
        inter1.fill_parents(organism, None)
        edge.add_intergenic_chain((inter1,))
        assert edge._organisms[organism]['intergenic'] == [(inter1,)]

    def test_add_intergenic_chain_tuple_with_genes(self, edge, organism):
        inter1 = Intergenic("int1")
        inter2 = Intergenic("int2")
        gene = Gene("g")
        inter1.fill_parents(organism, None)
        inter2.fill_parents(organism, None)
        gene.fill_parents(organism, None)
        chain = (inter1, gene, inter2)
        edge.add_intergenic_chain(chain)
        assert edge._organisms[organism]['intergenic'][-1] == chain

    def test_add_empty_intergenic_chain(self, edge):
        # Should not raise or change anything
        edge.add_intergenic_chain(())
        assert all(len(data['intergenic']) == 0 for data in edge._organisms.values())
