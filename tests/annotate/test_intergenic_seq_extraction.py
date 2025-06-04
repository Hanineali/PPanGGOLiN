""
from pickle import FALSE

import pytest
from unittest.mock import MagicMock, patch
from io import StringIO

# Import the function and dependencies (replace `your_module` with the actual module name)
from ppanggolin.annotate.synta import process_genes_intergenics_seq


@patch("ppanggolin.annotate.synta.create_intergenic")  # Mock intergenic region creation
@patch("ppanggolin.annotate.synta.get_dna_sequence")  # Mock DNA sequence retrieval

@pytest.mark.parametrize(
    "is_circular, genes, expected_calls",
    [
        # CASE 1: Linear Contig with One Intergenic Region
        (False, [(1, 35), (65, 80)], [((36, 64), "gene_1 | gene_2")]),

        # CASE 2: Circular Contig with Wrap-Around Intergenic
        (True, [(10, 25), (45, 60)], [((26, 44),"gene_1 | gene_2"), ([(61, 80),(1, 9)],"gene_2|gene_1")]),

        # CASE 3: Genes Touching (No Intergenic Region)
        (False, [(1, 50), (51, 80)], [((50, 51), "gene_1 | gene_2")]),

        # CASE 4: Overlapping Genes
        (False, [(1, 65), (62, 80)], [((62, 65), "gene_1 | gene_2")]),

        # CASE 5: Border Intergenic Only
        (False, [(20, 80)], [((1, 19), "|gene_1")]),

        # CASE 6: Empty Contig (No Genes)
        (False, [], []),

        # CASE 7: Circular Contig with Genes touching (No intergenic Region)
        (True, [(1, 25), (26, 80)], [((25, 26), 'gene_1 | gene_2'), ((80, 1), 'gene_2|gene_1')]),

        # CASE 8: Border end Intergenic
        (False, [(1, 55)], [((56,80), 'gene_1|')]),
    ],
)

def test_process_genes_intergenics_seq(mock_get_dna_sequence, mock_create_intergenic, is_circular, genes, expected_calls):
    """
    Test `process_genes_intergenics_seq` to ensure intergenic regions are correctly processed.

    Explore all cases:
    - Linear vs. Circular Contigs
    - Intergenic Regions Between Genes
    - Genes at the Start or End of Contigs
    - Overlapping Genes
    - Single Gene Contigs
    - Empty Contigs (No Genes)
    """
    # 1. Mock Contig
    mock_contig = MagicMock()
    mock_contig.name = "contig1"
    mock_contig.is_circular = is_circular

    # 2. Mock Contig Sequence
    contig_seq = "ATGCGTACGTAGCTAGCTAGCTGACTGACTGATCGATCGT" * 2  # Fake long DNA sequence

    # 3. Create Fake Genes
    mock_genes = []
    for idx, (start, stop) in enumerate(genes):
        mock_gene = MagicMock()
        mock_gene.ID = f"gene_{idx + 1}"
        mock_gene.start = start
        mock_gene.stop = stop
        mock_genes.append(mock_gene)

    # 4. Mock Organism
    mock_org = MagicMock()

    # 5. Mock DNA Extraction
    mock_get_dna_sequence.side_effect = lambda seq, gene: seq[gene.start - 1: gene.stop]

    # 6. Run Function
    process_genes_intergenics_seq(mock_contig, mock_genes, contig_seq, mock_org)

    # 7. Assertions - Check Intergenic Regions Created
    def normalize_coordinates(coord):
        """Ensure coordinates are correctly formatted as a tuple or tuple of tuples."""
        if isinstance(coord, list) and len(coord)==1:
            return coord[0]
        elif isinstance(coord, list) and len(coord) > 1:
            return coord

    actual_calls = [
        (normalize_coordinates(call.kwargs["coordinates"]), call.kwargs["intergenic_id"])
        for call in mock_create_intergenic.call_args_list
    ]

    # Print for Debugging
    print(f" CASE: {'Circular' if is_circular else 'Linear'}, Genes: {genes}")
    print(f" Expected: {expected_calls}")
    print(f" Actual: {actual_calls}")

    assert actual_calls == expected_calls, f"Mismatch in expected intergenic regions! Expected: {expected_calls}, Got: {actual_calls}"
    """"""""