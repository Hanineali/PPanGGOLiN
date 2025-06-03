#!/usr/bin/env python3

# default libraries
import logging
import argparse
from pathlib import Path

# installed libraries
from tqdm import tqdm

from ppanggolin import pangenome
from ppanggolin.genome import Organism, Gene, Intergenic, Contig
# local libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import read_pangenome, write_pangenome, erase_pangenome


def check_pangenome_former_graph(pangenome: Pangenome, force: bool = False):
    """
    Checks pangenome status and .h5 files for former neighbors graph, delete it if allowed or raise an error

    :param pangenome: Pangenome object
    :param force: Allow to force write on Pangenome file
    """
    if pangenome.status["neighborsGraph"] == "inFile" and not force:
        raise AttributeError(
            "You are trying to make a neighbors graph that is already built. "
            "If you REALLY want to do that, use --force "
            "(it will erase everything except annotation data !)"
        )
    elif pangenome.status["neighborsGraph"] == "inFile" and force:
        erase_pangenome(pangenome, graph=True)


def check_pangenome_for_neighbors_graph(pangenome, force, disable_bar=False):
    """
    Checks and read the pangenome for neighbors graph computing.

    :param pangenome: Pangenome object
    :param force: Allow to force write on Pangenome file
    :param disable_bar: Disable progress bar
    """
    check_pangenome_former_graph(pangenome, force)
    # TODO Check if possible to change for check_pangenome_info
    if pangenome.status["genomesAnnotated"] in [
        "Computed",
        "Loaded",
    ] and pangenome.status["genesClustered"] in ["Computed", "Loaded"]:
        pass  # nothing to do, can just continue.
    elif (
        pangenome.status["genomesAnnotated"] == "inFile"
        and pangenome.status["genesClustered"] == "inFile"
    ):
        read_pangenome(
            pangenome, annotation=True, gene_families=True, rna_families=True, intergenic_sequences=True,disable_bar=disable_bar
        )
    elif pangenome.status["genesClustered"] == "No" and pangenome.status[
        "genomesAnnotated"
    ] in ["inFile", "Computed", "Loaded"]:
        raise Exception(
            "You did not cluster the genes. See the 'ppanggolin cluster' if you want to do that."
        )
    else:
        # You probably can use readPangenome anyway.
        msg = (
            "Dev : You are probably writing a new workflow with a combination that I did not test."
            " You can probably use readPangenome instead of raising this Error. "
            "However please test it carefully.\n"
        )
        msg += (
            " User : I have no idea how you got there. You probably did something unexpected. "
            "Post an issue with what you did at https://github.com/labgem/PPanGGOLiN\n"
        )
        raise NotImplementedError(msg)


def remove_high_copy_number(pangenome, number, rnas = False):
    """Removes families present more than 'number' times from the pangenome graph

    :param pangenome: Pangenome object
    :param number: Maximum authorized repeat presence
    :param rnas: If True, also applies to RNA families
    """
    for fam in pangenome.gene_families:
        for gene_list in fam.get_org_dict().values():
            if len(gene_list) >= number:
                fam.removed = True

        # RNA families
    if rnas:
        for fam in pangenome.rna_families:
            for rna_list in fam.get_org_dict().values():
                if len(rna_list) >= number:
                    fam.removed = True

""" flag edge as linked fragment """
def flag_fragment_edge(edge, feature1, feature2):
    """
    If both features (RNA/Gene) belong to the same family and at least one is a fragment,
    mark this edge as linked_fragment = True.
    """
    if feature1.family == feature2.family:
        if feature1.is_fragment or feature2.is_fragment:
            edge.is_linked_fragment = True

def assign_edge_to_intergenic_data(items, edge):
    """
    Iterates over items which can be:
      - single Intergenic
      - removed Gene/RNA
      - a tuple of (Intergenic, Gene, Intergenic)
      - a nested structure with multiple items
    and assigns .edge = edge to any Intergenic found.
    """
    for element in items:
        if isinstance(element, tuple):
            # Potentially (Intergenic, Gene, Intergenic) or more
            for subelement in element:
                if isinstance(subelement, Intergenic):
                    subelement.edge = edge
        else:
            # Could be a single Intergenic, Gene, or RNA
            if isinstance(element, Intergenic):
                element.edge = edge


def compute_neighbors_graph(
    pangenome: Pangenome,
    remove_copy_number: int = 0,
    force: bool = False,
    disable_bar: bool = False,
):
    """
        Creates the Pangenome Graph. Will either load the information from the pangenome file if they are not loaded,
        or use the information loaded if they are.

        :param pangenome: Pangenome object
        :param remove_copy_number: Maximum authorized repeat presence of gene families. if zero no remove
        :param force: Allow to force write on Pangenome file
        :param disable_bar: Disable progress bar
        """
    check_pangenome_for_neighbors_graph(pangenome, force, disable_bar=disable_bar)

    if remove_copy_number > 0:
        remove_high_copy_number(pangenome, remove_copy_number)

    logging.getLogger("PPanGGOLiN").info("Computing the neighbors graph...")
    bar = tqdm(
        pangenome.organisms,
        total=pangenome.number_of_organisms,
        unit="genome",
        disable=disable_bar,
    )
    for org in bar:
        bar.set_description(f"Processing {org.name}")

        for contig in org.contigs:
            all_features = sorted(list(contig.genes) + list(contig.RNAs), key=lambda x: x.start)

            prev_conserved_feat = None
            temp_list = []

            for idx, feature in enumerate(all_features):
                try:
                    intergenic = None
                    if idx > 0:
                        prev_feat = all_features[idx - 1]
                        intergenic = contig.get_intergenic(
                            prev_feat, feature
                        )
                    # If the family is removed, store in temp_list
                    if feature.family.removed:
                        if intergenic:
                            temp_list.append(intergenic)
                        temp_list.append(feature)
                        continue

                    if prev_conserved_feat is not None:
                        # Create edge
                        feat_edge = pangenome.add_edge(prev_conserved_feat, feature)

                        # Attach any intergenic or removed features
                        if intergenic:
                            temp_list.append(intergenic)

                        if temp_list:
                            feat_edge.add_intergenic(tuple(temp_list))
                            assign_edge_to_intergenic_data(temp_list, feat_edge)

                        # Flag if this edge links two fragment features of the same family
                        flag_fragment_edge(feat_edge, prev_conserved_feat, feature)

                    # Reset
                    temp_list = []
                    prev_conserved_feat = feature

                except AttributeError:
                    raise AttributeError(
                        "A feature does not have a GeneFamily object associated."
                    )
                except Exception:
                    raise Exception("Unexpected error. Please report on our github.")

                # Handle circular contigs
            if contig.is_circular:
                # e.g., find first and last conserved features
                first_feat = next((f for f in all_features if not f.family.removed), None)
                last_feat = prev_conserved_feat

                if first_feat and last_feat:
                    circular_chain = []

                    # For leftover features after last_feat in the list
                    idx_last = all_features.index(last_feat)
                    if idx_last < len(all_features) - 1:
                        for feat in all_features[idx_last + 1:]:
                            try:
                                intergenic = contig.get_intergenic(last_feat, feat)
                                if intergenic:
                                    circular_chain.append(intergenic)
                                circular_chain.append(feat)
                                last_feat = feat
                            except AttributeError:
                                raise AttributeError(
                                    "A feature in circular chain has no Family object."
                                )
                            except Exception:
                                raise Exception(
                                    "Unexpected error during circular chain handling."
                                )

                    # Intergenic between last and first
                    intergenic = contig.get_intergenic(all_features[-1], all_features[0])
                    if intergenic:
                        circular_chain.append(intergenic)

                    # For features before first_feat in the list
                    idx_first = all_features.index(first_feat)
                    if idx_first > 0:
                        for feat in all_features[:idx_first]:
                            circular_chain.append(feat)
                            intergenic = contig.get_intergenic(feat, first_feat)
                            if intergenic:
                                circular_chain.append(intergenic)

                    # Create circular edge
                    try:
                        feat_edge = pangenome.add_edge(prev_conserved_feat, first_feat)
                        if circular_chain:
                            feat_edge.add_intergenic(tuple(circular_chain))
                            assign_edge_to_intergenic_data(circular_chain, feat_edge)

                        # Flag for fragment link if needed
                        flag_fragment_edge(feat_edge, prev_conserved_feat, first_feat)

                    except AttributeError:
                        raise AttributeError(
                            "An attribute error occurred while linking circular contig edges."
                        )
                    except Exception:
                        raise Exception(
                            "Unexpected error while creating circular edges."
                        )

    logging.getLogger("PPanGGOLiN").info(
        "Done making the neighbors graph with all features genes,rnas and intergenics.")
    pangenome.status["neighborsGraph"] = "Computed"
    pangenome.parameters["graph"] = {}
    if remove_copy_number > 0:
        pangenome.parameters["graph"]["remove_high_copy_number"] = remove_copy_number

def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    pangenome = Pangenome()
    pangenome.add_file(args.pangenome)
    compute_neighbors_graph(
        pangenome,
        args.remove_high_copy_number,
        args.force,
        disable_bar=args.disable_prog_bar,
    )
    write_pangenome(
        pangenome, pangenome.file, args.force, disable_bar=args.disable_prog_bar
    )

def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PPanGGOLiN in Command line

    :param sub_parser : sub_parser for graph command

    :return : parser arguments for graph command
    """
    parser = sub_parser.add_parser(
        "graph", formatter_class=argparse.RawTextHelpFormatter
    )
    parser_graph(parser)
    return parser

def parser_graph(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of graph command

    :param parser: parser for graph argument
    """
    required = parser.add_argument_group(
        title="Required arguments", description="Following arguments is required:"
    )
    required.add_argument(
        "-p", "--pangenome", required=False, type=Path, help="The pangenome .h5 file"
    )
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument(
        "-r",
        "--remove_high_copy_number",
        type=int,
        default=0,
        help="Positive Number: Remove families having a number of copy of gene in a single genome "
        "above or equal to this threshold in at least one genome "
        "(0 or negative values are ignored).",
    )


if __name__ == "__main__":
    """To test local change and allow using debugger"""
    from ppanggolin.utils import set_verbosity_level, add_common_arguments

    main_parser = argparse.ArgumentParser(
        description="Depicting microbial species diversity via a Partitioned PanGenome Graph Of Linked Neighbors",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser_graph(main_parser)
    add_common_arguments(main_parser)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
