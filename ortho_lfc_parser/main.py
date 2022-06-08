from pathlib import Path
from typing import Generator, List, Any, Dict
import re
import argparse
import logging

logging.basicConfig(format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

fix_protein_regex = re.compile(r'^([a-zA-Z\d]+[_.][a-zA-Z\d]+)')


def fix_protein_names(proteins: Generator[str, Any, None]) -> List[str]:
    """
    Fix the given protein names as:
    HmN_000011500.1 -> HmN_000011500
    MCOS_0000277701-mRNA-1 -> MCOS_0000277701
    MSTRG.861.1 -> MSTRG.861
    """
    fixed_proteins = []
    for protein in proteins:
        match = fix_protein_regex.search(protein)
        if match:
            fixed_name = match.group(1)
            fixed_proteins.append(fixed_name)

    return fixed_proteins


def process_csv(orthogroup_filepath: str, outdir: Path) -> Path:
    """
    Processes Orthogroups.csv keeping for each gene only species gene and removing everything else.
    e.g.:
    HmN_000011500.1 -> HmN_000011500
    MCOS_0000277701-mRNA-1 -> MCOS_0000277701
    Emu_MSTRG.861.1 -> MSTRG.861
    """
    logging.info(f"Processing {orthogroup_filepath}")
    processed_csv = outdir / "processed_orthogroups.csv"
    with open(orthogroup_filepath) as orthogroup_fh, \
            open(processed_csv, "w") as processed_orthogroup_fh:
        processed_orthogroup_fh.write(orthogroup_fh.readline())  # print the header

        for line in orthogroup_fh:
            line_split = line.strip("\r\n").split("\t")
            orthogroup = line_split[0]
            proteins_list = []
            for index, proteins_as_str in enumerate(line_split[1:]):
                proteins = (protein.strip() for protein in proteins_as_str.split(", "))
                fixed_proteins = fix_protein_names(proteins)
                proteins_list.append(", ".join(fixed_proteins))
            processed_orthogroup_fh.write("\t".join([orthogroup] + proteins_list) + "\n")

    return processed_csv

def get_organisms(orthogroup_file: str) -> List[str]:
    """
	parse the header of orthogroup_file to get the order of the organisms
	"""
    with open(orthogroup_file) as fh:
        header = fh.readline()
        organisms = header.strip().split()
    return organisms


get_everything_before_the_dot_RE = re.compile(r"^([a-zA-Z\d]+[_.][a-zA-Z\d]+)")


def parse_orthogroup_file(orthogroup_filepath: str, organisms: List[str]) -> Dict[str, Dict[str, List[str]]]:
    logging.info(f"Parsing {orthogroup_filepath}...")
    orthogroup2organism2gene_list = {}
    with open(orthogroup_filepath) as orthogroup_fh:
        orthogroup_fh.readline()  # skip header

        # populate orthogroup2organism2gene_list
        for line in orthogroup_fh:
            line_split = line.strip("\r\n").split("\t")
            orthogroup = line_split[0]
            orthogroup2organism2gene_list[orthogroup] = {}
            one_of_the_columns_is_empty = False  # a flag to indicate if one of the columns is empty - we want to process only when all columns are filled

            for idx, organism in enumerate(organisms):
                all_genes_as_str = line_split[idx + 1]
                if all_genes_as_str == "":
                    one_of_the_columns_is_empty = True
                    break
                else:
                    orthogroup2organism2gene_list[orthogroup][organism] = [geneWithPoint for geneWithPoint in
                                                                           all_genes_as_str.split(", ")]

            if one_of_the_columns_is_empty:  # we do not want this
                del orthogroup2organism2gene_list[orthogroup]

    return orthogroup2organism2gene_list


def read_differential_expression_file(differential_expression_file_name) -> Dict[str, float]:
    """
    Read the differential expression file (2 columns, gene and log2foldchange)
    @return: a dict{gene: log2foldchange}
    """
    logging.info(f"Reading DE file {differential_expression_file_name}...")
    gene2LFC = {}
    with open(differential_expression_file_name) as differential_expression_fh:
        for line in differential_expression_fh:
            line = line.strip()
            if line != "":
                if "," in line:
                    line_split = line.split(",")
                else:
                    line_split = line.split()
                gene_with_point = line_split[0]
                gene = gene_with_point
                gene2LFC[gene] = float(line_split[1])
    return gene2LFC


def write_combination_to_file(fh, all_events_fh, orthogroup, combination, gene2LFCs):
    line_to_write = [orthogroup]
    for gene in combination:
        line_to_write.append("_" + gene)
    for i, gene in enumerate(combination):
        line_to_write.append("\t" + str(gene2LFCs[i][gene]))
    line_to_write = "".join(line_to_write)
    print(line_to_write, file=fh)
    print(line_to_write, file=all_events_fh)


def process_and_print_combination(orthogroup, combination, gene2LFCs, all_positives_fh, all_negatives_fh,
                                  positives_and_negatives_fh, all_events_fh):
    nb_of_positives = 0
    nb_of_negatives = 0
    for i, gene in enumerate(combination):
        if gene2LFCs[i][gene] > 0:
            nb_of_positives += 1
        elif gene2LFCs[i][gene] < 0:
            nb_of_negatives += 1
    if nb_of_negatives == 0:
        write_combination_to_file(all_positives_fh, all_events_fh, orthogroup, combination, gene2LFCs)
    elif nb_of_positives == 0:
        write_combination_to_file(all_negatives_fh, all_events_fh, orthogroup, combination, gene2LFCs)
    else:
        write_combination_to_file(positives_and_negatives_fh, all_events_fh, orthogroup, combination, gene2LFCs)


def generate_all_combinations_recursively(orthogroup, organisms, organism2gene_list, combination, gene2LFCs, i,
                                          all_positives_fh, all_negatives_fh, positives_and_negatives_fh,
                                          all_events_fh, not_significatives_genes_fh):
    organism = organisms[i]

    # add all genes to the combination
    for gene in organism2gene_list[organism]:
        if gene not in gene2LFCs[i].keys():
            not_significatives_genes_fh.write(gene + "\n")
            continue  # this gene was not found as significative, skip

        combination.append(gene)
        if i < len(organisms) - 1:
            # we are still building the combination, call it recursively
            generate_all_combinations_recursively(orthogroup, organisms, organism2gene_list, combination, gene2LFCs,
                                                  i + 1, all_positives_fh, all_negatives_fh, positives_and_negatives_fh,
                                                  all_events_fh, not_significatives_genes_fh)
        else:
            # we have a combination, call the desired function to work on it
            process_and_print_combination(orthogroup, combination, gene2LFCs, all_positives_fh, all_negatives_fh,
                                          positives_and_negatives_fh, all_events_fh)
        combination.pop()  # remove it from the list


def generate_all_combinations_and_process(orthogroup2organism2gene_list, gene2LFCs, organisms, all_positives_fh,
                                          all_negatives_fh, positives_and_negatives_fh,
                                          all_events_fh, not_significatives_genes_fh):
    logging.info("Generating all combinations and processing...")
    for orthogroup, organism2gene_list in orthogroup2organism2gene_list.items():
        combination = []
        generate_all_combinations_recursively(orthogroup, organisms, organism2gene_list, combination, gene2LFCs, 0,
                                              all_positives_fh, all_negatives_fh, positives_and_negatives_fh,
                                              all_events_fh, not_significatives_genes_fh)


def orthogroup_parse(orthogroup_file, diff_expr_files, outdir):
    organisms = get_organisms(orthogroup_file)
    assert len(diff_expr_files) == len(diff_expr_files), "The number of organisms in the orthogroup file must match the" \
                                                         "number of differential expression files"

    orthogroup_2_organism_2_gene_list = parse_orthogroup_file(orthogroup_file, organisms)

    gene_2_LFCs = []
    for differential_expression_file_name in diff_expr_files:
        gene_2_LFCs.append(read_differential_expression_file(differential_expression_file_name))

    with open(outdir/"allPositives.tsv", "w") as all_positives_fh, \
            open(outdir/"allNegatives.tsv", "w") as all_negatives_fh, \
            open(outdir/"positivesAndNegatives.tsv", "w") as positives_and_negatives_fh, \
            open(outdir/"allEvents.tsv", "w") as all_events_fh, \
            open(outdir/"not_significatives_genes.tsv", "w") as not_significatives_genes_fh:
        generate_all_combinations_and_process(orthogroup_2_organism_2_gene_list, gene_2_LFCs, organisms,
                                              all_positives_fh, all_negatives_fh, positives_and_negatives_fh,
                                              all_events_fh, not_significatives_genes_fh)

    logging.info("All done!")

def get_args():
    parser = argparse.ArgumentParser(description='Orthogroup LFC parser.')
    parser.add_argument('--orthogroups', help='The orthogroups file', required=True)
    parser.add_argument('--outdir', help='Output directory', default="out", required=False)
    parser.add_argument('DE_files', type=str, nargs='+', help='Differential expression files')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    processed_csv = process_csv(args.orthogroups, outdir)
    orthogroup_parse(processed_csv, args.DE_files, outdir)


if __name__ == "__main__":
    main()
