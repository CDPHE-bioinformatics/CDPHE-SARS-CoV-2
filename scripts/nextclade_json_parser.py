"""
Parses a nextclade json file and outputs nextclade results and variant summary.
"""

import argparse
import json
import logging
import re
import sys

import pandas as pd


logger = logging.getLogger(__name__)


def parse_args(args: list) -> argparse.Namespace:
    """
    Parses the command line arguments.

    :param args: the command line arguments
    :returns: the parsed arguments
    """
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument(
        "--log_level",
        help="the level to log at",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
        default="INFO",
    )
    parser.add_argument("--nextclade_json", help="the nextclade json file path")
    parser.add_argument("--project_name", help="the project name")
    parser.add_argument("--workflow_version", help="the workflow version")
    return parser.parse_args(args)


def setup_logging(log_level: str):
    """
    Sets up the logging.

    :param log_level: the log level
    """
    log_format = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(
        level=log_level,
        stream=sys.stdout,
        format=log_format,
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def extract_hsn(sample_name: str) -> str:
    """
    Gets the hsn from the sample name. The hsn is the first 10 digits of the
    sample name and begins with 2.

    :param sample_name: the sample name
    :returns: the hsn
    """
    hsn = ""
    if sample_name:
        match = re.search(r"2[0-9]{9}", sample_name)
        if match:
            hsn = match.group()
    return hsn


def extract_sample_name(fasta_header: str) -> str:
    """
    Parses a fasta header to get the sample name. The fasta header is expected
    to have the format: >CO-CDPHE-<sample_name>

    :param fasta_header: the fasta header
    :returns: the sample name
    """
    sample_name = ""
    if fasta_header:
        match = re.search(r"CO-CDPHE-([0-9a-zA-Z_\\-\\.]+)", fasta_header)
        if match:
            sample_name = match.group(1)
        else:
            sample_name = fasta_header
    return sample_name


def get_results(json_data: str) -> pd.DataFrame:
    """
    Gets the results from the nextclade json data.

    :param json_data: the nextclade json data
    :returns: the results
    """
    logger.info("Processing results for %s samples.", len(json_data['results']))

    results = pd.DataFrame({
        "fasta_header": [],
        "sample_name": [],
        "hsn": [],
        "nextclade": [],
        "total_nucleotide_mutations": [],
        "total_nucleotide_deletions": [],
        "total_nucleotide_insertions": [],
        "total_AA_substitutions": [],
        "total_AA_deletions": []
    })

    for i in range(len(json_data["results"])):
        sample_json_data = json_data["results"][i]
        sample_fasta_header = sample_json_data["seqName"]
        sample_sample_name = extract_sample_name(sample_fasta_header)
        sample_hsn = extract_hsn(sample_sample_name)

        logger.debug("Processing sample %s.", sample_fasta_header)

        if "clade" in sample_json_data.keys():
            results.loc[len(results)] = {
                "fasta_header": sample_fasta_header,
                "sample_name": sample_sample_name,
                "hsn": sample_hsn,
                "nextclade": sample_json_data["clade"],
                "total_nucleotide_mutations": sample_json_data["totalSubstitutions"],
                "total_nucleotide_deletions": sample_json_data["totalDeletions"],
                "total_nucleotide_insertions": sample_json_data["totalInsertions"],
                "total_AA_substitutions": sample_json_data["totalAminoacidSubstitutions"],
                "total_AA_deletions": sample_json_data["totalAminoacidDeletions"]
            }

    return results


def get_variant_summary(json_data: str) -> pd.DataFrame:
    """
    Gets the variant summary from the nextclade json data.

    :param json_data: the nextclade json data
    :returns: the variant summary
    """
    variant_summary = pd.DataFrame({
        "fasta_header": [],
        "sample_name": [],
        "hsn": [],
        "variant_name": [],
        "gene": [],
        "codon_position": [],
        "refAA": [],
        "altAA": [],
        "start_nuc_pos": [],
        "end_nuc_pos": []
    })

    logger.info("Processing variant summary for %s samples.",
                len(json_data['results']))

    for i in range(len(json_data["results"])):
        sample_json_data = json_data["results"][i]
        sample_fasta_header = sample_json_data["seqName"]
        sample_sample_name = extract_sample_name(sample_fasta_header)
        sample_hsn = extract_hsn(sample_sample_name)

        logger.debug("Processing sample %s.", sample_fasta_header)

        if "aaDeletions" in sample_json_data.keys():
            for aa_deletion in sample_json_data["aaDeletions"]:
                # Get info required for variant_name
                gene = aa_deletion["gene"]
                ref_aa = aa_deletion["refAA"]
                alt_aa = "del"
                codon_position = aa_deletion["codon"] + 1

                # Compose variant_name
                variant_name = f"{gene}_{ref_aa}{codon_position}{alt_aa}"

                # Append to running summary
                variant_summary.loc[len(variant_summary)] = {
                    "fasta_header": sample_fasta_header,
                    "sample_name": sample_sample_name,
                    "hsn": sample_hsn,
                    "variant_name": variant_name,
                    "gene": gene,
                    "codon_position": codon_position,
                    "refAA": ref_aa,
                    "altAA": alt_aa,
                    "start_nuc_pos": aa_deletion["codonNucRange"]["begin"],
                    "end_nuc_pos": aa_deletion["codonNucRange"]["end"]
                }

        if "insertions" in sample_json_data.keys():
            for insertion in sample_json_data["insertions"]:
                # Get info required for variant_name
                gene = ""
                ref_aa = "ins"
                alt_aa = insertion["ins"]
                codon_position = insertion["pos"] + 1

                # Compose variant_name
                variant_name = f"{gene}_{ref_aa}{codon_position}{alt_aa}"

                # Append to running summary
                variant_summary.loc[len(variant_summary)] = {
                    "fasta_header": sample_fasta_header,
                    "sample_name": sample_sample_name,
                    "hsn": sample_hsn,
                    "variant_name": variant_name,
                    "gene": gene,
                    "codon_position": codon_position,
                    "refAA": ref_aa,
                    "altAA": alt_aa,
                    "start_nuc_pos": insertion["pos"] + 1,
                    "end_nuc_pos": insertion["pos"] + 1 + len(insertion["ins"])
                }

        if "aaSubstitutions" in sample_json_data.keys():
            for aa_sub in sample_json_data["aaSubstitutions"]:
                # Get info required for variant_name
                gene = aa_sub["gene"]
                ref_aa = aa_sub["refAA"]
                if aa_sub["queryAA"] == "*":
                    alt_aa = "stop"
                else:
                    alt_aa = aa_sub["queryAA"]
                codon_position = aa_sub["codon"] + 1

                # Compose variant_name
                variant_name = f"{gene}_{ref_aa}{codon_position}{alt_aa}"

                # Append to running summary
                variant_summary.loc[len(variant_summary)] = {
                    "fasta_header": sample_fasta_header,
                    "sample_name": sample_sample_name,
                    "hsn": sample_hsn,
                    "variant_name": variant_name,
                    "gene": gene,
                    "codon_position": codon_position,
                    "refAA": ref_aa,
                    "altAA": alt_aa,
                    "start_nuc_pos": aa_sub["codonNucRange"]["begin"],
                    "end_nuc_pos": aa_sub["codonNucRange"]["end"]
                }

    return variant_summary


def main(options: argparse.Namespace):
    """
    The main function.

    :param options: the options from the command line
    """
    setup_logging(log_level=options.log_level)

    # Open json file to read data
    nextclade_json_data = None
    try:
        with open(options.nextclade_json, encoding="utf-8") as json_file:
            nextclade_json_data = json.load(json_file)
    except FileNotFoundError as e:
        logger.exception(e)
        sys.exit(1)
    except json.decoder.JSONDecodeError as e:
        logger.exception(e)
        sys.exit(1)

    nextclade_results = get_results(json_data=nextclade_json_data)
    nextclade_results.to_csv(
        path_or_buf=f"{options.project_name}_nextclade_results_"
        f"{options.workflow_version}.csv",
        index=False,
    )

    nextclade_variant_summary = get_variant_summary(
        json_data=nextclade_json_data)
    nextclade_variant_summary.to_csv(
        path_or_buf=f"{options.project_name}_nextclade_variant_summary_"
        f"{options.workflow_version}.csv",
        index=False,
    )


if __name__ == "__main__":
    # Get options from command line args
    cli_options = parse_args(args=sys.argv[1:])

    # Inject options and run main
    main(options=cli_options)
