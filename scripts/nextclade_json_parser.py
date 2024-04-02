"""Parses a nextclade json file and outputs nextclade results and variant
summary.
"""

import argparse
import logging
import re
import sys

import pandas as pd

__version__ = "0.4.1"
__author__ = "CDPHE"
__copyright__ = "State of Colorado"
__license__ = "GPL-3.0-or-later"


log = logging.getLogger(__name__)


def parse_args(args: list) -> argparse.Namespace:
    """Parses the command line arguments."""
    parser = argparse.ArgumentParser(
        description="A tool to convert Nextclade JSON to CSV."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"cdphe-nextclade-json-parser {__version__}",
    )
    parser.add_argument(
        "--log_level",
        help="the level to log at",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
        default="INFO",
    )
    parser.add_argument(
        "--nextclade_json",
        help="the nextclade json file path",
    )
    parser.add_argument(
        "--project_name",
        help="the project name",
    )
    parser.add_argument(
        "--workflow_version",
        help="the workflow version",
    )
    return parser.parse_args(args)


def setup_logging(log_level: str) -> None:
    """Sets up the logging."""
    log_format = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(
        level=log_level,
        stream=sys.stdout,
        format=log_format,
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def extract_hsn(sample_name: str) -> str:
    """Extracts the hsn from the sample name.
    The hsn is the first 10 digits of the sample name and begins with 2.
    """
    hsn = ""
    if sample_name:
        match = re.search(r"2[0-9]{9}", sample_name)
        if match:
            hsn = match.group()
    return hsn


def extract_sample_name(fasta_header: str) -> str:
    """Extracts the sample name from the fasta header.
    The fasta header is expected to have the format: >CO-CDPHE-<sample_name>
    """
    sample_name = ""
    if fasta_header:
        match = re.search(r"CO-CDPHE-([0-9a-zA-Z_\\-\\.]+)", fasta_header)
        if match:
            sample_name = match.group(1)
        else:
            sample_name = fasta_header
    return sample_name


def parse_results(json_results: list) -> pd.DataFrame:
    """Parses the results summary from the nextclade json file."""
    # Normalize nextclade json results
    normalized_results = pd.json_normalize(json_results)

    # Rename columns and rearrange columns
    results = normalized_results.copy().rename(
        columns={
            "seqName": "fasta_header",
            "sample_name": "sample_name",
            "clade": "nextclade",
            "totalSubstitutions": "total_nucleotide_mutations",
            "totalDeletions": "total_nucleotide_deletions",
            "totalInsertions": "total_nucleotide_insertions",
            "totalAminoacidSubstitutions": "total_AA_substitutions",
            "totalAminoacidDeletions": "total_AA_deletions",
        }
    )

    # Get hsn and sample_name
    results["hsn"] = results["fasta_header"].apply(extract_hsn)
    results["sample_name"] = results["fasta_header"].apply(extract_sample_name)

    # Reorganize columns and return
    return results[
        [
            "fasta_header",
            "sample_name",
            "hsn",
            "nextclade",
            "total_nucleotide_mutations",
            "total_nucleotide_deletions",
            "total_nucleotide_insertions",
            "total_AA_substitutions",
            "total_AA_deletions",
        ]
    ].sort_values(by=["fasta_header"])


def extract_variant_components(row: pd.Series, mutation_type: str) -> pd.Series:
    """Extracts the values required to create a variant summary."""
    components = pd.Series(
        {
            "fasta_header": None,
            "sample_name": None,
            "hsn": None,
            "variant_name": None,
            "gene": None,
            "codon_position": None,
            "refAA": None,
            "altAA": None,
            "start_nuc_pos": None,
            "end_nuc_pos": None,
        }
    )

    # Extract common components
    components["fasta_header"] = row["seqName"]
    components["sample_name"] = extract_sample_name(row["seqName"])
    components["hsn"] = extract_hsn(row["seqName"])

    # Extract components that need special logic
    if mutation_type == "aa_deletion":
        components["gene"] = row["gene"]
        components["refAA"] = row["refAA"]
        components["altAA"] = "del"
        components["codon_position"] = int(row["codon"]) + 1
        components["start_nuc_pos"] = row["codonNucRange.begin"]
        components["end_nuc_pos"] = row["codonNucRange.end"]
    elif mutation_type == "insertion":
        components["gene"] = ""
        components["refAA"] = "ins"
        components["altAA"] = row["ins"]
        components["codon_position"] = int(row["pos"]) + 1
        components["start_nuc_pos"] = int(row["pos"]) + 1
        components["end_nuc_pos"] = int(row["pos"]) + 1 + len(row["ins"])
    elif mutation_type == "aa_substitution":
        components["gene"] = row["gene"]
        components["refAA"] = row["refAA"]
        if row["queryAA"] == "*":
            components["altAA"] = "stop"
        else:
            components["altAA"] = row["queryAA"]
        components["codon_position"] = int(row["codon"]) + 1
        components["start_nuc_pos"] = row["codonNucRange.begin"]
        components["end_nuc_pos"] = row["codonNucRange.end"]
    else:
        raise ValueError(f"Unknown mutation type: {mutation_type}.")

    # Use extracted components to compose variant name
    components["variant_name"] = (
        f"{components['gene']}_{components['refAA']}"
        f"{components['codon_position']}{components['altAA']}"
    )

    return components


def parse_variant_summary(json_results: list) -> pd.DataFrame:
    """Parses the variant summary from the nextclade json file."""
    # Process amino acid deletions
    aa_deletions = pd.json_normalize(json_results, "aaDeletions", meta=["seqName"])
    aa_deletions_summary = aa_deletions.apply(
        extract_variant_components, mutation_type="aa_deletion", axis=1
    )

    # Process insertions
    insertions = pd.json_normalize(json_results, "insertions", meta=["seqName"])
    insertions_summary = insertions.apply(
        extract_variant_components, mutation_type="insertion", axis=1
    )

    # Process aa substitutions
    aa_substitutions = pd.json_normalize(
        json_results, "aaSubstitutions", meta=["seqName"]
    )
    aa_substitutions_summary = aa_substitutions.apply(
        extract_variant_components, mutation_type="aa_substitution", axis=1
    )

    # Combine all mutations into one summary
    variant_summary = pd.concat(
        [
            df
            for df in [
                aa_deletions_summary,
                insertions_summary,
                aa_substitutions_summary,
            ]
            if not df.empty
        ]
    )

    # Sort and return
    return variant_summary.sort_values(by=["fasta_header", "start_nuc_pos"])


def compose_report_output_path(
    project_name: str, workflow_version: str, report_type: str
) -> str:
    """Composes the output path for the nextclade report file."""
    prefix = ""
    if project_name:
        prefix = f"{project_name}_"

    suffix = ""
    if workflow_version:
        suffix = f"_{workflow_version}"

    path = None
    if report_type == "results":
        path = f"{prefix}nextclade_results{suffix}.csv"
    elif report_type == "variant_summary":
        path = f"{prefix}nextclade_variant_summary{suffix}.csv"
    else:
        raise ValueError(f"Invalid report type: {report_type}")

    return path


def generate_report(
    json_results: list,
    project_name: str,
    workflow_version: str,
    report_type: str,
) -> None:
    """Generates the nextclade report by parsing and writing to file."""
    # Get summary
    report = None
    if report_type == "results":
        report = parse_results(json_results)
    elif report_type == "variant_summary":
        report = parse_variant_summary(json_results)
    log.info("Nextclade JSON %s summary created: %s rows.", report_type, len(report))

    # Get summary output path
    report_path = compose_report_output_path(
        project_name, workflow_version, report_type=report_type
    )

    # Write summary to csv
    report.to_csv(
        path_or_buf=report_path,
        index=False,
    )
    log.info(
        "Nextclade JSON %s report written to csv file: %s", report_type, report_path
    )


def main(options: dict) -> None:
    """The main function."""
    setup_logging(log_level=options["log_level"])
    log.info("Nextclade JSON Parser Start.")

    # Read nextclade json and capture results
    nextclade_json = pd.read_json(options["nextclade_json"], orient="index")
    nextclade_json_results = nextclade_json.loc["results"][0]  # pylint: disable=E1101
    log.info("Nextclade JSON results loaded: %s rows.", len(nextclade_json_results))

    # Take results, transform them into results_summary, and write to file
    generate_report(
        nextclade_json_results,  # type: ignore
        project_name=options["project_name"],
        workflow_version=options["workflow_version"],
        report_type="results",
    )
    log.info("Handling Nextclade JSON results complete.")

    # Take results, transform them into variant_summary, and write to file
    generate_report(
        nextclade_json_results,  # type: ignore
        project_name=options["project_name"],
        workflow_version=options["workflow_version"],
        report_type="variant_summary",
    )
    log.info("Handling Nextclade JSON variant summary complete.")

    log.info("Nextclade JSON Parser End.")


if __name__ == "__main__":
    # Get options from command line args
    cli_args = parse_args(args=sys.argv[1:])
    cli_options = vars(cli_args)

    # Inject options and run main
    main(options=cli_options)
