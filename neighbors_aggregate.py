#!/usr/bin/env python3

"""
Created on 01 Sep. 2025
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.0.0"

import argparse
import logging
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import numpy
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import seaborn as sns
from statannotations.Annotator import Annotator


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level og the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """

    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[level]

    if os.path.exists(path):
        os.remove(path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def get_domains(domains_file_path):
    """
    Extract the domains in order as they are set in the CSV file.
    
    :param domains_file_path: the path to the protein domains CSV file.
    :type domains_file_path: str
    :return: the ordered domains.
    :rtype: list
    """
    domains = None
    if domains_file_path:
        df = pd.read_csv(domains_file_path, sep=",")
        domains = list(df["domain"])
    return domains



def extract_colors(path, grouped):
    """
    Extract the colors by conditions for the boxplots and the dots.

    :param path: the path to the CSV file.
    :type path: str
    :param grouped: the grouped conditions.
    :type grouped: list
    :return: the colors to use.
    :rtype: dict
    """
    colors_by_condition = {"boxplots": {}, "dots": {}}
    first_item_group_not_stored = True
    df = pd.read_csv(path, sep=",", header=0)
    for _, row in df.iterrows():
        if grouped and row["condition"] in grouped:
            if first_item_group_not_stored:
                colors_by_condition["boxplots"]["/".join(grouped)] = row["boxplot color"]
                colors_by_condition["dots"]["/".join(grouped)] = row["dot color"]
                first_item_group_not_stored = False
        else:
            colors_by_condition["boxplots"][row["condition"]] = row["boxplot color"]
            colors_by_condition["dots"][row["condition"]] = row["dot color"]
    return colors_by_condition


def get_all_domains_with_contacts(neighbors_dict):
    """
    Get all the domains where a neighbor contact is present.

    :param neighbors_dict: the neighbors contact.
    :type neighbors_dict: dict
    :return: the domains where a neighbor contact is present.
    :rtype: list
    """
    domains_with_neighbors_contact = set()
    for condition in neighbors_dict:
        for smp in neighbors_dict[condition]:
            for domain in neighbors_dict[condition][smp]:
                domains_with_neighbors_contact.add(domain)
    return domains_with_neighbors_contact


def aggregate_neighbors(conditions, roi, md_time, dir_path, grouped):
    """
    Extract the number of contacts by region for each sample in a condition.

    :param conditions: the conditions dataframe.
    :type conditions: pandas.DataFrame
    :param roi: the region of interest.
    :type roi: str
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param grouped: the grouped conditions.
    :type grouped: list
    :return: the aggregated data for each frame and the conditions (in case ony condition is removed) and the domain of
    interest.
    :rtype: pandas.DataFrame, str
    """
    pattern_sample = re.compile(f"neighborhood_(.+)_{roi}.csv")
    data = {}
    conditions_to_remove = []
    for _, row_condition in conditions.iterrows():
        by_condition = []
        try:
            for fn in os.listdir(row_condition["path"]):
                if fn.startswith("neighborhood") and fn.endswith(".csv"):
                    by_condition.append(fn)
        except FileNotFoundError as exc:
            logging.error(exc, exc_info=True)
            sys.exit(1)
        if len(by_condition) == 0:
            conditions_to_remove.append(row_condition["condition"])
            logging.warning(f"Condition {row_condition['condition']}: no RMSD files, this condition is skipped.")
            continue
        logging.info(f"Aggregating {len(by_condition)} file{'s' if len(by_condition) > 1 else ''} data for condition: "
                     f"{row_condition['condition']}")

        # check if the condition belongs to the grouped conditions
        if grouped and row_condition.iloc[0] in grouped:
            condition = "/".join(grouped)
        else:
            condition = row_condition.iloc[0]
        if condition not in data:
            data[condition] = {}

        for item in sorted(by_condition):
            logging.info(f"\t\t- {item}")
            match_sample = pattern_sample.match(item)
            if match_sample:
                sample = match_sample.group(1)
            else:
                logging.error(f"No sample found with the pattern \"{pattern_sample.pattern}\" in the file {item}")
                sys.exit(1)
            data[condition][sample] = {}
            df_current = pd.read_csv(os.path.join(row_condition["path"], item), sep=",")

            # check there is only one domain for the residue 1 positions
            residue1_domains = df_current["residue 1 domain"].unique()
            if  len(residue1_domains) > 1:
                logging.error(f"For {sample}: more than one domain in the columns 'ROI partner domain' "
                              f"({', '.join(residue1_domains)}) of the neighbors contact CSV file.")
                sys.exit(1)

            # get the atom and the residue pairs contacts by domain
            for residue2_domain in df_current["residue 2 domain"].unique():
                df_res2_dom = df_current[df_current["residue 2 domain"] == residue2_domain]
                pairs_neighbors = 0
                for unique_residue1_position in df_res2_dom["residue 1 position"].unique():
                    df_by_pos1_in_res2_dom = df_res2_dom[df_res2_dom["residue 1 position"] == unique_residue1_position]
                    pairs_neighbors += df_by_pos1_in_res2_dom['residue 2 position'].nunique()
                data[condition][sample][residue2_domain] = {"by atom": len(df_res2_dom),
                                                            "by residue": pairs_neighbors}

    # complete missing data in some domains
    expected_domains = get_all_domains_with_contacts(data)
    for condition in data:
        for smp in data[condition]:
            for expected_domain in expected_domains:
                if expected_domain not in data[condition][smp]:
                    data[condition][smp][expected_domain] = {"by atom": 0, "by residue": 0}

    # reorganize the data
    reorganized_dict = {"sample": [], "conditions": [], "domains": [], "by atom": [], "by residue": []}
    for condition in data:
        for smp in data[condition]:
            for domain in data[condition][smp]:
                reorganized_dict["sample"].append(smp)
                reorganized_dict["conditions"].append(condition)
                reorganized_dict["domains"].append(domain)
                reorganized_dict["by atom"].append(data[condition][smp][domain]["by atom"])
                reorganized_dict["by residue"].append(data[condition][smp][domain]["by residue"])

    df_out = pd.DataFrame.from_dict(reorganized_dict)
    out_path = os.path.join(dir_path, f"neighbors_aggregated_{roi.lower().replace(' ', '-')}_{md_time}-ns.csv")
    df_out.to_csv(out_path, index=False)
    logging.info(f"Aggregated CSV file saved: {os.path.abspath(out_path)}")

    return df_out


def update_domains_order(labels, domains_ordered):
    """
    Update and order the domains by adding before, between and after annotations if some contacts are present outside
    the domains.

    :param labels: the labels.
    :type labels: list
    :param domains_ordered: the ordered domains on the protein.
    :type domains_ordered: list
    :return: the X axis labels ordered.
    :rtype: list
    """
    updated_labels = []
    # get the annotation before any domains
    for i in range(len(labels)):
        if labels[i].startswith("before"):
            updated_labels.append(labels[i])
            break
    labels[:] = [x for x in labels if not x.startswith("before")]
    logging.debug("Updating the labels:")
    logging.debug(f"\tinitial labels: {labels}")
    logging.debug(f"\tbefore domains:")
    logging.debug(f"\t\tnew labels:\t\t{updated_labels}")
    logging.debug(f"\t\tinitial labels:\t{labels}")
    # get the domains as in the ordered domains and add the between domains
    for dom in domains_ordered:
        domain_index_in_labels = {}
        logging.debug(f"\tdomain {dom}:")
        for i in range(len(labels)):
            if dom == labels[i]:
                domain_index_in_labels["dom"] = i
            elif labels[i].startswith(f"between {dom}"):
                domain_index_in_labels["between"] = i
        if "dom" in domain_index_in_labels:
            updated_labels.append(labels[domain_index_in_labels["dom"]])
        if "between" in domain_index_in_labels:
            updated_labels.append(labels[domain_index_in_labels["between"]])
        if "dom" in domain_index_in_labels:
            labels[:] = [x for x in labels if not x == dom]
        if "between" in domain_index_in_labels:
            labels[:] = [x for x in labels if not x.startswith(f"between {dom}")]
        logging.debug(f"\t\tnew labels:\t\t{updated_labels}")
        logging.debug(f"\t\tinitial labels:\t{labels}")
    # get the annotation after all the domains
    for i in range(len(labels)):
        if labels[i].startswith("after"):
            updated_labels.append(labels[i])
            break
    labels[:] = [x for x in labels if not x.startswith("after")]
    logging.debug(f"\tafter domains:")
    logging.debug(f"\t\tnew labels:\t\t{updated_labels}")
    logging.debug(f"\t\tinitial labels:\t{labels}")

    return updated_labels


def compute_stats(src, domains, out_dir, roi, md_time, level_of_interaction):
    """
    Test the different domains contacts with the region of interest between the different conditions.
    A Mann-Whitney U test is performed with the null hypothesis is that the condition 1 group is greater than the
    condition 2 group.

    :param src: the contacts dataframe.
    :type src: pandas.DataFrame
    :param domains: the domains.
    :type domains: list
    :param out_dir: the output directory path.
    :type out_dir: str
    :param roi: the region of interest.
    :type roi: str
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param level_of_interaction: the level of interaction 'by atom' or 'by residue'
    :type level_of_interaction: str
    """
    level_of_interaction_txt = level_of_interaction.split(" ")[1]
    logging.info(f"At the {level_of_interaction_txt} level:")
    logging.info("\tComputing Mann-Whitney U test with a null hypothesis group 1 is greater than group 2:")
    data = {"contact with": [], "group 1": [], "group 2": [], "p-value": [], "statistic": [], "test":[], "H0": [],
            "comment": []}
    # get the conditions as a list, for loop performed to keep the conditions' order
    conditions = []
    for condition in src["conditions"]:
        if condition not in conditions:
            conditions.append(condition)
    # extract the data and compute the statistic test
    for domain in domains:
        domain_rows = src[src["domains"] == domain]
        for i in range(0, len(conditions) - 1):
            for j in range(i + 1, len(conditions)):
                data["contact with"].append(domain)
                data["test"].append("Mann-Whitney U")
                data["group 1"].append(conditions[i])
                data["group 2"].append(conditions[j])
                try:
                    test = mannwhitneyu(x=domain_rows[level_of_interaction][domain_rows["conditions"] == conditions[i]],
                                        y=domain_rows[level_of_interaction][domain_rows["conditions"] == conditions[j]],
                                        alternative="greater")
                    data["p-value"].append(test.pvalue)
                    data["statistic"].append(test.statistic)
                    data["comment"].append("")
                except ValueError:
                    txt = "All the numbers are identical in the Mann-Whitney U test"
                    data["p-value"].append("N/A")
                    data["statistic"].append("N/A")
                    data["comment"].append(txt)
                    logging.warning(f"\t{txt} for the domain {domain} between {conditions[i]} and {conditions[j]}. "
                                    f"The test output is set to N/A.")
                data["H0"].append(f"{conditions[i]} is greater than {conditions[j]}")
    out_path = os.path.join(out_dir,
                            f"{interaction_level.replace(' ', '-')}_statistics_"
                            f"{roi.lower().replace(' ', '-')}_{md_time}-ns.csv")
    pd.DataFrame.from_dict(data).to_csv(out_path, index=False)
    logging.info(f"\t\t{level_of_interaction_txt.capitalize()} level neighbors statistics file saved: {out_path}")


def all_values_equals_correction(df, pairs_list, level):
    """
    If the values between two conditions are the sames, the Mann-Whitney test cannot be performed, a small variation is
    added to the first value of the first condition.

    :param df: the contacts dataframe.
    :type df: pandas.DataFrame
    :param pairs_list: the list of two tuples (domains and condition) to test with the Mann-Whitney test.
    :type pairs_list: list
    :param level: the level of interaction 'by atom' or 'by residue'
    :type level: str
    :return: the updated contacts dataframe.
    :rtype: pandas.DataFrame
    """
    for pairs in pairs_list:
        pair_1 = set(df[level][(df["domains"] == pairs[0][0]) & (df["conditions"] == pairs[0][1])])
        pair_2 = set(df[level][(df["domains"] == pairs[1][0]) & (df["conditions"] == pairs[1][1])])
        if len(pair_1) == 1 and len(pair_1) == len(pair_2):
            # add 0.00000001 to the first value of the first condition
            # to be able to perform the Mann-Withney test if all the values are the same, see:
            # https://stackoverflow.com/questions/54212583/change-1st-row-of-a-dataframe-based-on-a-condition-in-pandas
            # for explanations
            mask = (df["domains"] == pairs[0][0]) & (df["conditions"] == pairs[0][1])
            idx = mask.idxmax() if mask.any() else numpy.repeat(False, len(df))
            previous_value = df.loc[idx, level]
            df.loc[idx, level] = previous_value + 0.00000001
            logging.warning(f"\t\tDomain \"{pairs[0][0]}\" conditions \"{pairs[0][1]}\" and \"{pairs[1][1]}\" have the "
                            f"same values {previous_value}. The first value of the condition \"{pairs[0][1]}\" is set "
                            f"to {df.loc[idx, level]} to perform the Mann-Withney test.")
    return df


def boxplot_aggregated(src, roi, colors_plot, md_time, dir_path, fmt, domains, subtitle_arg, level_of_interaction):
    """
    Create a boxplot by conditions.

    :param src: the data source.
    :type src: pandas.DataFrame
    :param roi: region of interest, the region in contact with the other domains.
    :type roi: str
    :param colors_plot: the colors to use.
    :type colors_plot: dict
    :param md_time: the molecular dynamics duration.
    :type md_time: int
    :param dir_path: the output directory path.
    :type dir_path: str
    :param fmt: the plot output format.
    :type fmt: str
    :param domains: the updated and ordered list of domains.
    :type domains: list
    :param subtitle_arg: the subtitle of the plot.
    :type subtitle_arg: str
    :param level_of_interaction: the level of interaction 'by atom' or 'by residue'
    :type level_of_interaction: str
    """
    logging.info(f"\tPlotting the aggregated neighborhood contacts at the {level_of_interaction.split(' ')[1]}s level "
                 f"by condition:")
    plt.figure(figsize=(15, 15))
    # create the statistical pairs annotations
    boxplot_pairs = []
    conditions = list(set(src["conditions"]))
    for domain in domains:
        for i in range(0, len(conditions) - 1):
            for j in range(i + 1, len(conditions)):
                boxplot_pairs.append(((domain, conditions[i]), (domain, conditions[j])))
    # for a domain, if all the values for both tested conditions are equals, add a little variation to the first value
    # of the first condition to be able to perform a Mann-Whitney test
    src = all_values_equals_correction(src, boxplot_pairs, level_of_interaction)
    # creating the plotting parameters
    plotting_parameters = {
       "data": src,
        "x": "domains",
        "y": level_of_interaction,
        "hue": "conditions",
        "order": domains
    }

    # create the plot
    with sns.plotting_context():
        ax = sns.boxplot(**plotting_parameters, palette=colors_plot["boxplots"])
        sns.stripplot(**plotting_parameters, size=8, marker="o", linewidth=2, dodge=True, palette=colors_plot["dots"])

        # annotate with the statistical test
        annotator = Annotator(ax, boxplot_pairs, **plotting_parameters)
        annotator.configure(test="Mann-Whitney", text_format="star", hide_non_significant=True)
        annotator.apply_test(alternative="greater")
        annotator.annotate()

        # add separators between conditions
        [ax.axvline(x + 0.5, alpha=0.2) for x in ax.get_xticks()]

        # modify the ticks labels for the X axis by adding new lines every 3 words
        modified_x_labels = [re.sub(r'(\w+ \w+ \w+)( )',
                                    r'\1\n', x_label.get_text()) for x_label in ax.get_xticklabels()]
        # set the number of ticks for the X axis to avoid a matplotlib warning
        ax.set_xticks([num_tick for num_tick in range(len(modified_x_labels))])
        ax.set_xticklabels(modified_x_labels, rotation=45, horizontalalignment="right")

        # remove extra legend handles and add the count of samples by condition
        handles, labels = ax.get_legend_handles_labels()
        custom_labels = []
        number_of_labels = len(colors_plot["boxplots"])
        for label in labels[:number_of_labels]:
            sample_set = set(src[src["conditions"] == label]["sample"])
            custom_labels.append(f"{label} ({len(sample_set)})")
        ax.legend(handles[:3], custom_labels, title="Condition")

        plt.suptitle(f"Neighbors {level_of_interaction.split(' ')[1]}s contacts by domain with the {roi} at {md_time} "
                     f"ns of molecular dynamics", fontsize="large",
                     fontweight="bold")
        subtitle = "Mann-Withney H0: first condition greater than the second."
        if subtitle_arg:
            subtitle = f"{subtitle_arg}, {subtitle}"
        plt.title(subtitle)
        plt.xlabel("Domains", fontweight="bold")
        plt.ylabel(f"Number of contacts", fontweight="bold")
        plot = ax.get_figure()
        out_path_plot = os.path.join(dir_path, f"{level_of_interaction.replace(' ', '-')}_neighbors_"
                                               f"aggregated_{roi.lower().replace(' ', '-')}_{md_time}-ns."
                                               f"{fmt}")
        plot.savefig(out_path_plot)
    logging.info(f"\t\tBoxplot at {level_of_interaction.split(' ')[1]}s level aggregated neighbors by condition: "
                 f"{os.path.abspath(out_path_plot)}")


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.

    Aggregate the neighbors contacts by domains in one plot to compare between various conditions.

    The input is a comma separated file without header which first column is the condition, the second column the path 
    of the directory containing the contacts analysis files and the third column the color in hexadecimal format. i.e:

    insertions,tests/inputs/insertions,#fc030b
    WT,tests/inputs/WT,#0303fc

    The output boxplots of neighborhood contacts (at both the atomic and residue levels) for each condition. In 
    addition, a file containing the results of the statistical tests is produced.
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="the path to the output directory.")
    parser.add_argument("-t", "--md-time", required=True, type=int,
                        help="the molecular dynamics duration in nanoseconds.")
    parser.add_argument("-r", "--region-of-interest", required=True, type=str,
                        help="the name of the studied region of interest. The plot neighbors results file last part, "
                             "in example for the testing dataset files (data/plot_neighbors_outputs), the \"HVR\" part"
                             "of the file name.")
    parser.add_argument("-d", "--domains", required=True, type=str,
                        help="a sample CSV domains annotation file, to set the order of the protein domains on the X "
                             "axis. If this option is not used, the domains will be displayed randomly.")
    parser.add_argument("-g", "--group", required=False, nargs="+", type=str,
                        help="a list of conditions, separated by spaces, to group as they appear in the first column "
                             "of the input file. The color used will be the color of the first condition.")
    parser.add_argument("-s", "--subtitle", required=False, type=str,
                        help="Free text used as a subtitle for the boxplots.")
    parser.add_argument("-x", "--format", required=False, default="svg",
                        choices=["eps", "jpg", "jpeg", "pdf", "pgf", "png", "ps", "raw", "svg", "svgz", "tif", "tiff"],
                        help="the output plots format: 'eps': 'Encapsulated Postscript', "
                             "'jpg': 'Joint Photographic Experts Group', 'jpeg': 'Joint Photographic Experts Group', "
                             "'pdf': 'Portable Document Format', 'pgf': 'PGF code for LaTeX', "
                             "'png': 'Portable Network Graphics', 'ps': 'Postscript', 'raw': 'Raw RGBA bitmap', "
                             "'rgba': 'Raw RGBA bitmap', 'svg': 'Scalable Vector Graphics', "
                             "'svgz': 'Scalable Vector Graphics', 'tif': 'Tagged Image File Format', "
                             "'tiff': 'Tagged Image File Format'. Default is 'svg'.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help="the path for the log file. If this option is skipped, the log file is created in the "
                             "output directory.")
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If the option is skipped, log level is INFO.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("input", type=str,
                        help="the path to the CSV (comma separated without header) file which first column is the "
                             "condition, the second column the path of the directory containing the plots_contacts "
                             "script CSV output files and the third column the color.")
    args = parser.parse_args()

    # create output directory if necessary
    os.makedirs(args.out, exist_ok=True)
    # create the logger
    if args.log:
        log_path = args.log
    else:
        log_path = os.path.join(args.out, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
    create_log(log_path, args.log_level)

    logging.info(f"version: {__version__}")
    logging.info(f"CMD: {' '.join(sys.argv)}")
    logging.info(f"MD simulation time: {args.md_time} ns")

    ordered_domains = get_domains(args.domains)
    data_conditions = pd.read_csv(args.input, sep=",", header=0)
    colors = extract_colors(args.input, args.group)
    df_contacts = aggregate_neighbors(data_conditions, args.region_of_interest, args.md_time, args.out, args.group)
    updated_ordered_domains = update_domains_order(list(set(df_contacts["domains"])), ordered_domains)

    for interaction_level in ["by atom", "by residue"]:
        compute_stats(df_contacts, updated_ordered_domains, args.out, args.region_of_interest, args.md_time,
                      interaction_level)
        boxplot_aggregated(df_contacts, args.region_of_interest, colors, args.md_time, args.out, args.format,
                           updated_ordered_domains, args.subtitle, interaction_level)
