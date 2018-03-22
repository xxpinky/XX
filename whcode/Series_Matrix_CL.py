import os
import re
import sys
import logging

import pandas as pd
import numpy as np

from config import Config


array_logger = logging.getLogger("ArrayLogger")
array_logger.setLevel(logging.DEBUG)
GPL = "E:\\Data\\xx_data\\CHIP-NORMAL"


def read_gpl(gpl="GPL570", probe=0, gene=10):
    """
    gpl, probe, gene: platform id, probe and gene col number
    common references: 570-10, 6864-6, 6884-12, 5175-9, 17586-7(SPLIT[1])
    10295-5, 96-10, 4133-9, 10558-5, 14550-6, 6244-9(SPLIT[1])
    """

    probe_gene = {}
    with open(os.path.join(GPL, "{}-tbl-1.txt".format(gpl))) as f:
        for line in f:
            items = line.split(sep="\t")
            p = items[probe].strip()
            g = items[gene].split(sep="///")[0].strip()

            if g:
                probe_gene[p] = g
            else:
                probe_gene[p] = None
    return probe_gene


def extract_features(fname, feature):

    with open(fname, encoding="utf-8") as f:
        for line in f:
            if line.startswith("!Sample_{}".format(feature)):
                return [i.replace('"', "") for i in
                        line.strip().split(sep="\t")[1:]]

    array_logger.critical("No such feature")
    sys.exit()
    return


def clean_series_matrix(fname='GSE13670_GPL570.txt',
                        feature="characteristics_ch1"):
    """缺失值以基因在样本中中值填充，基因表达值由探针中值得到,
    返回：表达谱的dataframe， 如果表达大于100则取log2，
    注意基因名未过滤，可能包含---等"""

    name, gpl = fname[:-4].split(sep="_")
    gpl_config = getattr(Config, gpl)
    p2g = read_gpl(gpl, gpl_config["probe"], gpl_config["gene"])

    data = pd.read_csv(fname, header=0, comment="!", sep="\t", index_col=0)
    data = data.T.fillna(data.median(axis=1)).T
    data.columns = extract_features(fname, feature)
    data["Symbol"] = data.index.to_series().map(p2g)
    data.set_index("Symbol", inplace=True)
    data = data.loc[data.index.dropna().unique()]
    data = data.groupby(level=0).median()

    data_max = data.abs().max().max()
    if data_max > 100:
        data = data.apply(np.log2)

    data.name = name
    data.to_csv("{}_cleaned.csv".format(name))
    return data


def divide_series_matrix(matrix, case_pattern, control_pattern, merge=True):
    """Divive matrix into case control or only return fold change
    Pattern is python regex pattern"""

    case = []
    control = []

    for i in matrix.columns:

        match1 = re.search(case_pattern, i)
        if match1:
            case.append(i)
        else:
            match2 = re.search(control_pattern, i)
            if match2:
                control.append(i)

    print("Matched case: {}".format(case))
    print("Matched contro: {}".format(control))
    case = matrix[case]
    control = matrix[control]

    if merge:
        fd = case.median(axis=1) - control.median(axis=1)
        fd.to_csv("{}_fd.csv".format(matrix.name))
    else:
        control.to_csv("{}_control.csv".format(matrix.name))
        case.to_csv("{}_case.csv".format(matrix.name))


if __name__ == "__main__":

    fname, feature, p1, p2, merge = sys.argv[1:]
    # fname, feature, p1, p2, merge = "data/GSE13670_GPL570.txt", "title", "treated", "control" "1 combine-foldchange/0 ngf sample"
    data = clean_series_matrix(fname, feature)
    divide_series_matrix(data, r"{}".format(p1), r"{}".format(p2), int(merge))
