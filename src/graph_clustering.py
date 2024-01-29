import networkx as nx
from tqdm import tqdm
from pathlib import Path
import extract_scores
import numpy as np
import json

def get_connected_components(graph, strong=True):
    if strong:
        components = sorted(nx.strongly_connected_components(graph), key=len, reverse=True)
    else:
        components = sorted(nx.connected_components(graph), key=len, reverse=True)
    key_to_component = {}
    for i, c in enumerate(components):
        for k in c:
            key_to_component[k] = i
    return key_to_component


def label_protein_pocket_clusters(score_files, score_name, thresholds, strong=True):
    thresholds = sorted(thresholds)
    if strong:
        graph = nx.DiGraph()
    else:
        graph = nx.Graph()
    columns = extract_scores.INFO_COLUMNS
    for suffix in ["_weighted_sum", "_weighted_max", "_max"]:
        for s in extract_scores.SCORE_NAMES:
            columns.append(f"{s}{suffix}")
            if suffix != "_weighted_sum":
                columns.append(f"{s}{suffix}_mapping")
    for score_file in tqdm(score_files):
        with open(score_file) as f:
            for i, line in enumerate(f):
                if i == 0 or not len(line):
                    continue
                parts = line.strip().split("\t")
                parts = dict(zip(columns, parts))
                score = float(parts[score_name])
                if score < thresholds[0] or np.isnan(score) or str(score) == "nan":
                    continue
                graph.add_edge(parts["query_pocket"], parts["target_pocket"], score=score)
    key_to_component = {thresholds[0]: get_connected_components(graph, strong=strong)}
    for threshold in tqdm(thresholds[1:]):
        graph.remove_edges_from([e for e in graph.edges(data=True) if e[2]["score"] < threshold])
        key_to_component[threshold] = get_connected_components(graph, strong=strong)
    return key_to_component

def label_df(pocket_df, cluster_dir):
    for cluster_file in tqdm(Path(cluster_dir).iterdir()):
        with open(cluster_file) as f:
            key_to_component = json.load(f)
        score_name, strong = cluster_file.stem.split("__")
        for key in key_to_component:
            col = f"{score_name}__{key}__{strong}__component"
            pocket_df[col] = pocket_df["pocket_ID"].apply(lambda x: key_to_component[key].get(x, None))
            max_component = int(pocket_df[col].max()) + 1
            pocket_df.loc[pocket_df[col].isnull(), col] = [int(max_component) + i for i in range(len(pocket_df[pocket_df[col].isnull()]))]
    return pocket_df

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--score_name", type=str, required=True)
    parser.add_argument("--score_dir", type=Path, required=True)
    parser.add_argument("--type", type=str, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--thresholds", type=float, nargs="+", default=[0.5, 0.7, 0.99])
    args = parser.parse_args()
    score_dir = Path(args.score_dir)
    score_files = list(score_dir.iterdir())
    if args.type == "strong":
        key_to_components = label_protein_pocket_clusters(score_files, args.score_name, args.thresholds)
        with open(args.output_dir / f"{args.score_name}__strong.json", "w") as f:
            json.dump(key_to_components, f)
    else:
        key_to_components = label_protein_pocket_clusters(score_files, args.score_name, args.thresholds, strong=False)
        with open(args.output_dir / f"{args.score_name}__weak.json", "w") as f:
            json.dump(key_to_components, f)

if __name__ == "__main__":
    main()