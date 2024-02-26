import networkx as nx
from tqdm import tqdm
from pathlib import Path
import json
from multiprocessing import Pool

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


def get_edges(args):
    score_file, score_name, threshold = args
    edges = []
    with open(score_file) as f:
        for i, line in enumerate(f):
            if not len(line):
                continue
            if i == 0:
                columns = line.strip().split("\t")
                continue
            parts = dict(zip(columns, line.strip().split("\t")))
            if parts[score_name] == "nan":
                continue
            score = float(parts[score_name])
            if score < threshold:
                continue
            edges.append((parts["query_pocket"], parts["target_pocket"], score))
    return edges

def label_protein_pocket_clusters(score_dir, score_name, thresholds, strong=True, num_threads=8):
    thresholds = sorted(thresholds)
    if strong:
        graph = nx.DiGraph()
    else:
        graph = nx.Graph()
    with Pool(num_threads) as p:
        for edge_list in tqdm(p.imap(get_edges, [(score_file, score_name, thresholds[0]) for score_file in score_dir.iterdir()]), total=1060):
            graph.add_weighted_edges_from(edge_list, weight="score")
    key_to_component = {thresholds[0]: get_connected_components(graph, strong=strong)}
    for threshold in tqdm(thresholds[1:]):
        graph.remove_edges_from([e for e in graph.edges(data=True) if e[2]["score"] < threshold])
        key_to_component[threshold] = get_connected_components(graph, strong=strong)
    return key_to_component

def get_present(args):
    score_file, score_name = args
    present = set()
    with open(score_file) as f:
        for i, line in enumerate(f):
            if not len(line):
                continue
            if i == 0:
                columns = line.strip().split("\t")
                continue
            parts = dict(zip(columns, line.strip().split("\t")))
            if parts[score_name] == "nan":
                continue
            present.add(parts["query_pocket"])
    return present

def label_present(score_dir, score_name, num_threads=8):
    present = set()
    with Pool(num_threads) as p:
        for present_set in tqdm(p.imap(get_present, [(score_file, score_name) for score_file in score_dir.iterdir()]), total=1060):
            present |= present_set
    return present

def label_df(pocket_df, cluster_dir, num_threads=8):
    for cluster_file in tqdm(Path(cluster_dir).iterdir()):
        with open(cluster_file) as f:
            key_to_component = json.load(f)
        score_name, strong = cluster_file.stem.split("__")
        present = label_present(cluster_dir, score_name, num_threads=num_threads)
        pocket_df[f"{score_name}__no_hit"] = pocket_df["system_ID"].apply(lambda x: x not in present)
        for key in key_to_component:
            col = f"{score_name}__{key}__{strong}__component"
            pocket_df[col] = pocket_df["system_ID"].apply(lambda x: key_to_component[key].get(x, None))
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
    parser.add_argument("--thresholds", type=float, nargs="+", default=[0.3, 0.5, 0.7, 0.99])
    parser.add_argument("--num_threads", type=int, default=8)
    args = parser.parse_args()
    key_to_components = label_protein_pocket_clusters(Path(args.score_dir), args.score_name, args.thresholds, strong=args.type == "strong", num_threads=args.num_threads)
    with open(args.output_dir / f"{args.score_name}__{args.type}.json", "w") as f:
        json.dump(key_to_components, f)

if __name__ == "__main__":
    main()