import logging
import os as _os
import pickle as _pickle
from collections import Counter as _Counter

# 파이프라인 전역 logging 수준. 모든 SKYPE 스크립트가 'from skype_utils import *' 를 가장 먼저
# 실행하므로, 여기서 basicConfig 를 호출해두면 각 스크립트의 자체 basicConfig 보다 먼저 root
# logger 가 설정되어 이 레벨이 전체에 적용된다(자체 basicConfig 는 핸들러가 이미 있어 무시됨).
# 디버그 출력(logging.debug)을 보려면 logging.DEBUG 로 바꾼다.
LOG_LEVEL = logging.INFO
logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=LOG_LEVEL,
    datefmt='%m/%d/%Y %I:%M:%S %p',
    force=True
)

K = 1000
M = 1000 * K

HARD_PATH_COUNT_BASELINE = 500 * K
chrY_MINIMUM_RATIO = 0.7

CHROMOSOME_COUNT = 23

PIPELINE_MODE_KARYOTYPE = "karyotype"
PIPELINE_MODE_VARIANT = "variant"
PIPELINE_MODE_PKL = "pipeline_mode.pkl"
PIPELINE_MODE_NCLOSE_LIMIT = 1000
TYPE4_INDEL_GRAPH_EDGE_PKL = "type4_indel_graph_edges.pkl"

VCF_TYPE4_EVENTS_PKL = "vcf_type4_events.pkl"
VCF_TYPE4_OUTLIER_INDEX_PKL = "vcf_type4_outlier_index.pkl"
VCF_TYPE4_GRAPH_NODE_PREFIX = "vcf_type4_graph_"
VCF_MODE_SUMMARY_JSON = "vcf_mode_summary.json"
VCF_MODE_SUMMARY_TSV = "vcf_mode_summary.tsv"
VCF_TYPE4_MIN_SPAN = 100 * K
INDEL_MERGE_TOLERANCE = 10 * K
RAW_VIRTUAL_INV_MIN_VAF = 0.2


def invert_vcf_strand(strand):
    if strand == "+":
        return "-"
    if strand == "-":
        return "+"
    raise ValueError(f"Unknown strand: {strand}")


def vcf_strands_from_skype_dirs(dir_a, dir_b):
    """
    Compose VCF breakpoint STRANDS from SKYPE path directions.

    SKYPE stores traversal as A_dir => B_dir. VCF-style STRANDS stores
    the joined breakend sides, so B's side is the opposite of traversal.
    """
    if dir_a not in {"+", "-"} or dir_b not in {"+", "-"}:
        raise ValueError(f"Unknown SKYPE path directions: {dir_a}, {dir_b}")
    return dir_a + invert_vcf_strand(dir_b)


def skype_dirs_from_vcf_strands(strands):
    """
    Invert VCF breakpoint STRANDS into SKYPE path directions.

    Examples:
      ++ -> +,-
      +- -> +,+
      -+ -> -,-
      -- -> -,+
    """
    if not isinstance(strands, str) or len(strands) != 2:
        raise ValueError(f"Malformed VCF STRANDS value: {strands}")
    if strands[0] not in {"+", "-"} or strands[1] not in {"+", "-"}:
        raise ValueError(f"Malformed VCF STRANDS value: {strands}")
    return strands[0], invert_vcf_strand(strands[1])


def skype_dirs_from_bnd_alt(alt):
    """
    Infer SKYPE path directions from one VCF BND ALT string.

    This mirrors the ALT forms emitted by 31_depth_analysis.py:
      t[p[ -> +,+
      t]p] -> +,-
      [p[t -> -,+
      ]p]t -> -,-
    """
    if not alt:
        return None
    first_bracket = None
    for i, char in enumerate(alt):
        if char in "[]":
            first_bracket = (i, char)
            break
    if first_bracket is None:
        return None

    bracket_idx, bracket = first_bracket
    dir_a = "+" if bracket_idx > 0 else "-"
    dir_b = "+" if bracket == "[" else "-"
    return dir_a, dir_b


def make_pipeline_mode_config(
    requested_mode=PIPELINE_MODE_KARYOTYPE,
    nclose_node_count=None,
    nclose_limit=PIPELINE_MODE_NCLOSE_LIMIT,
    vcf_input=False,
    vcf_input_path=None,
):
    vcf_input = bool(vcf_input)
    if requested_mode not in {PIPELINE_MODE_KARYOTYPE, PIPELINE_MODE_VARIANT}:
        raise ValueError(f"Unknown pipeline mode: {requested_mode}")

    forced_variant = (
        requested_mode == PIPELINE_MODE_KARYOTYPE
        and nclose_node_count is not None
        and nclose_node_count > nclose_limit
    )
    effective_mode = PIPELINE_MODE_VARIANT if vcf_input or forced_variant else requested_mode

    return {
        "requested_mode": requested_mode,
        "mode": effective_mode,
        "karyotype_mode": effective_mode == PIPELINE_MODE_KARYOTYPE,
        "variant_mode": effective_mode == PIPELINE_MODE_VARIANT,
        "forced_variant": forced_variant,
        "nclose_node_count": nclose_node_count,
        "nclose_limit": nclose_limit,
        "vcf_input": vcf_input,
        "vcf_input_path": vcf_input_path,
    }


def save_pipeline_mode(
    prefix,
    requested_mode=PIPELINE_MODE_KARYOTYPE,
    nclose_node_count=None,
    nclose_limit=PIPELINE_MODE_NCLOSE_LIMIT,
    vcf_input=False,
    vcf_input_path=None,
):
    config = make_pipeline_mode_config(
        requested_mode=requested_mode,
        nclose_node_count=nclose_node_count,
        nclose_limit=nclose_limit,
        vcf_input=vcf_input,
        vcf_input_path=vcf_input_path,
    )
    with open(_os.path.join(prefix, PIPELINE_MODE_PKL), "wb") as f:
        _pickle.dump(config, f)
    return config


def load_pipeline_mode(prefix):
    path = _os.path.join(prefix, PIPELINE_MODE_PKL)
    vcf_summary_path = _os.path.join(prefix, VCF_MODE_SUMMARY_JSON)
    if not _os.path.isfile(path):
        return make_pipeline_mode_config(vcf_input=_os.path.isfile(vcf_summary_path))

    with open(path, "rb") as f:
        config = _pickle.load(f)

    mode = config.get("mode", config.get("requested_mode", PIPELINE_MODE_KARYOTYPE))
    vcf_input = bool(config.get("vcf_input", False)) or _os.path.isfile(vcf_summary_path)
    loaded_config = make_pipeline_mode_config(
        requested_mode=mode,
        nclose_node_count=config.get("nclose_node_count"),
        nclose_limit=config.get("nclose_limit", PIPELINE_MODE_NCLOSE_LIMIT),
        vcf_input=vcf_input,
        vcf_input_path=config.get("vcf_input_path"),
    )
    loaded_config.update({
        "requested_mode": config.get("requested_mode", mode),
        "forced_variant": bool(config.get("forced_variant", False)),
    })
    return loaded_config


def pipeline_mode_is_karyotype(config):
    return config.get("mode") == PIPELINE_MODE_KARYOTYPE


def pipeline_mode_is_variant(config):
    return config.get("mode") == PIPELINE_MODE_VARIANT


def pipeline_mode_is_vcf_input(config):
    return bool(config.get("vcf_input", False))


def describe_pipeline_mode(config):
    text = f"Pipeline mode: {config.get('mode', PIPELINE_MODE_KARYOTYPE)}"
    requested = config.get("requested_mode")
    if requested and requested != config.get("mode"):
        text += f" (requested={requested})"
    if config.get("forced_variant"):
        text += (
            f" forced by nclose_node_count={config.get('nclose_node_count')} "
            f"> {config.get('nclose_limit')}"
        )
    if config.get("vcf_input"):
        text += " (input=vcf)"
    return text


def _type4_indel_event_key(edge):
    type4_tuple = edge.get("type4_tuple")
    if type4_tuple is not None:
        return tuple(type4_tuple)
    return (
        edge.get("contig_name"),
        edge.get("indel_kind"),
        int(edge.get("base_st", edge.get("span_st", 0))),
        int(edge.get("base_nd", edge.get("span_nd", 0))),
    )


def _normalise_type4_indel_event(edge):
    indel_kind = edge.get("indel_kind")
    if indel_kind == "insertion":
        event_type = "i"
    elif indel_kind == "deletion":
        event_type = "d"
    else:
        return None

    chrom = edge.get("base_chrom") or edge.get("chrom")
    if chrom is None:
        return None

    st = int(edge.get("base_st", edge.get("span_st", 0)))
    nd = int(edge.get("base_nd", edge.get("span_nd", 0)))
    lo, hi = sorted((st, nd))
    event_key = _type4_indel_event_key(edge)

    return {
        "event_key": event_key,
        "event_type": event_type,
        "indel_kind": indel_kind,
        "chrom": chrom,
        "st": lo,
        "nd": hi,
        "span_len": hi - lo,
        "contig_name": edge.get("contig_name"),
        "type4_tuple": tuple(edge.get("type4_tuple", ())),
    }


def load_type4_indel_graph_event_index(prefix):
    edge_path = _os.path.join(prefix, TYPE4_INDEL_GRAPH_EDGE_PKL)
    if not _os.path.isfile(edge_path):
        return {}, {}

    with open(edge_path, "rb") as f:
        edges = _pickle.load(f)

    event_by_key = {}
    edge_to_event_key = {}
    for edge in edges:
        event = _normalise_type4_indel_event(edge)
        if event is None:
            continue
        event_key = event["event_key"]
        event_by_key.setdefault(event_key, event)
        edge_to_event_key[(tuple(edge["src"]), tuple(edge["dst"]))] = event_key

    return event_by_key, edge_to_event_key


def count_type4_indel_graph_edges_in_path(index_path, edge_to_event_key):
    usage = _Counter()
    if not edge_to_event_key:
        return usage

    compact_path = [tuple(list(node)[:2]) for node in index_path]
    for left, right in zip(compact_path, compact_path[1:]):
        event_key = edge_to_event_key.get((left, right))
        if event_key is not None:
            usage[event_key] += 1
    return usage


def summarize_type4_indel_graph_usage(prefix, path_records, weights, import_index_path_func):
    event_by_key, edge_to_event_key = load_type4_indel_graph_event_index(prefix)
    event_weight = _Counter()
    path_event_usage = {}

    if not event_by_key:
        return event_by_key, event_weight, path_event_usage

    for path_idx, record in enumerate(path_records):
        if path_idx >= len(weights):
            break
        paf_loc = record[0]
        index_path = import_index_path_func(paf_loc)
        usage = count_type4_indel_graph_edges_in_path(index_path, edge_to_event_key)
        if not usage:
            continue
        path_event_usage[paf_loc] = usage
        path_weight = float(weights[path_idx])
        for event_key, count in usage.items():
            event_weight[event_key] += path_weight * count

    return event_by_key, event_weight, path_event_usage


def add_weighted_indel_event(event_by_key, event_type, chrom, st, nd, weight,
                             source=None, tolerance=INDEL_MERGE_TOLERANCE):
    if event_type not in {"i", "d"}:
        raise ValueError(f"Unknown indel event type: {event_type}")

    lo, hi = sorted((int(st), int(nd)))
    for key, event in event_by_key.items():
        if (
            event["event_type"] == event_type
            and event["chrom"] == chrom
            and abs(event["st"] - lo) <= tolerance
            and abs(event["nd"] - hi) <= tolerance
        ):
            event["weight"] += float(weight)
            if source is not None:
                event.setdefault("sources", []).append(str(source))
            return key

    key = (event_type, chrom, lo, hi)
    event_by_key[key] = {
        "event_type": event_type,
        "chrom": chrom,
        "st": lo,
        "nd": hi,
        "weight": float(weight),
        "sources": [str(source)] if source is not None else [],
    }
    return key


def indel_event_source_label(event, fallback="INDEL"):
    unique_sources = []
    for source in event.get("sources", []):
        if source not in unique_sources:
            unique_sources.append(source)
    if not unique_sources:
        return fallback
    if len(unique_sources) <= 3:
        return "|".join(unique_sources)
    return "|".join(unique_sources[:3] + [f"+{len(unique_sources) - 3}"])
