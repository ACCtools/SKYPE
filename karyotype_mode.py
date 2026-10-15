"""Karyotype/variant output mode shared by stages 25, 30 and 31.

Karyotype mode is the default. It adds the greedy-filtered (``_filter``) and
karyotype-clustered (``_cluster``) fits next to the unsuffixed raw fit. Variant
mode, VCF input and graphs above the NClose node limit keep only the raw fit.
"""
from __future__ import annotations

import json
import logging
import os

KARYOTYPE_MODE_JSON = "karyotype_mode.json"
NCLOSE_NODE_LIMIT = 1000
STAGE_SUFFIXES = ("_filter", "_cluster")

MODE_KARYOTYPE = "karyotype"
MODE_VARIANT = "variant"


def decide_karyotype_mode(requested_mode, nclose_node_count, vcf_input,
                          nclose_limit=NCLOSE_NODE_LIMIT):
    if requested_mode not in {MODE_KARYOTYPE, MODE_VARIANT}:
        raise ValueError(f"Unknown output mode: {requested_mode}")
    forced_variant = (
        requested_mode == MODE_KARYOTYPE
        and (bool(vcf_input) or int(nclose_node_count) > nclose_limit)
    )
    mode = MODE_VARIANT if forced_variant else requested_mode
    return {
        "requested_mode": requested_mode,
        "mode": mode,
        "karyotype_mode": mode == MODE_KARYOTYPE,
        "forced_variant": forced_variant,
        "vcf_input": bool(vcf_input),
        "nclose_node_count": int(nclose_node_count),
        "nclose_limit": int(nclose_limit),
    }


def save_karyotype_mode(prefix, config):
    path = os.path.join(prefix, KARYOTYPE_MODE_JSON)
    with open(path, "w") as handle:
        json.dump(config, handle, indent=2)
        handle.write("\n")
    return path


def load_karyotype_mode(prefix):
    """Return the stage-25 decision; older results without it are variant-only."""
    path = os.path.join(prefix, KARYOTYPE_MODE_JSON)
    if not os.path.isfile(path):
        logging.warning(
            "%s is missing; writing only unsuffixed outputs. "
            "Rerun SKYPE from stage 22 to restore karyotype outputs.", path
        )
        return {"mode": MODE_VARIANT, "karyotype_mode": False}
    with open(path) as handle:
        return json.load(handle)


def karyotype_mode_enabled(prefix):
    return bool(load_karyotype_mode(prefix).get("karyotype_mode", False))


def describe_karyotype_mode(config):
    text = f"Output mode: {config['mode']}"
    if config.get("forced_variant"):
        reason = (
            "VCF input" if config.get("vcf_input")
            else f"{config['nclose_node_count']} NClose nodes > {config['nclose_limit']}"
        )
        text += f" (karyotype requested; forced by {reason})"
    return text


def suffixed_output_names(suffix):
    """Every file written only for one karyotype-mode stage."""
    return [
        f"weight{suffix}.npy",
        f"predict_B{suffix}.npy",
        f"virtual_sky{suffix}.png",
        f"virtual_sky{suffix}.pdf",
        f"karyotype{suffix}.txt",
        f"total_cov{suffix}.png",
        f"total_cov{suffix}.pdf",
        f"cn_data{suffix}.pkl",
        f"SV_call_result{suffix}.vcf",
        f"SV_call_result{suffix}.bnd_weights.tsv",
        f"SV_call_result{suffix}.nclose_bnds.tsv",
        f"SKYPE_result{suffix}.bed",
        f"nclose_report{suffix}.tsv",
        f"nclose_sources{suffix}.tsv",
        f"structure_report{suffix}.tsv",
        f"structure_nclose_usage{suffix}.tsv",
    ]


def remove_suffixed_outputs(prefix):
    """Drop karyotype-only results left by an earlier karyotype-mode run."""
    removed = []
    for suffix in STAGE_SUFFIXES:
        for name in suffixed_output_names(suffix):
            path = os.path.join(prefix, name)
            if os.path.isfile(path):
                os.remove(path)
                removed.append(path)
    return removed
