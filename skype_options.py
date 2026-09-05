"""Option routing and the stage-01 settings handoff for native SKYPE."""

import json
import shlex
from pathlib import Path

OPTIONS_FILE = "skype_options.json"


def normalize_extra_args(extra_args):
    if not extra_args:
        return []
    return shlex.split(extra_args) if isinstance(extra_args, str) else list(extra_args)


_STAGE_OPTION_SPECS = {
    "--exclude-nclose-list-loc": ("01", "--exclude-nclose-list-loc", 1),
    "--exclude_nclose_list_loc": ("01", "--exclude-nclose-list-loc", 1),
    "--skip-bam-analysis": ("01", "--skip-bam-analysis", 0),
    "--skip_bam_analysis": ("01", "--skip-bam-analysis", 0),
    "--check-nclose-count": ("01", "--check-nclose-count", 0),
    "--check_nclose_count": ("01", "--check-nclose-count", 0),
    "--nclose-count-vaf-threshold": (
        "01", "--nclose-count-vaf-threshold", 1,
    ),
    "--nclose_count_vaf_threshold": (
        "01", "--nclose-count-vaf-threshold", 1,
    ),
    "--disable-alt-ctg-simple": ("01", "--disable-alt-ctg-simple", 0),
    "--disable_alt_ctg_simple": ("01", "--disable-alt-ctg-simple", 0),
    "--vcf-filter-pass": ("01", "--vcf-filter-pass", "+"),
    "--vcf_filter_pass": ("01", "--vcf-filter-pass", "+"),
    "--debug-force-nclose": ("01", "--debug-force-nclose", 2),
    "--debug_force_nclose": ("01", "--debug-force-nclose", 2),
    "--verbose": ("10", "--verbose", 0),
    "--add-indel-graph": ("10", "--add-indel-graph", 0),
    "--add_indel_graph": ("10", "--add-indel-graph", 0),
    "--limit-combinations": ("10", "--limit-combinations", 1),
    "--limit_combinations": ("10", "--limit-combinations", 1),
}

_REPEATABLE_OPTIONS = {
    "--debug-force-nclose",
}

_ORCHESTRATOR_OWNED_OPTIONS = {
    "-t",
    "--thread",
    "-d",
    "--graph-depth",
    "--graph_depth",
    "--progress",
    "--alt",
    "--original-paf-loc",
    "--original_paf_loc",
    "--vcf-input",
    "--vcf_input",
    "--main-stat-path",
    "--censat-bed-path",
}


def split_stage_options(option_skype):
    """Route SKYPE options to the stage that consumes them.

    Inputs and resource paths remain owned by ACCtools so an extra option
    cannot silently replace the files used to construct the handoff.
    """

    tokens = normalize_extra_args(option_skype)
    stage_args = {"01": [], "10": []}
    seen = set()
    index = 0
    while index < len(tokens):
        raw_token = tokens[index]
        option, separator, inline_value = raw_token.partition("=")
        if option in _ORCHESTRATOR_OWNED_OPTIONS:
            raise ValueError(
                f"{option} is controlled by ACCtools in the native pipeline"
            )
        if option not in _STAGE_OPTION_SPECS:
            raise ValueError(
                f"Unknown --option_skype argument for the native pipeline: "
                f"{option}"
            )
        stage, canonical, arity = _STAGE_OPTION_SPECS[option]
        if canonical in seen and canonical not in _REPEATABLE_OPTIONS:
            raise ValueError(
                f"Duplicate --option_skype argument: {canonical}"
            )
        seen.add(canonical)

        if arity == 0:
            if separator:
                raise ValueError(
                    f"{option} does not accept a value"
                )
            values = []
        elif isinstance(arity, int) and arity > 1:
            if separator:
                raise ValueError(
                    f"{option} requires {arity} separate values"
                )
            if index + arity >= len(tokens):
                raise ValueError(
                    f"{option} requires {arity} values"
                )
            values = tokens[index + 1:index + 1 + arity]
            if any(value.startswith("-") for value in values):
                raise ValueError(
                    f"{option} requires {arity} values"
                )
            index += arity
        elif separator:
            if not inline_value:
                raise ValueError(f"{option} requires a value")
            values = [inline_value]
        elif arity == 1:
            if index + 1 >= len(tokens):
                raise ValueError(f"{option} requires a value")
            index += 1
            values = [tokens[index]]
        else:
            values = []
            while index + 1 < len(tokens):
                next_token = tokens[index + 1]
                next_option = next_token.partition("=")[0]
                if next_option.startswith("-"):
                    break
                index += 1
                values.append(tokens[index])
            if not values:
                raise ValueError(
                    f"{option} requires at least one value"
                )

        stage_args[stage].append(canonical)
        stage_args[stage].extend(values)
        index += 1

    return stage_args["01"], stage_args["10"]


def save_stage_options(prefix, stage01, stage10):
    path = Path(prefix) / OPTIONS_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({"01": stage01, "10": stage10}, indent=2) + "\n",
        encoding="utf-8",
    )


def load_stage_options(prefix, stage):
    path = Path(prefix) / OPTIONS_FILE
    if not path.exists():
        return []
    return json.loads(path.read_text(encoding="utf-8")).get(str(stage), [])
