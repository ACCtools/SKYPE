"""Helpers for annotating existing VCF records without rewriting bad input."""


def upsert_vcf_info_fields(info_text, updates):
    update_items = list(updates.items())
    update_keys = {key for key, _value in update_items}
    items = []
    if info_text not in ("", "."):
        for item in info_text.split(";"):
            if not item:
                continue
            key = item.split("=", 1)[0]
            if key not in update_keys:
                items.append(item)
    items.extend(f"{key}={value}" for key, value in update_items)
    return ";".join(items) if items else "."


def write_vcf_record_with_fallback(writer, vcf_record, raw_line, info_updates):
    """Annotate an input record, preserving its raw form if vcfpy rejects it."""
    for info_id in info_updates:
        vcf_record.INFO.pop(info_id, None)
    vcf_record.INFO.update(info_updates)

    try:
        writer.write_record(vcf_record)
    except (TypeError, ValueError) as error:
        columns = raw_line.split("\t")
        if len(columns) < 8:
            raise ValueError(
                "Cannot write VCF fallback record with fewer than 8 columns"
            ) from error
        columns[7] = upsert_vcf_info_fields(columns[7], info_updates)
        print(*columns, sep="\t", file=writer.stream)
        return error
    return None
