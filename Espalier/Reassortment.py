##########################
##
## Espalier: A Python package for tree reconciliation and reconstructing ARGs using maximum agreement forests.
##
## Copyright 2021-2022 David A. Rasmussen (drasmus@ncsu.edu)
##
## If using Espalier or this code, please cite:
##
##      Rasmussen, D.A. and Guo, F. Espalier: Efficient tree reconciliation and ARG reconstruction using maximum agreement forests. 2022.
##
############################

import csv
import os
import random
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional, Sequence

import dendropy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

ANALYSIS_MODE_RECOMBINATION = "recombination"
ANALYSIS_MODE_REASSORTMENT = "reassortment"
ANALYSIS_MODES = {
    ANALYSIS_MODE_RECOMBINATION,
    ANALYSIS_MODE_REASSORTMENT,
}
DATE_FORMAT = "%Y-%m-%d"


@dataclass(frozen=True)
class HeaderIdentity:
    full_label: str
    isolate_key: str
    analysis_label: str
    sample_date: Optional[datetime]


@dataclass
class SegmentData:
    seq_file: str
    tree_file: Optional[str]
    segment_name: str
    records_by_full_label: Dict[str, SeqRecord]
    records_by_analysis_label: Dict[str, SeqRecord]
    identities_by_full_label: Dict[str, HeaderIdentity]
    identities_by_analysis_label: Dict[str, HeaderIdentity]
    identities_by_tree_label: Dict[str, HeaderIdentity]


@dataclass
class PreparedAnalysisInputs:
    analysis_mode: str
    seq_files: List[str]
    tree_files: Optional[List[str]]
    derived_dir: Optional[str] = None
    normalized_alignments_dir: Optional[str] = None
    normalized_trees_dir: Optional[str] = None
    label_map_path: Optional[str] = None
    common_isolates: Optional[List[str]] = None


def sanitize_analysis_label(label: str) -> str:
    normalized = re.sub(r"\s+", "-", label.strip())
    if not normalized:
        raise ValueError("Isolate label cannot be empty after normalization")
    return normalized


def parse_date_from_label(label: str) -> Optional[datetime]:
    fields = label.split("|")
    if not fields:
        return None
    try:
        return datetime.strptime(fields[-1].strip(), DATE_FORMAT)
    except ValueError:
        return None


def get_iso_week(sample_date: datetime) -> str:
    iso_cal = sample_date.isocalendar()
    return f"{iso_cal[0]}-W{iso_cal[1]:02d}"


def parse_header_identity(label: str, isolate_field: int = 0) -> HeaderIdentity:
    full_label = label.lstrip(">").strip()
    fields = [field.strip() for field in full_label.split("|")]
    if isolate_field >= len(fields):
        raise ValueError(
            f"Cannot extract isolate field {isolate_field} from FASTA label '{full_label}'"
        )
    isolate_key = fields[isolate_field]
    if not isolate_key:
        raise ValueError(f"Empty isolate key extracted from FASTA label '{full_label}'")
    return HeaderIdentity(
        full_label=full_label,
        isolate_key=isolate_key,
        analysis_label=sanitize_analysis_label(isolate_key),
        sample_date=parse_date_from_label(full_label),
    )


def _tree_label_aliases(label: str) -> List[str]:
    aliases = []
    for candidate in (label, label.replace("_", " "), label.replace(" ", "_")):
        if candidate not in aliases:
            aliases.append(candidate)
    return aliases


def find_matching_tree(seq_file: str, tree_dir: str) -> Optional[str]:
    base = os.path.basename(seq_file)
    patterns = [
        base + ".tre",
        base + ".rooted.tre",
        base.replace(".fasta", ".tre").replace(".aln", ".tre"),
        base.replace(".fasta", ".rooted.tre").replace(".aln", ".rooted.tre"),
        base.replace("_aligned_common.fasta", "_tree_common_binary.tre"),
        base.replace("_aligned_common.fasta", "_tree_common.tre"),
    ]
    for pattern in patterns:
        tree_file = os.path.join(tree_dir, pattern)
        if os.path.exists(tree_file):
            return tree_file

    prefix = base.split(".")[0]
    for candidate in sorted(os.listdir(tree_dir)):
        if candidate.startswith(prefix) and candidate.endswith(".tre"):
            return os.path.join(tree_dir, candidate)

    return None


def load_segment_data(seq_file: str, isolate_field: int = 0, tree_file: Optional[str] = None) -> SegmentData:
    records_by_full_label: Dict[str, SeqRecord] = {}
    records_by_analysis_label: Dict[str, SeqRecord] = {}
    identities_by_full_label: Dict[str, HeaderIdentity] = {}
    identities_by_analysis_label: Dict[str, HeaderIdentity] = {}
    identities_by_tree_label: Dict[str, HeaderIdentity] = {}

    for record in SeqIO.parse(seq_file, "fasta"):
        full_label = record.description.strip()
        identity = parse_header_identity(full_label, isolate_field=isolate_field)
        if full_label in records_by_full_label:
            raise ValueError(f"Duplicate FASTA header '{full_label}' in {seq_file}")
        existing_identity = identities_by_analysis_label.get(identity.analysis_label)
        if existing_identity is not None:
            raise ValueError(
                "Duplicate isolate key "
                f"'{identity.isolate_key}' in {seq_file}; "
                f"headers '{existing_identity.full_label}' and '{full_label}' collide "
                f"after normalization to '{identity.analysis_label}'"
            )
        records_by_full_label[full_label] = record
        records_by_analysis_label[identity.analysis_label] = record
        identities_by_full_label[full_label] = identity
        identities_by_analysis_label[identity.analysis_label] = identity
        for alias in _tree_label_aliases(full_label):
            existing_tree_identity = identities_by_tree_label.get(alias)
            if (
                existing_tree_identity is not None
                and existing_tree_identity.full_label != identity.full_label
            ):
                raise ValueError(
                    "Ambiguous tree label alias "
                    f"'{alias}' in {seq_file}; "
                    f"headers '{existing_tree_identity.full_label}' and "
                    f"'{full_label}' would both match"
                )
            identities_by_tree_label[alias] = identity

    return SegmentData(
        seq_file=seq_file,
        tree_file=tree_file,
        segment_name=os.path.basename(seq_file),
        records_by_full_label=records_by_full_label,
        records_by_analysis_label=records_by_analysis_label,
        identities_by_full_label=identities_by_full_label,
        identities_by_analysis_label=identities_by_analysis_label,
        identities_by_tree_label=identities_by_tree_label,
    )


def _load_tree(tree_file: str) -> dendropy.Tree:
    return dendropy.Tree.get(
        path=tree_file,
        schema="newick",
        rooting="default-rooted",
        preserve_underscores=True,
    )


def resolve_tree_identity(segment_data: SegmentData, tree_label: str) -> HeaderIdentity:
    label = tree_label.strip()
    identity = segment_data.identities_by_tree_label.get(label)
    if identity is None:
        raise ValueError(
            f"Tree tip '{label}' from {segment_data.tree_file} does not match "
            f"any alignment record in {segment_data.seq_file}"
        )
    return identity


def validate_tree_matches_segment(tree_file: str, segment_data: SegmentData) -> None:
    tree = _load_tree(tree_file)
    tree_tip_labels = []
    resolved_alignment_labels = set()
    for leaf in tree.leaf_node_iter():
        tree_label = leaf.taxon.label.strip()
        tree_tip_labels.append(tree_label)
        identity = segment_data.identities_by_tree_label.get(tree_label)
        if identity is not None:
            resolved_alignment_labels.add(identity.full_label)

    alignment_labels = set(segment_data.records_by_full_label)

    missing_in_alignment = sorted(
        label for label in tree_tip_labels if label not in segment_data.identities_by_tree_label
    )
    missing_in_tree = sorted(alignment_labels - resolved_alignment_labels)
    if missing_in_alignment or missing_in_tree:
        parts = [f"Tree/alignment mismatch for {segment_data.segment_name}"]
        if missing_in_alignment:
            parts.append(f"tree-only labels: {missing_in_alignment[:5]}")
        if missing_in_tree:
            parts.append(f"alignment-only labels: {missing_in_tree[:5]}")
        raise ValueError("; ".join(parts))


def _write_normalized_alignment(segment_data: SegmentData, selected_labels: Sequence[str], output_file: str) -> None:
    selected_set = set(selected_labels)
    normalized_records: List[SeqRecord] = []
    for full_label, record in segment_data.records_by_full_label.items():
        identity = segment_data.identities_by_full_label[full_label]
        if identity.analysis_label not in selected_set:
            continue
        normalized_records.append(
            SeqRecord(
                record.seq,
                id=identity.analysis_label,
                name=identity.analysis_label,
                description=identity.analysis_label,
            )
        )
    SeqIO.write(normalized_records, output_file, "fasta")


def _write_normalized_tree(segment_data: SegmentData, selected_labels: Sequence[str], output_file: str) -> None:
    tree = _load_tree(segment_data.tree_file)
    selected_set = set(selected_labels)
    for leaf in tree.leaf_node_iter():
        identity = resolve_tree_identity(segment_data, leaf.taxon.label)
        leaf.taxon.label = identity.analysis_label

    taxa_to_remove = [taxon for taxon in tree.taxon_namespace if taxon.label not in selected_set]
    if taxa_to_remove:
        tree.prune_taxa(taxa_to_remove)
    tree.suppress_unifurcations()
    tree.write(path=output_file, schema="newick", suppress_annotations=True, suppress_rooting=True)


def prepare_analysis_inputs(
    seq_files: Sequence[str],
    tree_files: Optional[Sequence[str]],
    output_dir: str,
    analysis_mode: str = ANALYSIS_MODE_RECOMBINATION,
    isolate_field: int = 0,
) -> PreparedAnalysisInputs:
    if analysis_mode not in ANALYSIS_MODES:
        raise ValueError(f"Unsupported analysis mode '{analysis_mode}'")

    if analysis_mode == ANALYSIS_MODE_RECOMBINATION:
        return PreparedAnalysisInputs(
            analysis_mode=analysis_mode,
            seq_files=list(seq_files),
            tree_files=list(tree_files) if tree_files else None,
        )

    if tree_files is not None and len(tree_files) != len(seq_files):
        raise ValueError("seq_files and tree_files must have the same length in reassortment mode")

    derived_dir = os.path.join(output_dir, "derived")
    normalized_alignments_dir = os.path.join(derived_dir, "normalized_alignments")
    normalized_trees_dir = os.path.join(derived_dir, "normalized_trees")
    os.makedirs(normalized_alignments_dir, exist_ok=True)
    os.makedirs(normalized_trees_dir, exist_ok=True)

    segment_data: List[SegmentData] = []
    for idx, seq_file in enumerate(seq_files):
        tree_file = tree_files[idx] if tree_files else None
        data = load_segment_data(seq_file, isolate_field=isolate_field, tree_file=tree_file)
        if tree_file:
            validate_tree_matches_segment(tree_file, data)
        segment_data.append(data)

    common_isolates = set(segment_data[0].records_by_analysis_label)
    for data in segment_data[1:]:
        common_isolates &= set(data.records_by_analysis_label)
    if not common_isolates:
        raise ValueError("No shared isolates found across segments in reassortment mode")

    selected_labels = sorted(common_isolates)
    normalized_seq_files: List[str] = []
    normalized_tree_files: List[str] = []
    label_map_path = os.path.join(derived_dir, "label_map.tsv")

    with open(label_map_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "segment",
                "original_full_label",
                "isolate_key",
                "analysis_label",
                "selected",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for data in segment_data:
            normalized_seq_file = os.path.join(
                normalized_alignments_dir,
                os.path.basename(data.seq_file),
            )
            _write_normalized_alignment(data, selected_labels, normalized_seq_file)
            normalized_seq_files.append(normalized_seq_file)

            if data.tree_file:
                normalized_tree_file = os.path.join(
                    normalized_trees_dir,
                    os.path.basename(data.tree_file),
                )
                _write_normalized_tree(data, selected_labels, normalized_tree_file)
                normalized_tree_files.append(normalized_tree_file)

            for full_label, identity in data.identities_by_full_label.items():
                writer.writerow(
                    {
                        "segment": data.segment_name,
                        "original_full_label": full_label,
                        "isolate_key": identity.isolate_key,
                        "analysis_label": identity.analysis_label,
                        "selected": "1" if identity.analysis_label in common_isolates else "0",
                    }
                )

    return PreparedAnalysisInputs(
        analysis_mode=analysis_mode,
        seq_files=normalized_seq_files,
        tree_files=normalized_tree_files if normalized_tree_files else None,
        derived_dir=derived_dir,
        normalized_alignments_dir=normalized_alignments_dir,
        normalized_trees_dir=normalized_trees_dir,
        label_map_path=label_map_path,
        common_isolates=selected_labels,
    )


def select_isolates_by_week(segment_data: Sequence[SegmentData], seed: int = 42) -> List[str]:
    if not segment_data:
        raise ValueError("No segment data supplied for weekly isolate selection")

    common_isolates = set(segment_data[0].records_by_analysis_label)
    for data in segment_data[1:]:
        common_isolates &= set(data.records_by_analysis_label)
    if not common_isolates:
        raise ValueError("No shared isolates found across segments")

    weeks: Dict[str, List[str]] = defaultdict(list)
    excluded_without_dates: List[str] = []
    for isolate in sorted(common_isolates):
        sample_dates = {
            data.identities_by_analysis_label[isolate].sample_date
            for data in segment_data
        }
        if None in sample_dates:
            if len(sample_dates) > 1:
                raise ValueError(f"Inconsistent dates for isolate '{isolate}' across segments")
            excluded_without_dates.append(isolate)
            continue
        if len(sample_dates) != 1:
            raise ValueError(f"Inconsistent dates for isolate '{isolate}' across segments")
        sample_date = next(iter(sample_dates))
        weeks[get_iso_week(sample_date)].append(isolate)

    rng = random.Random(seed)
    selected_isolates = [rng.choice(weeks[week]) for week in sorted(weeks)]
    if not selected_isolates:
        raise ValueError(
            "No shared isolates with parseable dates were available for weekly subsampling"
        )
    return selected_isolates
