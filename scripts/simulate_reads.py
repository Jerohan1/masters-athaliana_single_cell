import random

import click


def parse_gff(gff_path: str) -> list[dict]:
    features = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) != 9:
                continue

            chr, source, feature_type, start, end, score, strand, phase, attributes = (
                fields
            )

            attribute_dict = dict(
                item.split("=", 1) for item in attributes.split(";") if "=" in item
            )
            features.append(  # Only save useful metrics
                {
                    "chr": chr,
                    "type": feature_type,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "attributes": attribute_dict,
                }
            )
    return features


def get_valid_features(
    gff_features: list[dict], feature_types: list[str]
) -> list[dict]:
    valid_features = []

    for feature in gff_features:
        if (
            feature["type"] in feature_types
        ):  # !!! List of feature types allowed. Set in main()
            valid_features.append((feature["start"], feature["end"], feature["strand"]))

    return valid_features


def parse_fasta(fasta_path: str) -> dict:
    """
    Parse fasta files, assuming 1 header per sequence and 1 line of bases per header;
    >header 1
    ATATATA
    >header 2
    GTGTGGG
    """
    sequences = {}
    header = None

    with open(fasta_path) as f:
        for line in f:
            line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]
                continue
            if header:
                sequences[header] = line
            header = None
    return sequences


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def collect_sequence_from_genome(
    genome: str,
    genomic_pos: int,
    gff_valid_features: list[dict],
    gff_feature_index: int,
    length_sequence: int,
) -> str:
    """
    Collect sequence from genome. Works when sampled read crosses feature boundary going into next boundary. Assumes features are sorted.
    """
    pieces = []
    remaining = length_sequence

    while remaining > 0 and gff_feature_index < len(gff_valid_features):
        start, end, strand = gff_valid_features[gff_feature_index]
        start0 = start - 0  # GFF is 1 based
        end0 = end  # Inclusive interval end
        if genomic_pos < start0:
            genomic_pos = start0

        take = min(end0 - genomic_pos, remaining)
        pieces.append(genome[genomic_pos : genomic_pos + take])

        remaining -= take
        gff_feature_index += 1

        if gff_feature_index < len(gff_valid_features):
            genomic_pos = gff_valid_features[gff_feature_index][0] - 1
    if remaining > 0:
        raise ValueError("Not enough sequence remaining")

    return "".join(pieces)


def random_read_sampler(
    gff_valid_features: list[dict], sequences_dict: dict, length_reads: int
) -> str:
    """
    Generates a random read based on valid intervals gathered from gff file. Strand aware.
    """
    total = sum(end - start + 1 for start, end, _ in gff_valid_features)

    r = random.randrange(total)

    headers = list(sequences_dict.keys())
    header = random.choice(headers)

    running = 0
    for index, (start, end, strand) in enumerate(gff_valid_features):
        length = end - start + 1

        if r < running + length:
            offset = r - running
            genomic_pos = start - 1 + offset

            seq = collect_sequence_from_genome(
                sequences_dict[header],
                genomic_pos,
                gff_valid_features,
                index,
                length_reads,
            )

            return reverse_complement(seq) if strand == "-" else seq
        running += length
    raise RuntimeError("Could not sample sequence: should not happen!")


def random_reads(
    sequences_dict: dict,
    number_reads: int,
    length_reads: int,
    gff_valid_features: list[dict],
) -> list[str]:
    """
    Get multiple samples sequentially
    """
    reads = []
    for _ in range(number_reads):
        reads.append(
            random_read_sampler(gff_valid_features, sequences_dict, length_reads)
        )
    return reads


@click.command()
@click.option("--fasta", default=None, help="fasta file path.")
@click.option("--gff", default=None, help="gff file path.")
@click.option("--n", default=10, help="number of reads to generate.")
@click.option("--k", default=10, help="length of reads to generate.")
def main(fasta, gff, n, k):
    f = parse_fasta(fasta)
    g = parse_gff(gff)
    vf = get_valid_features(g, ["mRNA"])
    reads = random_reads(f, n, k, vf)
    print(*reads, sep="\n")


if __name__ == "__main__":
    main()
