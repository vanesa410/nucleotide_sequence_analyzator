separator = "-" * 60


def get_sequence() -> str:
    """
    Gets the sequence from the user and cleans it by removing spaces and other whitespace characters.

    Returns:
        string: The sequence that is going to be analyzed.
    """
    user_sequence = input("Enter your DNA / RNA sequence in FASTA format: ")
    lines = user_sequence.splitlines()
    FASTA_header = lines[0]
    sequence_lines = [line for line in lines if not line.startswith(">")]
    sequence = "".join(sequence_lines)  # sequence without gaps and headers
    return sequence


def sequence_type(sequence: str) -> str:
    """
    Analyses the sequence by the bases and decides if its DNA or RNA.

    Returns:
        str: A molecule type of the sequence (DNA / RNA)
    """
    if "T" in sequence.upper():
        return "DNA"

    elif "U" in sequence.upper():
        return "RNA"
    else:
        return "Unknown"


def is_valid(sequence: str) -> bool:
    """
    Validates the sequence, it has to contain only A,G,C,T,U and it cannot contain both thymine (T) and uracil (U).

    Returns:
        bool
    """
    valid_bases = "AGCTURYSWKMBDHVN"
    gaps = ".-"
    upper_sequence = sequence.upper()

    # Check for invalid characters
    for base in upper_sequence:
        if base not in valid_bases and base not in gaps:
            return False

    # Sequence cannot contain both T and U
    if "U" in upper_sequence and "T" in upper_sequence:
        return False
    return True


def base_counter(sequence: str) -> dict[str, int]:
    """
    Counts the bases in a nucleotide sequence.

    Returns:
        dict: A dictionary with the results.
    """

    adenine_count = 0
    guanine_count = 0
    cytosine_count = 0
    thymine_count = 0
    uracil_count = 0
    r_count = 0
    y_count = 0
    s_count = 0
    w_count = 0
    k_count = 0
    m_count = 0
    b_count = 0
    d_count = 0
    h_count = 0
    v_count = 0
    n_count = 0
    gap_count = 0

    for base in sequence.upper():
        if base == "A":
            adenine_count += 1
        elif base == "G":
            guanine_count += 1
        elif base == "C":
            cytosine_count += 1
        elif base == "T":
            thymine_count += 1
        elif base == "U":
            uracil_count += 1
        elif base == "R":
            r_count += 1  # R... A or G
        elif base == "Y":
            y_count += 1  # Y... C or T
        elif base == "S":
            s_count += 1  # S... G or C
        elif base == "W":
            w_count += 1  # W... A or T
        elif base == "K":
            k_count += 1  # K... G or T
        elif base == "M":
            m_count += 1  # M... A or C
        elif base == "B":
            b_count += 1  # B... C or G or T
        elif base == "D":
            d_count += 1  # D... A or G or T
        elif base == "H":
            h_count += 1  # H... A or C or T
        elif base == "V":
            v_count += 1  # V... A or C or G
        elif base == "N":
            n_count += 1  # N... any base
        elif base == "-" or base == ".":
            gap_count += 1  # gap
        else:
            raise ValueError(f"Invalid base: {base}")
    total_length = sum(
        1 for base in sequence.upper()
    )  # length of the sequence including gaps
    total_base_count = sum(1 for base in sequence.upper() if base not in "-.")
    ambiguous_base_count = sum(1 for base in sequence.upper() if base in "RYSWKMBDHVN")
    known_base_count = sum(
        (
            adenine_count,
            guanine_count,
            cytosine_count,
            thymine_count,
            uracil_count,
        ),
    )

    return {
        "Adenine": adenine_count,
        "Guanine": guanine_count,
        "Cytosine": cytosine_count,
        "Thymine": thymine_count,
        "Uracil": uracil_count,
        "Gaps": gap_count,
        "Known base count": known_base_count,
        "Ambiguous base count": ambiguous_base_count,
        "Total base count": total_base_count,
        "Total length": total_length,
    }


def gc_content(sequence: str) -> tuple[float, float]:
    """
    Counts the GC content in a nucleotide sequence.

    Returns:
        tuple: A float with GC fraction and GC percentage.
    """
    base_counts = base_counter(sequence)
    if base_counts["Known base count"] == 0:
        raise ValueError("Sequence does not contain any known bases")
    else:
        gc_fraction = round(
            (
                (base_counts["Guanine"] + base_counts["Cytosine"])
                / base_counts["Known base count"]
            ),
            2,
        )
        gc_percentage = round((gc_fraction * 100), 2)
        return gc_fraction, gc_percentage


def purines_pyrimidines(
    sequence: str,
) -> dict[str, dict[str, float | int]]:
    """
    Counts the purines (A + G) and pyrimidines (C + T/U).

    Returns:
        dict: A dictionary with two keys, 'Purines' and 'Pyrimidines', each containing
              a count and a percentage of the respective nucleotide group.
    """
    base_counts = base_counter(sequence)
    purines_count = base_counts["Adenine"] + base_counts["Guanine"]
    pyrimidines_count = (
        base_counts["Cytosine"] + base_counts["Thymine"] + base_counts["Uracil"]
    )
    if base_counts["Known base count"] == 0:
        raise ValueError("Sequence does not contain any known bases.")
    else:
        purines_percentage = round(
            (purines_count / base_counts["Known base count"] * 100), 2
        )
        pyrimidines_percentage = round(
            (pyrimidines_count / base_counts["Known base count"] * 100), 2
        )

    return {
        "Purines": {
            "Purines count": purines_count,
            "Purines percentage": purines_percentage,
        },
        "Pyrimidines": {
            "Pyrimidines count": pyrimidines_count,
            "Pyrimidines percentage": pyrimidines_percentage,
        },
    }


# TRANSCRIPTION


def transcription(sequence: str) -> str:
    """
    Transcribes DNA into mRNA (T -> U).
    If input is already RNA, returns it unchanged.
    If type is Unknown, returns a warning string.

    Returns:
        str
    """
    seq_type = sequence_type(sequence)
    if seq_type == "DNA":
        mRNA_sequence = (sequence.upper()).replace("T", "U")
        return mRNA_sequence
    elif seq_type == "RNA":
        return sequence.upper()
    else:
        raise ValueError("Unknown sequence type")


def get_results(sequence: str) -> dict:
    """
    Returns all of the results.
    Returns:
        dict
    """
    if not is_valid(sequence):
        raise ValueError("Invalid sequence")

    seq_type = sequence_type(sequence)
    base_counts = base_counter(sequence)
    gc_fraction, gc_percentage = gc_content(sequence)
    base_types = purines_pyrimidines(sequence)
    transcript = transcription(sequence)
    results = {
        "sequence_type": seq_type,
        "base_counts": base_counts,
        "gc_fraction": gc_fraction,
        "gc_percentage": gc_percentage,
        "base_types": base_types,
        "transcript": transcript,
    }
    return results


def print_results(results: dict) -> None:
    """
    Prints the results.
    Returns:
        None
    """
    print(separator)
    print(f"SEQUENCE ANALYSIS\n")
    print(f"Type of the sequence: {results['sequence_type']}")
    print(
        f"Length of the sequence including gaps: {results['base_counts']['Total length']}"
    )
    print(f"Total base count: {results['base_counts']['Total base count']} bases")
    print(
        f"Known base count (A, G, C, T, U): {results['base_counts']['Known base count']} bases"
    )
    print(
        f"Ambiguous base count: {results['base_counts']['Ambiguous base count']} bases"
    )
    print(
        f"GC percentage (counts only with the known bases): {results['gc_percentage']} %"
    )
    print(f"{separator}\nAdenine: {results['base_counts']['Adenine']} bases")
    print(f"Guanine: {results['base_counts']['Guanine']} bases")
    print(f"Cytosine: {results['base_counts']['Cytosine']} bases")
    print(f"Thymine: {results['base_counts']['Thymine']} bases")
    print(f"Uracil: {results['base_counts']['Uracil']} bases\n{separator}")

    print("Purines and pyrimidines in the sequence (counts only with the known bases)")
    for group, values in results["base_types"].items():
        for key, value in values.items():
            if "percentage" in key:
                print(f"{key}: {value} %")
            else:
                print(f"{key}: {value}")
    print(f"{separator}\nTranscript: {results['transcript']}")


if __name__ == "__main__":
    sequence = get_sequence()
    results = get_results(sequence)
    print_results(results)
