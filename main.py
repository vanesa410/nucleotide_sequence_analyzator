separator = "-" * 60


def get_sequence() -> str:
    """
    Gets the sequence from the user and cleans it by removing spaces and other whitespace characters.

    Returns:
        string: The sequence that is going to be analyzed.
    """
    user_sequence = input("Enter your DNA / RNA sequence in FASTA format: ")
    sequence = "".join(user_sequence.split())  # sequence without gaps
    # dodělat odstranění FASTA hlaviček
    return sequence


def is_valid(sequence: str) -> bool:
    """
    Validates the sequence, it cannot contain both thymine (T) and uracil (U).

    Returns:
        bool
    """
    upper_sequence = sequence.upper()
    # dodělat další podmínky
    if "U" in upper_sequence and "T" in upper_sequence:
        return False
    return True


def base_counter(sequence: str) -> tuple[dict[str, int], int]:
    """
    Counts the bases in a nucleotide sequence.

    Returns:
        tuple: A dictionary with base counts and an integer with the total base count.
    """
    adenine_count = 0
    guanine_count = 0
    cytosine_count = 0
    thymine_count = 0
    uracil_count = 0
    total_base_count = 0

    for base in sequence:
        if base.upper() == "A":
            adenine_count += 1
        elif base.upper() == "G":
            guanine_count += 1
        elif base.upper() == "C":
            cytosine_count += 1
        elif base.upper() == "T":
            thymine_count += 1
        elif base.upper() == "U":
            uracil_count += 1
        else:
            print(f"Invalid base: {base}")  # dokončit
    total_base_count = sum(
        (adenine_count, guanine_count, cytosine_count, thymine_count, uracil_count)
    )
    return {
        "Adenine": adenine_count,
        "Guanine": guanine_count,
        "Cytosine": cytosine_count,
        "Thymine": thymine_count,
        "Uracil": uracil_count,
    }, total_base_count


def gc_content(sequence: str) -> tuple[float, float]:
    """
    Counts the GC content in a nucleotide sequence.

    Returns:
        tuple: A float with GC fraction and GC percentage.
    """
    base_counts, total_base_count = base_counter(sequence)
    gc_fraction = round(
        ((base_counts["Guanine"] + base_counts["Cytosine"]) / total_base_count), 2
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
    base_counts, total_base_count = base_counter(sequence)
    purines_count = base_counts["Adenine"] + base_counts["Guanine"]
    purines_percentage = round((purines_count / total_base_count * 100), 2)
    pyrimidines_count = (
        base_counts["Cytosine"] + base_counts["Thymine"] + base_counts["Uracil"]
    )
    pyrimidines_percentage = round((pyrimidines_count / total_base_count * 100), 2)

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


def get_results(sequence: str) -> None:
    """
    Main function that prints the results.
    """
    if is_valid(sequence):
        print(separator)
        base_counts, total_base_count = base_counter(sequence)
        gc_fraction, gc_percentage = gc_content(sequence)
        print(f"Total base count - length of the sequence: {total_base_count} bases")
        print(f"Count of bases: {base_counts}")
        print(f"GC percentage: {gc_percentage} %\n{separator}")

        for group, values in purines_pyrimidines(sequence).items():
            for key, value in values.items():
                if "percentage" in key:
                    print(f"{key}: {value} %")
                else:
                    print(f"{key}: {value}")

    else:
        print("Invalid sequence. The sequence contains uracil and thymine")


if __name__ == "__main__":
    sequence = get_sequence()
    get_results(sequence)
