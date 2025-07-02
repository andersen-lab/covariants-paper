import pandas as pd
import os 


MONTHS = {
    "JAN": "01",
    "FEB": "02",
    "MAR": "03",
    "APR": "04",
    "MAY": "05",
    "JUN": "06",
    "JUL": "07",
    "AUG": "08",
    "SEP": "09",
    "OCT": "10",
    "NOV": "11",
    "DEC": "12",
}

START_DATE = pd.to_datetime("2023-01-01")
END_DATE = pd.to_datetime("2024-12-31")

def clean_sample_name(sample):
    if ".R" in sample:
        sample = sample.split(".R")[0]
        if 'rerun' in sample:
            sample += "_rerun"
    elif "_R" in sample:
        sample = sample.split("_R")[0]
        if 'rerun' in sample:
            sample += "_rerun"
    else:
        sample = sample.split("__")[0]
    return sample

def parse_search_metadata(sample):
    """
    Get collection date and location from SEARCH sample name.
    """
    if "ENC" in sample:
        loc = "Encina"
        date = sample[sample.index("ENC") + 3 : sample.index("ENC") + 8]
        submission_date = sample[: sample.index("ENC")]
    elif "SB" in sample:
        loc = "South Bay"
        date = sample[sample.index("SB") + 2 : sample.index("SB") + 7]
        submission_date = sample[: sample.index("SB")]
    elif "PL" in sample:
        loc = "Point Loma"
        date = sample[sample.index("PL") + 2 : sample.index("PL") + 7]
        submission_date = sample[: sample.index("PL")]

    try:
        collection_month = MONTHS[date[:3]]
        collection_day = date[3:]
    except KeyError:
        print(f"Error: {sample}")
        return None, None
    
    submission_month = submission_date[:2]
    submission_year = submission_date[-3:-1]

    if submission_month == "01" and collection_month == "12":
        collection_year = f"20{int(submission_year)-1}"
    else:
        collection_year = f"20{submission_year}"

    collection_date = f"{collection_year}-{collection_month}-{collection_day}"

    return collection_date, loc

def main():

    cryptics = pd.read_csv("../cryptic_variants.tsv", sep="\t")
    cryptics["sample"] = cryptics["sample"].apply(clean_sample_name)

    cryptics["collection_date"], cryptics["location"] = zip(
        *cryptics["sample"].map(parse_search_metadata)
    )
    cryptics["collection_date"] = pd.to_datetime(cryptics["collection_date"])

    # Filter to study range
    cryptics = cryptics[
        (cryptics["collection_date"] >= START_DATE)
        & (cryptics["collection_date"] <= END_DATE)
    ]

    # Drop rows with missing collection date or location
    cryptics = cryptics.dropna(subset=["collection_date", "location"])

    cryptics = cryptics.sort_values("collection_date")

    # Save combined metadata
    cryptics.to_csv("../cryptic_variants_metadata.tsv", index=False, sep="\t")

if __name__ == '__main__':
    main()