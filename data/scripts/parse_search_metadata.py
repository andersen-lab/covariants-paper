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

START_DATE = "2023-01-01"
END_DATE = "2024-12-31"

def parse_metadata(sample):
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
        month = MONTHS[date[:3]]
        day = date[3:]
    except KeyError:
        print(f"Error: {sample}")
        return None, None
    
    submission_month = submission_date[:2]
    submission_year = submission_date[-3:-1]

    if submission_month == "01" and month == "12":
        year = f"20{int(submission_year)-1}"
    else:
        year = f"20{submission_year}"

    collection_date = f"{year}-{month}-{day}"

    return collection_date, loc

def main():
    BAM_DIR = "../bam/"
    combined_metadata = pd.DataFrame(columns=["sample", "collection_date", "location"])
    samples = []
    # Read metadata from bam filename
    for file in os.listdir(BAM_DIR):
        if file.endswith(".bam"):
            if any(x in file for x in ["rerun", "separate", "combined"]):
                continue
            samples.append(file.split(".trimmed")[0])

    combined_metadata["sample"] = samples            
    combined_metadata["collection_date"], combined_metadata["location"] = zip(
        *combined_metadata["sample"].map(parse_metadata)
    )
    combined_metadata["collection_date"] = pd.to_datetime(combined_metadata["collection_date"])

    # Filter to study range
    combined_metadata = combined_metadata[
        (combined_metadata["collection_date"] >= START_DATE)
        & (combined_metadata["collection_date"] <= END_DATE)
    ]

    combined_metadata = combined_metadata.sort_values("collection_date")

    # Save combined metadata
    combined_metadata.to_csv("../search_metadata.csv", index=False)

if __name__ == '__main__':
    main()