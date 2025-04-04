import os
import pandas as pd

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


def parse_metadata(sample):

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

    month = MONTHS[date[:3]]
    day = date[3:]

    
    submission_month = submission_date[:2]
    submission_year = submission_date[-3:-1]

    if submission_month == "01" and month == "12":
        year = f"20{int(submission_year)-1}"
    else:
        year = f"20{submission_year}"

    collection_date = f"{year}-{month}-{day}"

    return collection_date, loc

def main():
    metadata = pd.DataFrame(columns=["sample", "collection_date", "location"])

    metadata["sample"] = [f for f in os.listdir("covariants")]

    metadata["collection_date"], metadata["location"] = zip(
        *metadata["sample"].map(parse_metadata)
    )
    metadata["collection_date"] = pd.to_datetime(metadata["collection_date"])
    metadata = metadata.sort_values("collection_date")
    metadata.to_csv("data/SEARCH_metadata.tsv", index=False, sep="\t")

if __name__ == '__main__':
    main()