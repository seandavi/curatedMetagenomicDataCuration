import polars as pl
import pathlib


def download_sra_run_members(replace: bool = False) -> None:
    """Download SRA Run Members file"""
    if replace or not pathlib.Path("SRA_Run_Members.tsv.gz").exists():
        pathlib.Path("SRA_Run_Members.tsv.gz").unlink(missing_ok=True)
        pl.read_csv(
            "https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Run_Members.tab",
            separator="\t",
            has_header=True,
            null_values="-",
            infer_schema_length=1000000,
        ).write_csv("SRA_Run_Members.tsv.gz")


def load_sra_run_members() -> pl.DataFrame:
    """Load SRA Run Members"""
    return (
        pl.read_csv(
            "SRA_Run_Members.tsv.gz",
            separator="\t",
            has_header=True,
            null_values="-",
            infer_schema_length=1000000,
        )
        .with_columns(
            pl.col("Spots").cast(pl.Int64),
            pl.col("Bases").cast(pl.Int64),
        )
        .filter(pl.col("Status") == "live")
        .drop("Status", "Member_Name")
    )


def load_curated_data(f: pathlib.Path) -> pl.DataFrame:
    """Load curated data"""
    return (
        pl.read_csv(
            f,
            separator="\t",
            has_header=True,
            null_values="-",
            infer_schema_length=10000,
        )
        .rename(str.lower)
        .rename(str.strip)
    )


def add_study_and_sample_info(
    df: pl.DataFrame, df1: pl.DataFrame, f: pathlib.Path
) -> None:
    """Add study and sample info to a curated dataframe"""
    df2 = pl.DataFrame()
    if "ncbi_accession" in df.columns:
        df2 = (
            df.select("sample_id", "ncbi_accession")
            .with_columns(pl.col("ncbi_accession").str.split(";"))
            .explode("ncbi_accession")
            .join(df1, right_on="Run", left_on="ncbi_accession", how="left")
            .group_by("sample_id")
            .agg(
                pl.col("ncbi_accession"),
                pl.col("Spots").sum().alias("spots"),
                pl.col("Bases").sum().alias("bases"),
                pl.col("BioSample").first().alias("biosample"),
                pl.col("Experiment").first().alias("sra_experiment"),
                pl.col("Sample").first().alias("sra_sample"),
                pl.col("Study").first().alias("sra_study"),
            )
            .join(df.drop("ncbi_accession"), on="sample_id", how="right")
        )
    elif "sample_id" in df.columns and df.item(1, "sample_id").startswith("SAM"):
        df2 = (
            df.select("sample_id")
            .join(df1, right_on="BioSample", left_on="sample_id", how="left")
            .group_by("sample_id")
            .agg(
                pl.col("Spots").sum().alias("spots"),
                pl.col("Bases").sum().alias("bases"),
                pl.col("sample_id").first().alias("biosample"),
                pl.col("Experiment").first().alias("sra_experiment"),
                pl.col("Sample").first().alias("sra_sample"),
                pl.col("Study").first().alias("sra_study"),
            )
            .join(df, on="sample_id", how="right")
        )
    df2.write_ndjson(f.with_suffix(".ndjson"))
    assert df2.shape[0] == df.shape[0]


def main() -> None:
    """Main function"""
    download_sra_run_members(replace=False)
    df1 = load_sra_run_members()
    files = pathlib.Path("inst/curated").glob("**/*.tsv")
    for f in files:
        print(f)
        df = load_curated_data(f)
        add_study_and_sample_info(df, df1, f)


if __name__ == "__main__":
    main()
