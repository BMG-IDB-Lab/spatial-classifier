#!/usr/bin/env python3
"""Parse command-line settings; replace the example prints with your own code."""

import argparse
from pathlib import Path


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Example command-line interface for an analysis script.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Example: python cli_template.py --input data.h5ad "
            "--output results --threads 8 --genes Maf Igf1 --verbose"
        ),
    )

    parser.add_argument(
        "-i", "--input", type=Path, required=True,
        help="Path to the input file (required).",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("results"),
        help="Output directory to use in your analysis code.",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=4,
        help="Number of worker threads to use in your analysis code.",
    )
    parser.add_argument(
        "--genes", nargs="+", default=None, metavar="GENE",
        help="Optional space-separated list of genes.",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable detailed messages in your analysis code.",
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # Access parsed values as args.input, args.output, args.threads, etc.
    # This template only prints settings; add your analysis here.
    print(f"Input:   {args.input}")
    print(f"Output:  {args.output}")
    print(f"Threads: {args.threads}")
    print(f"Genes:   {args.genes}")
    if args.verbose:
        print("Verbose mode is enabled.")


if __name__ == "__main__":
    main()
