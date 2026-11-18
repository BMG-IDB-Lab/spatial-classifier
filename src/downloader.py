#!/usr/bin/env python3
import argparse
from http.client import HTTPException
from pathlib import Path, PurePosixPath
import shutil
import sys
import time
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET


BASE_URL = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/"
SAMPLES = tuple("Zhuang-ABCA-{}".format(i) for i in range(1, 5))
SAMPLE_PREFIXES = {
    sample: {
        "expression": "expression_matrices/{}/20230830/".format(sample),
        "metadata": "metadata/{}/20241115/".format(sample),
        "ccf": "metadata/{}-CCF/20230830/".format(sample),
    }
    for sample in SAMPLES
}
# Retain the original constant for callers of the single-sample downloader.
COMPONENTS = SAMPLE_PREFIXES["Zhuang-ABCA-2"]
LOCAL_DIRS = {
    "expression": "expression_matrices",
    "metadata": "metadata",
    "ccf": "ccf_coordinates",
}
PARCELLATION_KEY = (
    "metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/"
    "cell_metadata_with_parcellation_annotation.csv"
)
TAXONOMY_PREFIX = "metadata/WMB-taxonomy/20231215/"
TAXONOMY_DIR = "WMB_taxonomy"
NS = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}
TIMEOUT = 60
ATTEMPTS = 3
SOURCE_ERRORS = (OSError, HTTPException, ET.ParseError, ValueError)


class DownloadError(RuntimeError):
    """At least one requested component could not be listed or downloaded."""


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                    argparse.RawDescriptionHelpFormatter):
    pass


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Download Zhuang-ABCA-1, -2, -3 and -4 (raw .h5ad expression matrices, "
            "metadata, CCF coordinates), plus the requested parcellation CSV "
            "and mouse WMB-taxonomy/20231215 files into WMB_taxonomy/. "
            "Each sample has its own folder. Uses only Python's standard library."
        ),
        formatter_class=HelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 download_zhuang_abca2_https.py --output data\n"
            "  python3 download_zhuang_abca2_https.py --dry-run\n"
            "  python3 download_zhuang_abca2_https.py --samples Zhuang-ABCA-1 Zhuang-ABCA-3\n"
            "  python3 download_zhuang_abca2_https.py --components metadata ccf --skip-parcellation\n\n"
            "Rerunning skips existing files with the expected byte count. "
            "Interrupted individual files restart from the beginning. "
            "Incomplete files use a .part suffix; successful downloads are "
            "checked for the expected size before being renamed."
        ),
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("data"),
        help="Root directory for sample folders and a separate parcellation folder.",
    )
    parser.add_argument(
        "--samples", nargs="+", choices=SAMPLES, default=list(SAMPLES),
        metavar="SAMPLE",
        help="Sample names Zhuang-ABCA-1 through Zhuang-ABCA-4; downloads all by default.",
    )
    parser.add_argument(
        "--skip-parcellation", action="store_true",
        help="Skip the separate MERFISH-C57BL6J-638850-CCF parcellation table.",
    )
    parser.add_argument(
        "--skip-taxonomy", action="store_true",
        help="Skip the mouse WMB taxonomy files otherwise saved in WMB_taxonomy/.",
    )
    parser.add_argument(
        "--components", nargs="+", choices=tuple(COMPONENTS),
        default=list(COMPONENTS), metavar="COMPONENT",
        help=("Components: expression, metadata, ccf. Expression includes only "
              "*-raw.h5ad; log2 matrices are skipped. Downloads all components by default."),
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="List all file URLs and sizes without downloading their contents.",
    )
    return parser.parse_args(argv)


def object_url(key):
    return BASE_URL + quote(key, safe="/")


def list_files(prefix):
    """Read public S3 XML listings, including continuation pages."""
    token = None
    while True:
        params = {"list-type": "2", "prefix": prefix}
        if token is not None:
            params["continuation-token"] = token
        request = Request(BASE_URL + "?" + urlencode(params))
        with urlopen(request, timeout=TIMEOUT) as response:
            root = ET.fromstring(response.read())
        for item in root.findall("s3:Contents", NS):
            key = item.findtext("s3:Key", namespaces=NS)
            size = int(item.findtext("s3:Size", namespaces=NS))
            parts = PurePosixPath(key).parts
            if not key.startswith(prefix) or key.startswith("/") or ".." in parts:
                raise ValueError("Unexpected path in remote listing: " + str(key))
            if not key.endswith("/"):
                yield key, size
        if root.findtext("s3:IsTruncated", namespaces=NS) != "true":
            break
        token = root.findtext("s3:NextContinuationToken", namespaces=NS)
        if not token:
            raise ValueError("Remote listing is incomplete: missing continuation token.")


def download_file(key, expected_size, output, relative_path=None):
    """Transfer one object; optionally map it to a sample-specific local path."""
    destination = Path(output) / (relative_path if relative_path is not None else key)
    if destination.is_file() and destination.stat().st_size == expected_size:
        print("Skip completed: " + str(destination), flush=True)
        return destination

    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    for attempt in range(1, ATTEMPTS + 1):
        print(
            "Download ({}/{}): {} ({:.2f} MB)".format(
                attempt, ATTEMPTS, key, expected_size / 1e6
            ),
            flush=True,
        )
        try:
            request = Request(object_url(key), headers={"Accept-Encoding": "identity"})
            with urlopen(request, timeout=TIMEOUT) as response:
                with partial.open("wb") as target:
                    shutil.copyfileobj(response, target, length=1024 * 1024)
            actual_size = partial.stat().st_size
            if actual_size != expected_size:
                raise OSError(
                    "Incomplete transfer: expected {} bytes, received {}.".format(
                        expected_size, actual_size
                    )
                )
            partial.replace(destination)
            print("Saved: " + str(destination), flush=True)
            return destination
        except (OSError, HTTPException) as error:
            if attempt == ATTEMPTS:
                raise
            print("Transfer failed: {}. Retrying...".format(error), file=sys.stderr)
            time.sleep(attempt)


def _select(values, allowed, label):
    if values is None:
        selected = list(allowed)
    elif isinstance(values, str):
        selected = [values]
    else:
        selected = list(dict.fromkeys(values))
    if not selected:
        raise ValueError("Select at least one " + label + ".")
    unknown = [name for name in selected if name not in allowed]
    if unknown:
        raise ValueError("Unknown {}: {}".format(label, ", ".join(map(str, unknown))))
    return selected


def _download_zhuang(output_dir, samples, components, include_parcellation,
                     dry_run, legacy_layout=False, include_taxonomy=False):
    selected_samples = _select(samples, SAMPLES, "samples")
    selected_components = _select(components, COMPONENTS, "components")
    output = Path(output_dir).expanduser().resolve()
    result = {
        "samples": {
            sample: {component: [] for component in selected_components}
            for sample in selected_samples
        },
        "parcellation": None,
        "taxonomy": [],
    }
    files = []
    failures = []
    for sample in selected_samples:
        for component in selected_components:
            prefix = SAMPLE_PREFIXES[sample][component]
            print("Listing {} / {}...".format(sample, component), flush=True)
            try:
                component_files = list(list_files(prefix))
                if component == "expression":
                    component_files = [
                        (key, size) for key, size in component_files
                        if key.endswith("-raw.h5ad")
                    ]
                if not component_files:
                    label = "raw .h5ad expression files" if component == "expression" else "files"
                    raise ValueError("No {} found in {}".format(label, prefix))
                release = prefix.rstrip("/").rsplit("/", 1)[-1]
                for key, size in component_files:
                    if legacy_layout:
                        relative = Path(key)
                    else:
                        relative = (
                            Path(sample) / LOCAL_DIRS[component] / release
                            / key[len(prefix):]
                        )
                    files.append((sample, component, key, size, relative))
            except SOURCE_ERRORS as error:
                failures.append("{} / {}".format(sample, component))
                print("Could not list {} / {}: {}".format(sample, component, error),
                      file=sys.stderr)

    if include_parcellation:
        print("Listing separate parcellation table...", flush=True)
        try:
            matches = [item for item in list_files(PARCELLATION_KEY)
                       if item[0] == PARCELLATION_KEY]
            if len(matches) != 1:
                raise ValueError("Requested parcellation table was not found.")
            key, size = matches[0]
            relative = Path("parcellation") / key[len("metadata/"):]
            files.append((None, "parcellation", key, size, relative))
        except SOURCE_ERRORS as error:
            failures.append("parcellation")
            print("Could not list parcellation: {}".format(error), file=sys.stderr)

    if include_taxonomy:
        print("Listing WMB (mouse) taxonomy files...", flush=True)
        try:
            taxonomy_files = list(list_files(TAXONOMY_PREFIX))
            if not taxonomy_files:
                raise ValueError("No taxonomy files found in " + TAXONOMY_PREFIX)
            for key, size in taxonomy_files:
                relative = Path(TAXONOMY_DIR) / key[len(TAXONOMY_PREFIX):]
                files.append((None, "taxonomy", key, size, relative))
        except SOURCE_ERRORS as error:
            failures.append("taxonomy")
            print("Could not list taxonomy: {}".format(error), file=sys.stderr)

    print(
        "Found {} files, {:.2f} GB in total.".format(
            len(files), sum(item[3] for item in files) / 1e9
        ),
        flush=True,
    )
    for sample, component, key, size, relative in files:
        destination = output / relative
        if dry_run:
            print("{} bytes  {} -> {}".format(size, object_url(key), destination))
        else:
            try:
                download_file(key, size, output, relative_path=relative)
            except (OSError, HTTPException) as error:
                failures.append(key)
                print("Failed: {}: {}".format(key, error), file=sys.stderr)
                continue
        if component == "taxonomy":
            result["taxonomy"].append(destination)
        elif sample is None:
            result["parcellation"] = destination
        else:
            result["samples"][sample][component].append(destination)

    if failures:
        raise DownloadError(
            "Could not complete: {}. Rerun to retry.".format(", ".join(failures))
        )
    if dry_run:
        print("Dry run complete; no data files were downloaded.")
    else:
        print("Download complete: " + str(output))
    return result


def download_zhuang_abca(output_dir="data", samples=None, components=None,
                         include_parcellation=True, dry_run=False,
                         include_taxonomy=True):
    return _download_zhuang(
        output_dir, samples, components, include_parcellation, dry_run,
        include_taxonomy=include_taxonomy,
    )


def download_zhuang_abca2(output_dir="abc_atlas", components=None, dry_run=False):

    result = _download_zhuang(
        output_dir, ("Zhuang-ABCA-2",), components,
        include_parcellation=False, dry_run=dry_run, legacy_layout=True,
    )
    return result["samples"]["Zhuang-ABCA-2"]


def main(argv=None):
    """Command-line wrapper; the reusable function above never parses argv."""
    args = parse_args(argv)
    try:
        download_zhuang_abca(
            output_dir=args.output,
            samples=args.samples,
            components=args.components,
            include_parcellation=not args.skip_parcellation,
            dry_run=args.dry_run,
            include_taxonomy=not args.skip_taxonomy,
        )
    except DownloadError as error:
        print(str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted. Rerun the command to retry.", file=sys.stderr)
        sys.exit(130)
