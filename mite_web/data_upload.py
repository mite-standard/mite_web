"""CLI-interface to upload submitted data as Issues on `mite_data`, intended to be used by mite-bot.

This script needs a functioning and authenticated installation of the GitHub CLI.

Copyright (c) 2024-present Mitja Maximilian Zdouc, PhD

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import requests
from jsondiff import diff
from pydantic import BaseModel, DirectoryPath


class IssueManager(BaseModel):
    """Methods to upload submitted data as issues to GitHub

    Attributes:
        version: the version of mite_data that was downloaded
        record_url: pointing to the newest MITE Zenodo entry
        src: the location of the data
        mite_current_zip: location to download current mite_data
        mite_current_unzip: location of the unzipped mite_data
        mite_current_data: location of the current mite_data json files
        reviewers: the people to assign issues to
    """

    version: str | None = None
    record_url: str = "https://zenodo.org/api/records/13294303"
    mite_current_zip: DirectoryPath | None = None
    mite_current_unzip: DirectoryPath | None = None
    mite_current_data: DirectoryPath | None = None
    src: DirectoryPath
    reviewers: list = ["@mmzdouc", "@jorgecnavarrom" "@adafede"]

    def run(self):
        """Parses files in the source directory and pushes them to GitHub"""

        self.download_current_data()

        for filename in self.src.iterdir():
            if filename.suffix != ".json":
                continue

            with open(filename) as infile:
                json_dict = json.load(infile)

            title = ""
            if json_dict.get("accession") == "MITE9999999":
                title = f"NEW entry {json_dict.get('enzyme', {}).get('name')} ({filename.stem})"
                self.body_new(data=json_dict)
            else:
                title = f"EXISTING entry {json_dict.get("accession")} ({filename.stem})"
                self.body_existing(data=json_dict)

            subprocess.run(
                [
                    "gh",
                    "issue",
                    "create",
                    "--repo",
                    "mite-standard/mite_data",
                    "--title",
                    f"{title}",
                    "--body-file",
                    f"{self.src.joinpath("temp.txt")}",
                    "--label",
                    "review",
                ],
                check=True,
            )

        self.cleanup_current_data()

    def download_current_data(self) -> None:
        """Downloads the current MITE dataset"""
        self.mite_current_zip = self.src.joinpath("record.zip")
        self.mite_current_unzip = self.src.joinpath("record")

        response_metadata = requests.get(self.record_url)
        if response_metadata.status_code != 200:
            raise RuntimeError(
                f"Error fetching 'mite_data' record metadata: {response_metadata.status_code}"
            )

        record_metadata = response_metadata.json()
        self.version = record_metadata["metadata"]["version"]
        files_url = record_metadata["files"][0]["links"]["self"]

        response_data = requests.get(files_url)
        if response_data.status_code != 200:
            raise RuntimeError(
                f"Error downloading 'mite_data' record: {response_data.status_code}"
            )

        with open(self.mite_current_zip, "wb") as f:
            f.write(response_data.content)

        shutil.unpack_archive(
            filename=self.mite_current_zip,
            extract_dir=self.mite_current_unzip,
            format="zip",
        )

        if not self.mite_current_unzip.exists():
            raise NotADirectoryError(
                f"Could not find the unzipped directory {self.mite_current_unzip}."
            )

        matching_dirs = list(self.mite_current_unzip.glob("mite-standard-mite_data-*"))
        if not matching_dirs:
            raise RuntimeError(
                f"Could not determine data storage location in downloaded directory."
            )

        subdir = matching_dirs[0]
        self.mite_current_data = self.mite_current_unzip.joinpath(subdir).joinpath(
            "mite_data/data"
        )

    def cleanup_current_data(self) -> None:
        """Removes the MITE dataset again"""
        os.remove(self.mite_current_zip)
        shutil.rmtree(self.mite_current_unzip)

    def body_new(self, data: dict) -> None:
        """Creates issue body text for a new entry.

        Arguments:
            data: the newly submitted data
        """
        body = f"""
A submission for a new entry was performed via the MITE Submission portal.
Please find the file in the code block below.

### Labels
https://github.com/mite-standard/mite_data/labels/review

### Review requested
{", ".join(self.reviewers)}

## TODO Reviewers

See the (Reviewer Instructions)[https://github.com/mite-standard/mite_data/wiki/How-to-Review-Entries] in the Wiki.

- [x] Automated validation using `mite_extras` performed
- [ ] Factual correctness of submission checked
- [ ] Database ID crosslinks verified, added where possible
- [ ] References properly formatted
- [ ] Automated validation check ID "BBBBBB..." replaced with real reviewer ID

Please propose and discuss changes in the comments.

## Submitted Data

```
{json.dumps(data, indent=4)}
```

*This action was performed by `mite-bot`*
"""

        with open(self.src.joinpath("temp.txt"), "w") as outfile:
            outfile.write(body)

    def body_existing(self, data: dict) -> None:
        """Creates issue body text for an existing entry.

        Arguments:
            data: the newly submitted data
        """
        with open(
            self.mite_current_data.joinpath(f'{data.get("accession")}.json')
        ) as infile:
            original_data = json.load(infile)

        diff_string = diff(a=original_data, b=data, marshal=True)

        body = f"""
A submission for an existing entry was performed via the MITE Submission portal.
Please find the file in the code block below.

### Labels
https://github.com/mite-standard/mite_data/labels/review

### Review requested
{", ".join(self.reviewers)}

## TODO Reviewers

See the (Reviewer Instructions)[https://github.com/mite-standard/mite_data/wiki/How-to-Review-Entries] in the Wiki.

- [x] Automated validation using `mite_extras` performed
- [ ] Factual correctness of submission checked
- [ ] Database ID crosslinks verified, added where possible
- [ ] References properly formatted
- [ ] Automated validation check ID "BBBBBB..." replaced with real reviewer ID

Please propose and discuss changes in the comments.

## Proposed changes (compared to current version in mite_data v{self.version})

```
{json.dumps(diff_string, indent=4)}
```

## Submitted Data

```
{json.dumps(data, indent=4)}
```

*This action was performed by `mite-bot`*
"""

        with open(self.src.joinpath("temp.txt"), "w") as outfile:
            outfile.write(body)


def setup_cli(args: list) -> argparse.Namespace:
    """Define command line interface options.

    Arguments:
        args: a list of CLI arguments

    Returns:
        argparse.Namespace object with command line parameters
    """
    parser = argparse.ArgumentParser(
        description=f"CLI to upload submitted data to mite_data.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input_dir",
        type=str,
        required=True,
        help="Specifies a directory containing submitted files for processing.",
    )

    return parser.parse_args(args)


def main() -> None | SystemExit:
    """Run main part of code"""
    args = setup_cli(sys.argv[1:])

    try:
        issue_manager = IssueManager(src=Path(args.input_dir).resolve())
        issue_manager.run()
    except Exception as e:
        print(f"An exception has occurred: {e}")
        return sys.exit(1)


if __name__ == "__main__":
    main()
