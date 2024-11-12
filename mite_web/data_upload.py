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
import subprocess
import sys
from pathlib import Path

from pydantic import BaseModel, DirectoryPath


class IssueManager(BaseModel):
    """Methods to upload submitted data as issues to GitHub

    Attributes:
        src: the location of the data
        reviewers: the people to assign issues to
    """

    src: DirectoryPath
    reviewers: list = ["@mmzdouc", "@adafede"]

    def run(self):
        """Parses files in the source directory and pushes them to GitHub"""

        for filename in self.src.iterdir():
            with open(filename) as infile:
                json_dict = json.load(infile)

            title = ""
            if json_dict.get("accession") == "MITE9999999":
                title = f"NEW entry ({filename.stem})"
            else:
                title = f"EXISTING entry {json_dict.get("accession")} ({filename.stem})"

            body = f"""
A submission was performed via the MITE Submission portal.
Please find the file in the code block below.

### Labels
https://github.com/mite-standard/mite_data/labels/review

### Review requested
{", ".join(self.reviewers)}

## TODO Reviewers

- [x] Automated validation using `mite_extras` performed
- [ ] Factual correctness of submission checked
- [ ] Database ID crosslinks verified, added where possible
- [ ] References properly formatted
- [ ] Automated validation check ID "BBBBBB..." replaced with real reviewer ID

Please propose and discuss changes in the comments.

## Submitted Data

```
{json.dumps(json_dict, indent=4)}
```

*This action was performed by `mite-bot`*
"""
            with open(self.src.joinpath("temp.txt"), "w") as outfile:
                outfile.write(body)

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
