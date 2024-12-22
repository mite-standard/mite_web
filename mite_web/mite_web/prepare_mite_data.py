"""Populates the app with data from mite_data

This script runs as part of the Dockerfile during container-building.

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

import sys

from mite_web.prep_data.aux_file_manager import AuxFileManager, SummaryManager
from mite_web.prep_data.download_manager import DownloadManager
from mite_web.prep_data.html_json_manager import HtmlJsonManager


def main() -> None | SystemExit:
    """Prepares mite_data for use in mite_web"""
    try:
        download_manager = DownloadManager()
        download_manager.run()

        json_manager = HtmlJsonManager()
        json_manager.run()

        summary_manager = SummaryManager()
        summary_manager.run()

        aux_manager = AuxFileManager()
        aux_manager.run()
    except Exception as e:
        print(f"An exception has occurred: {e}")
        return sys.exit(1)


if __name__ == "__main__":
    main()
