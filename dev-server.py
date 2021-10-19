#!/usr/bin/env python3
import os
import sys

try:
    from livereload import Server, shell
except Exception:
    print(
        """Livereaload not found.
Please make sure your virtual env is set up correctly.
From the root directory of this repository:"""
    )
    print("\t$ virtualenv env")
    print("\t$ source env/bin/activate")
    print("\t$ pip install -r requirements.txt")
    sys.exit(1)

make_clean = "make -C docs clean"
make_html = "make -C docs html"

if __name__ == "__main__":
    os.system(make_clean)
    os.system(make_html)
    server = Server()
    server.watch("docs/*.rst", shell(make_html), delay=1)
    server.watch("docs/*.md", shell(make_html), delay=1)
    server.watch("docs/*.py", shell(make_html), delay=1)
    server.watch("docs/_static/*", shell(make_html), delay=1)
    server.watch("docs/_templates/*", shell(make_html), delay=1)
    server.serve(root="docs/_build/html")
