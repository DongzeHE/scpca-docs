# scpca-docs

This repository contains user facing documentation for the ScPCA Portal.

## Spell Check

This repo contains a github workflow that performs spell check on all markdown files.
A list of words to be ignored by spell check is listed in `components/dictionary.txt` and all words are to be listed in alphabetical order.

## Running Locally

### To install:

#### First download the repo
```shell
 $ git clone https://github.com/AlexsLemonade/scpca-docs
 $ cd scpca-docs
```

#### Set up `virtualenv`
```shell
 $ virtualenv env
 $ source env/bin/activate
 $ pip install -r requirements.txt
```

### To run locally:
You must set up `virtualenv` before running locally.
To run a local server that watches for changes please run:

```shell
 $ ./dev-server.py
```
