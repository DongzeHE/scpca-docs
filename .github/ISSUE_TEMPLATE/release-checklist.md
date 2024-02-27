---
name: Release checklist
about: Prepare for a new release of scpca-docs
title: Prepare for scpca-docs release
labels: release
assignees: ''

---

## Steps for a new release of `scpca-docs`
<!-- Please include any notes about what will be included in this release -->
<!-- e.g., you might want to mention if the release is related to a specific portal feature, like addition of AnnData objects -->
<!-- Update the title of this issue to reflect the contents of the release (e.g., Prepare for scpca-docs release: AnnData) -->

### Preparing for the release
- [ ] Does this release require an update to the CHANGELOG because there will be changes to downloadable files? If so, please create an issue tracking the CHANGELOG update and mark it as blocking this issue.
- [ ] Are all other issues planned for this release resolved? If any issues are unresolved, mark this issue as blocked by those on ZenHub.
- [ ] Optional: If not all changes in `development` are ready to be released, create a feature branch off `main` and cherry pick commits from `development` to that feature branch.
- [ ] File a PR from the `development` branch (or the new feature branch) to the `main` branch. This should include all of the changes that will be associated with the next release.
  - [ ] If a CHANGELOG entry was required, add the date to the entry's header as part of this PR.

### Creating a release
- [ ] On the [releases page](https://github.com/AlexsLemonade/scpca-docs/releases), choose `Draft a new release`.
- [ ] In `Choose a tag`, use the current date as the release tag (`YYYY.MM.DD`), then click `Create a new tag: YYYY.MM.DD on publish`.
- [ ] Write a description of the major changes in this release. You may want to start with the auto-generated release notes to save time.
- [ ] Optional: If not all issues have been addressed, save a draft to return to later.
- [ ] Publish the release!
