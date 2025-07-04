# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!--
The rules for this file:
  * entries are sorted newest-first.
  * summarize sets of changes - don't reproduce every git log comment here.
  * don't ever delete anything.
  * keep the format consistent:
    * do not use tabs but use spaces for formatting
    * 79 char width
    * YYYY-MM-DD date format (following ISO 8601)
  * accompany each entry with github issue/PR number (Issue #xyz)
-->

## [1.1.0] - 2025-07-04

### Authors
* @rsexton2
* @copilot
* @orbeckst

### Added
* Support for combining contact timeseries from multiple repeat runs through
  new `CombineContacts` class and `python -m basicrta.combine` CLI interface.
  Enables pooled analysis of binding kinetics data with metadata preservation
  and trajectory source tracking (Issue #16)

### Changed
* package has final paper citation



## [1.0.0] - 2025-05-24

### Authors
* @rsexton2
* @ianmkenney
* @rjoshi44
* @orbeckst

### Added
* added option processing for label-cutoff to cluster.py (PR #13)

### Fixed
* Fix package detection and installation (PR #12)
* fix citation in reST docs (PR #7)
* update codcov action in workflow (PR #9)

### Removed
* no testing on Windows, temporarily exclude windows-latest from CI (PR #11)

## [0.2.0] - 2024-11-14

### Authors
* @rsexton2

### Summary
Feature-complete release.

### Added
* Workflow executable through command-line
* updated docs/tutorial

