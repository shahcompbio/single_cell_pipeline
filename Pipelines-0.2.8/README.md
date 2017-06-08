# Overview

Pipelines is a set of helper scripts for working with the ruffus pipeline library. It adds some utility classes for making directories, running jobs on clusters and handling transfers file transfers via ssh. 

# License

Pipelines is licensed under the GPL v3, see the LICENSE.txt file for details.

# Versions

## 0.2.8

- Updated ClusterJobManger to support log file cleanup and output paths to log files.

## 0.2.7

- Added additional function in the utils for parsing file names as key/value pairs.

- Reduced the default memory requested by the cluster job manager.

## 0.2.6

- Removed the ability to launch asynchronous jobs and wait with job managers

- Added exception handling for DRMAA communicaton errors due to unreliable qmasters

## 0.2.5.2

- Modified the cluster job manager to not use parallel environment when only using 1 core

## 0.2.5

- Added some helper functions to cut down on boiler plate in ruffus pipelines.

### 0.2.5.1

- Bug fix in ruffus_helpers

## 0.2.4

- Updated cluster job manager to output information about commands which fail along with job id

## 0.2.3

- Allowed for specification of file permissions for io functions

## 0.2.2

- Add the ability to pass job names to ClusterJobManager

- Added the ability to transfer files from local host to remote host

## Dependencies

### Required

* None

### Optional

- [drmaa](https://pypi.python.org/pypi/drmaa/0.7.6) >= 0.7.6 - Required if using ClusterJobManger

- [spur](https://pypi.python.org/pypi/spur/0.3.7) >= 0.3.7 - Required if using RemoteJobManager

- [ruffus](https://pypi.python.org/pypi/ruffus/2.4) >= 2.4 - Required if using ruffus module
