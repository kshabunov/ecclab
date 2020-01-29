# ECCLab

This is a new repo for the set of source files that I have been developing and used for
error correcting code research with [Dr. Ilya Dumer](https://www.itsoc.org/profiles/ilyadumer) since 1998.

Initially they were mostly "research" and "proof of concept" quality for private use.
Lately, we decided to revamp them, mostly for better readability, and put here for reference
and future development. Currently, this is work in progress.

## Build
Run `make` in the project root directory. It will create `./build` subdirectory
(if it does not exist yet) and put executables there.

## Running programs

### Typical usage

The simulation programs are run from the command line:
```
executable simulation_parameters_file [options]
```
Common options are:
* `-nr` - don't randomize.
* `-si` - set saving interval in seconds.
* `-ri` - set return (status update) interval in seconds.

The simulation parameters files are plain text with the following structure:
* Empty lines are skipped.
* Lines starting with `#` are considered to be comments and ignored.
* Each parameter is a key/value pair separated by one or more spaces.
* Parameter may have a group value surrounded by curly brackets (`{}`).

For example:
```
% Include file:
include file.txt

% Single value parameter:
param value

% Single boolean value parameter (can be "on" or "off):
boolean_param on

% Group value parameter. Items may have several tokens separated by spaces:
group_param {
  item1
  item2
}
```

### Common simulation parameters

* `res_file filename` - file for saving simulation results. String. Required.
* `SNR_val_trn {...}` - SNR values to simulate at. Group value,
 each item contains SNR value and the number of trials (see an example below). Required.
* `random_codeword on|off` - Use random information sequence for every simulation trial or not.
 Boolean. Off by default.
* `ml_lb on|off` - Estimate ML lower bound or not. Boolean. Off by default.

Example:
```
res_file result.spf
SNR_val_trn {
   1      50000
   1.25   50000
   1.5    50000
}
random_codeword on
ml_lb off
```

### `dtrm_glp_gc`
Simulate recursive list decoding for RM/SubRM codes, AWGN channel.
This algorithm was applied later to Polar codes with much success
and now is known as Successive Cancellation List (SCL) decoding.
For pure RM codes, we may also use up to mCr permutations.

Specific simulation parameters:
* `RM_m m` - RM m parameter. Integer, positive, non-zero. Required.
* `RM_r r` - RM r parameter. Integer, positive, non-zero. Required.
* `order_node_mask bin_str` - mask marking with zeroes "border" nodes that should be set as zeroes
 thus specifying the RM subcode (better explanation TBD later). Binary string. Required.
* `border_node_lsize L` - list size. Integer, positive, non-zero. Required.

## References
* I. Dumer, K. Shabunov, "Soft-decision decoding of ReedMuller codes: Recursive lists",
 IEEE Transactions On Information Theory, vol. 52, no. 3, March 2006. [IEEE Xplore](https://ieeexplore.ieee.org/document/1603792)
* I. Dumer, K. Shabunov, "Recursive and permutation decoding for Reed-Muller codes",
 Proc. 2002 IEEE Int. Symp. Inform. Theory, pp. 146, 2002. [IEEE Xplore](https://ieeexplore.ieee.org/document/1023418)
* I. Dumer, K. Shabunov, "Near-optimum decoding for subcodes of Reed-Muller codes",
 2001 IEEE Intern. Symp. Info. Theory, pp. 329, June 24–29, 2001. [IEEE Xplore](https://ieeexplore.ieee.org/document/936192/)
* I. Dumer, K. Shabunov, "Recursive decoding of Reed-Muller codes",
 Proc. 2000 IEEE Int. Symp. Inform. Theory, pp. 63, 2000. [IEEE Xplore](https://ieeexplore.ieee.org/document/866353)