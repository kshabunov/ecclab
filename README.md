# ECCLab

This is a new repo for the set of source files that I have been developing and used for
error correcting code research with [Dr. Ilya Dumer](https://www.itsoc.org/profiles/ilyadumer) since 1998.
The main area of the research addressed efficient recursive decoding for Reed-Muller (RM) codes
and their subcodes with frozen bits.
 
Originally, these were research programs for private use. Lately, we decided to revamp them,
mostly for a better readability. These programs can be used to run simulations of RM codes
and their subcodes on the binary symmetric channels and AWGN channels.
The programs can also be used for reference and developments. This work is in  progress.

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

### Simulation parameters file format

These are plain text with the following structure:
* Empty lines are skipped.
* Lines starting with `#` are considered to be comments and ignored.
* Each parameter is a key/value pair separated by one or more spaces.
* Parameter may have a group value surrounded by curly brackets (`{}`).

Parameter files may include other parameter files using `include some.file` directive.

For example:
```
% Include file:
include another.file

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

### `dtrm0_bg`
Simulate recursive decoding for RM codes, AWGN channel. RM code is recursively decomposed
down to (0, m) and (m, m) nodes.

Specific simulation parameters:
* `RM_m m` - RM m parameter. Integer, positive, non-zero. Required.
* `RM_r r` - RM r parameter. Integer, positive, non-zero. Required.

### `dtrm_glp_bg`
Simulate recursive list decoding for RM codes and their subcodes (SubRM), AWGN channel.
This algorithm was applied later to Polar codes with much success
and now is known as Successive Cancellation List (SCL) decoding.
For pure RM codes, we may also use up to mCr permutations.

Specific simulation parameters:
* `RM_m m` - RM m parameter. Integer, positive, non-zero. Required.
* `RM_r r` - RM r parameter. Integer, positive, non-zero. Required.
* `order_node_mask mask_str` - mask marking border nodes that should be set as zeroes
 thus specifying the RM subcode (see explanation below). Binary string. Required.
* `border_node_lsize L` - list size. Integer, positive, non-zero. Required.
* `permutations {...}` - permutations set. Group value, each line
 defines one permutation (see example below). No additional permutations,
 except (0 1 2 ...), are used by default.

`order_node_mask mask_str` is a binary string
with each bit symbol corresponding to a single (0, m) or (m, m) node on the border of
the RM code decomposition triangle. The order of the bits is the order by which the recursive
decoder visits these nodes. If the bit is set as `0` the corresponding node is excluded.
For example, here is the mask for RM(2, 4) with
the first border node (which is the leftmost (0, 2)) eliminated:
```
border_node_mask 011111
```
This results in a subcode of RM(2, 4) with k = 10.

`permutations` is a group parameter. Each line consists of m integers from 0 to m - 1 separated
by spaces. They specify the order of the boolean polynomial variables that define
the codewords (relative to the original order). We can use the permutations that put any r variables
to the last r positions because they result in permutations of the information symbols.
For example, here are all 6 possible permutations for RM(2, 4):
```
permutations {
  2 3 0 1
  1 3 0 2
  1 2 0 3
  0 3 1 2
  0 2 1 3
  0 1 2 3
}
```

## References
* I. Dumer, K. Shabunov, "Soft-decision decoding of Reed-Muller codes: Recursive lists",
 IEEE Transactions On Information Theory, vol. 52, no. 3, March 2006. [IEEE Xplore](https://ieeexplore.ieee.org/document/1603792)
* I. Dumer, K. Shabunov, "Recursive and permutation decoding for Reed-Muller codes",
 Proc. 2002 IEEE Int. Symp. Inform. Theory, pp. 146, 2002. [IEEE Xplore](https://ieeexplore.ieee.org/document/1023418)
* I. Dumer, K. Shabunov, "Near-optimum decoding for subcodes of Reed-Muller codes",
 2001 IEEE Intern. Symp. Info. Theory, pp. 329, June 24–29, 2001. [IEEE Xplore](https://ieeexplore.ieee.org/document/936192/)
* I. Dumer, K. Shabunov, "Recursive decoding of Reed-Muller codes",
 Proc. 2000 IEEE Int. Symp. Inform. Theory, pp. 63, 2000. [IEEE Xplore](https://ieeexplore.ieee.org/document/866353)