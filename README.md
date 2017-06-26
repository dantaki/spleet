# spleet
Validate SV with Split-Reads

## Install

#### Clone from GitHub

```
git clone --recursive https://github.com/dantaki/spleet.git
```

#### Compile with [CMake](https://cmake.org/)

```
cd spleet/ && mkdir build && cd build && cmake .. && make
```

##### Binary executable found under `spleet/build/src/spleet`

## Usage 

Spleet currently supports validation of deletion and duplication SV

`spleet --help`

```
eeeee   eeeee   e         eeee   eeee   eeeee
8   "   8   8   8         8      8        8
8eeee   8eee8   8e        8eee   8eee     8e
   88   88      88        88     88       88
8ee88   88      88eee     88ee   88ee     88

spleet         Validate SV with Split-Reads
Version: 1.0	Author: Danny Antaki <dantaki@ucsd.edu>

Usage: spleet -i <in.bam> -r <sv.bed> -x <FLOAT> -q <INT> -o <output.txt>

Options:
    -i        Input: BAM filename
    -r        Input: SV bed file
    -x        Minimum reciprocal overlap [0.5]
    -q        Mapping quality threshold [10]
    -o        Output: filename
```
## Input

* `-i    Input: BAM filename`
   *  BAM file must be indexed with `samtools index`

* `-r    Input: SV bed file`
   * The first four tab-delimited columns must be `CHROM START END [DEL|DUP]`   

## Output

| CHR | START | END | LENGTH | OVERLAP | READNAME | ALIGNMENT | STRANDS | SV | TYPE |
| --- | ----- | --- | ------ | ------- | -------- | ----------| ------- | --- | --- | 
| 22 | 582 | 777 | 196 | 0.999992 | SPLEET-READ | PRIMARY | +\|+ | 22:581-776 | DEL | 

## Acknowledgements

spleet uses [BamTools](https://github.com/pezmaster31/bamtools)

> *BamTools*
> Copyright (c) 2009-2010 Derek Barnett, Erik Garrison, Gabor Marth, Michael Stromberg

Notice: BamTools requires [zlib](http://zlib.net/).

## License

MIT License

Copyright (c) 2017 Danny Antaki

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

## Contact

dantaki@ucsd.edu
