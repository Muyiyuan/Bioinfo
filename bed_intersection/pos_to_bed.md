#### Description

<p>Generate bed file based on pos file</p>

#### Usage

```shell
usage: pos_to_bed.py [-h] -p POS -b BED [-e {True,False}]

generate bed file according pos file

optional arguments:
  -h, --help            show this help message and exit
  -p POS, --pos POS     pos file
  -b BED, --bed BED     bed name
  -e {True,False}, --equal {True,False}
                        start = end
```

#### Note
1. Pos file format (Tab-delimited)
> chr1  1
>
> chr1  2
>
> ... ...
>
> chrx  10
2. Argument "-e"
<p>The dufault vaule of "-e": "True".</p>
<p>If you do not want to save the bed sub region that start postion is equal to end position, you can use "-e False"</p>
