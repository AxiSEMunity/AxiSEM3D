# Running AxiSEM3D


AxiSEM3D reads a set of input files from an input directory and writes output to a separate directory. 

* `-- input <dir>`  Input directory: default ./input
* `-- output <dir>` Output directory: default ./output

To see the options available run:

```bash
axisem3d --help
```

AxiSEM3D expects the following input files:

<span style="margin-left:48px;">inparam.model.yaml<br></span>
<span style="margin-left:48px;">inparam.nr.yaml<br></span>
<span style="margin-left:48px;">inparam.source.yaml<br></span>
<span style="margin-left:48px;">inparam.output.yaml<br></span>
<span style="margin-left:48px;"> inparam.advanced.yaml<br></span>
<span style="margin-left:48px;"> and/or a user specified station file such as
GSN.txt or USArray.txt<br></span>

 .yaml files are human readable.


---

```{toctree}
---
maxdepth: 1
---
example.md
y_model.md
y_source.md
y_nr.md
y_output.md
y_advanced.md
```
