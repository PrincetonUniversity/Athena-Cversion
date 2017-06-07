---
title: The <Job> Block
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Job Block

Parameters in this block control properties of the jobs.  Example:


	<job>
	problem_id      = Brio-Wu    # problem ID: basename of output filenames
	maxout          = 3          # Output blocks number from 1 -> maxout
	num_domains     = 1          # number of Domains in Mesh


**problem_id:**
string added as basename of output filenames, see [Overview of Outputs]({{site.baseurl}}/AthenaDocsUGOutputOverview).  Usually same
as *problem-name* used in input file.  There is no maximum length
for this name (except that set by the maximum length of a
line in the input file, which is 256 characters).

**maxout:** Specifies how many output blocks will be read from
input file.  Output
blocks `<output1>` through `<outputN>` where `N = maxout`
are scanned for valid output descriptions.  Missing output
blocks are permitted.  Output blocks with `N>maxout` are ignored.

**num_domains:**  Specifies how many grid domains will be read from the input file.
Domian blocks `<domain1>` through `<domainN>` where `N=num_domains`
will be read.  Domain blocks with `N>num_domains` will be ignored.
