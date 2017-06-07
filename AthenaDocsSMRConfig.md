---
title: Configuring and Compiling with SMR
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/SMR Configuring and Compiling

To configure Athena to include SMR, simply enable the SMR feature:

	% configure --enable-smr

Then compile the code as normal:

	% make all

No special libraries are needed, so no special links are required.  

SMR may be combined with almost any of the other
features and packages in Athena (see the [User Guide]({{site.baseurl}}/AthenaDocsUG) for a complete list of configurable options).
However, be warned, since SMR is under continual development, **not all features and/or packages will work properly with SMR**.
