---
title: Code Style
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Code Style

The GNU development tools (gmake, gcc, autoconf) have been used
to develop the code.  It *should* work on any platform that properly supports
these tools (undoubtedly, however, there will be some platforms on which it will fail).
To date we have tested the code on Linux, Solaris, and MacOS X.

Coding style is very personal and has been the source of many animated
discussions at Cafe Nero and Small World Coffee.  We have only agreed
on the following:

 * indent by 2 as in the K&R style.
 * try to limit line widths to 80 characters, to prevent line wraps in windows of that width.
