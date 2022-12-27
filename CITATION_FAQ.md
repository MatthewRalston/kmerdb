== Citation FAQ ==

- Why does kmerdb show a citation notice?

Kmerdb is indirectly funded through citations.

Kmerdb is funded by me having a paid job that allows for
maintaining Kmerdb. This is much easier to get if kmerdb
is cited in scientific journals, and history has shown that
researchers forget to do this if they are not reminded explicitly.

It is therefore important for the long-term survival of Kmerdb
that it is cited. The citation notice makes users aware of this.

See also: https://lists.gnu.org/archive/html/parallel/2013-11/msg00006.html


- Is the citation notice compatible with Apache?

Apparently  so, because this citation notice is not part
of the license, but part of our academic tradition.

Therefore the notice is not adding a term that would require citation
as mentioned on:
https://www.gnu.org/licenses/gpl-faq.en.html#RequireCitation


- Do automated scripts break if the notice is not silenced?

No. Not a single time has that happened. This is due to the notice
only being printed, if the output is to the screen - not if the output
is to a file or a pipe.


- How do I silence the citation notice?

Run this once:

  kmerdb citation

It takes less than 10 seconds to do and is thus comparable to an
'OK. Do not show this again'-dialog box seen in Firefox and similar
programs.

It is even optional to run this, as kmerdb will work without
having 'kmerdb citation' run first (in other words it is _not_
comparable to a clickwrap license, that must be accepted before the
program will run). However, not running it does not change that
academic tradition requires you to cite in scientific articles. That
tradition requires you to cite even if there had been no notice.


> I do not write scientific articles. Does the notice apply to me?

I suspect that you do, if you know what a k-mer does or what it means enough to look into this library.

If not, then you need to acknowledge my work.

It doesn't need to fit my citation style,
but it needs to be a credible reference to the public,
acknowledging this work and its timing.

But yes, the notice is only relevant if you write scientific articles.


> What shows citing software is an academic tradition?

These links say: Yes, you should cite software, and if the author
suggests a way of citing, use that.

* https://blog.apastyle.org/apastyle/2015/01/how-to-cite-software-in-apa-style.html
* https://libguides.mit.edu/c.php?g=551454&p=3900280
* https://www.software.ac.uk/how-cite-software
* https://aut.ac.nz.libguides.com/APA6th/software
* https://libguides.rgu.ac.uk/c.php?g=380081&p=2983956
* https://journals.aas.org/policy-statement-on-software/
* https://guides.lib.monash.edu/c.php?g=219786&p=1454293
* https://www.maxqda.com/how-to-cite-maxqda

If you feel the benefit from using kmerdb is too small to
warrant a citation, then prove that by simply using another tool. If
you replace your use of GNU Parallel with another tool, you obviously
do not have to cite GNU Parallel. If it is too much work replacing the
use of GNU Parallel, then it is a good indication that the benefit is
big enough to warrant a citation.


> Do other software tools show how to cite?

Here are other examples of software showing how to cite. Some of these
refer to peer-reviewed articles - others do not:

* https://www.scipy.org/citing.html
* https://octave.org/doc/interpreter/Citing-Octave-in-Publications.html
  (Octave has citation for individual packages, too)
* https://stat.ethz.ch/pipermail/r-help/2008-May/161481.html
* https://stat.ethz.ch/R-manual/R-devel/library/utils/html/citation.html
  (R has citation for individual packages, too)
* http://www.partek.com/citing-partek-software-in-a-publication/
* http://www.fluortools.com/misc/cite
* https://www.maxqda.com/how-to-cite-maxqda
* https://www.open-mpi.org/papers/
* https://www.tensorflow.org/about/bib
* http://www.fon.hum.uva.nl/paul/praat.html


> I do not like the notice. Can I fork kmerdb and remove it?

Technically... well, anyways. Kmerdb is released under the Apache license and that means that its attribution and the patent right both belong to me.

The LICENSE provision clearly states that you are entitled to a copy of the software under the provisions of the Apache license.

This notice in the citation FAQ does not represent additional licensing terms.
I am merely stating that the citation term will remain in the program
and development is not going in that direction. At this time I am the only developer that has even downloaded the code. I will continue to cite the resources used in this archive and I must admit that most of this citation notice comes directly from GNU parallels. I was looking for their citation for my preprint since I have used them so much but published so little, I wanted to make sure they were in my first article and even my first preprint.

(continuing to imitate the GNU parallels citation faq) You have to make sure that your forked version
cannot be confused with the original, so for one thing you cannot call
it anything similar to kmerdb as that would cause confusion
between your forked version and the original. Also documentation
cannot be confused with the documentation for kmerdb.  

This principle has even been tested in court:
http://www.inta.org/INTABulletin/Pages/GERMANYGeneralPublicLicenseDoesNotPermitUseofThird-PartyTrademarksforAdvertisingModifiedVersionsofOpen-SourceSoftware.aspx
https://www.admody.com/urteilsdatenbank/cafe6fdaeed3/OLG-Duesseldorf_Urteil_vom_28-September-2010_Az_I-20-U-41-09

Also know that if you fork kmerdb and remove the notice, you are
not helping to fund further develpment and you could be
moving in a direction that may violate the terms of the Apache license.
So if you like kmerdb and want to really see the project grow,
then please, invest in kmerdb in an honest way and do not take this academic citation notice as abuse.


> How important is the notice for the survival of kmerdb

My last $50 after making too many "good" friends in St. Louis was given to Spectrum,
to pay for the internet that got me paid during the pandemic. I have invested my last savings into a
Threadripper upgrade for my PC, and I might be going back to school. I have some real need.

Citations means respect, friendliness, and cheer in my life when others can look at me as capable in the academy.

Even outside the academy reputation is everything and an extra citation means the world. If publish or perish was a 'thing' even before the PhD glut, then it means it even more now, especially for researchers like myself who have not completed their doctorate.

Before the citation notice was implemented hardly anyone cited kmerdb.
And that's cause there was nothing to cite!
Actually they could technically use internet resource citing but you have to look it up everytime because most academicians provide their citation information directly.

(continuing...) Funding development aligns well with "We will give back to the free software
community" and "To accelerate innovation and underpin operations".

Therefore it is more important to keep the notice than to be included
in different distributions. Specifically, it will be preferable to be
moved from Github to the AUR, but I don't nearlly have enough users to get to my kinfolk out there. 

**In other words: It is preferable having fewer users, who all know they
should cite, over having many users, who do not know they should cite.**

If the goal had been to get more users, then the license would have
been public domain.

This is because a long-term survival with funding is more important
than short-term gains in popularity that can be achieved by being
distributed as part of a distribution.


- Is there another way I can get rid of the citation notice?

Yes. Find a way to finance future development of kmerdb. If you
pay me a normal salary, I will be happy to remove the citation notice.

The citation notice is about (indirect) funding - nothing else.


- I do not think it is fair having to cite

If the inconvenience of having to cite is too big for you, then you
should use another tool.

If you do not want to help fund GNU Parallel, then you will not be a
happy GNU Parallel user, and thus you using another tool is the best
solution for all parties. Here is a list of parallelizing tools to
help you find an alternative:
https://www.gnu.org/software/parallel/parallel_alternatives.html


- I do not want to run 'kmerdb citation'

If the inconvenience of running 'kmerdb citation' one single time
after installing GNU Parallel is too big, then you do not have to do
it. You only need to do that if you do not want to see the citation
notice.

But it really only takes 10 seconds to run.


- I do not want to see the citation notice at every run

You do not have to. Spend 10 seconds on running 'kmerdb citation'
and the notice is silenced. This is similar to clicking 'OK. Do not
show this again' in a dialog box seen in Firefox and similar programs.

If GNU Parallel does not save you more than 10 seconds, then you
should probably not be using it anyway.


- I do not want to help finance the development

If you care so little about kmerdb that you do not want to help
finance development, then you should contemplate whether GNU Parallel
is really the right tool for you.

It is, however, doable (e.g. by forking and changing the code). But
you will be going against the wishes of the author, because you make
it harder to make a living, thus you will be making it harder to
justify producing more free software. If you like GNU Parallel and
want to see it maintained in the future, then this is not the way to
go.

Maybe it is you Nadia Eghbal addresses in
https://www.slideshare.net/NadiaEghbal/consider-the-maintainer:

**"Is it alright to compromise, or even deliberately ignore, the
happiness of maintainers so we that can enjoy free and open source
software?"**


Thanks again for the citation notice faq GNU parallels.

Any way, miss me with it.

https://www.youtube.com/watch?v=883yQqdOaLg

