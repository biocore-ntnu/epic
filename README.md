# epic

diffuse domain chip seq caller based on SICER.

This is a pre-alpha release. Please do aggressively report issues, quirks, complaints and anything that just feels slightly off to the issue tracker. Also please ask questions and make docrequests - there are loads of neat stuff I have not documented.

I have been using epic with great success, but no others have looked at it so you are likely to encounter some issues.

#### License

MIT

#### Requirements

Python data science stack and a fairly recent version of Pandas (0.17 =<).
Python 2.7 or 3+

#### Quickstart

```
pip install bioepic
# you only need git clone to get the test data
git clone https://github.com/endrebak/epic.git
epic -i control epic/examples/test.bed epic/examples/control.bed
```

#### Difference from the original SICER

Note that the island enriched threshold computation produces results that are < ~10% more conservative than in the original SICER.
This is due to numerics (summing many very small numbers is done in both implementations, albeit slightly differently).

This gives a different cutoff than the original, but produces virtually identical results (since epic and SICER produces the same candidate island list, with the same order, but epic selects slightly fewer islands from this list).

#### Why another ChIP-Seq domain caller?

MACS2 is great for narrow peaks, but epic performs better on diffuse domains. For medium size domains, such as PolII, our tests indicate that both perform about equally well, but epic uses only a fraction of the time.

#### Why not SICER?

SICER is a wonderful piece of software, but advances in the Python data science libraries has made it possible to implement it much more efficiently. Furthermore, SICER was not made to handle the mountains of data we have now; it simply cannot run on very large datasets due to (sensible) restrictions in the original implementation.

#### Credit

Chongzhi Zang, Dustin E. Schones, Chen Zeng, Kairong Cui, Keji Zhao and Weiqun Peng for the original SICER. Please consider citing their paper (*in addition* to our pre-print) if you use epic. And if you use any (helper) scripts in SICER that are not included in epic you should of course cite the SICER paper!

Most of the improvements in epic were possible due to Python Science libraries that were not available when SICER was originally written. Thanks to the Pandas developers!

Endre Bakken Stovner for the implementation of epic.

#### Why the name epic?

It stands for electronic pic [sic] caller or epigenome cartographer, whichever you prefer. Or perhaps it isn't just another bogus bioinformatics acronym. Hope you find the name fitting.

#### Other pieces of software you might prefer

* [SICER](http://home.gwu.edu/~wpeng/Software.htm) - great diffuse domain ChIP-Seq caller (which epic is based on.)
* [SICERpy](https://github.com/dariober/SICERpy) - a wrapper around SICER for convenience/parallelism.
* [csaw](https://github.com/LTLA/csaw) - R package. Uses an approach to island finding that complements epic very well. Requires more statistical sophistication and programming skill to use.
* [MACS2](https://github.com/taoliu/MACS) - my preferred peak caller.

#### Bug in the original SICER?

I would appreciate it if anyone can send me the output islands with FDR from running the original SICER on the example data they provide. I seem to remember finding a bug in the original, but cannot be bothered to get SICER running again. (It wasn't critical or anything just giving slightly inaccurate counts for islands.)
